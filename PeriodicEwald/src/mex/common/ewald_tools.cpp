#include "ewald_tools.h"

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double* psrc, double* ptar, double Lx, double Ly, int nsrc,
        int ntar, int nside_x, int nside_y, int* particle_offsets_src,
        int* box_offsets_src,int* nsources_in_box, int* particle_offsets_tar,
        int* box_offsets_tar,int* ntargets_in_box){
    
    //The total number of boxes.
    int number_of_boxes = nside_x*nside_y;
    int j;
    int* in_box_src = new int[nsrc];
    int* in_box_tar = new int[ntar];
    int* particlenum_in_box_src = new int[number_of_boxes+1];
    int* particlenum_in_box_tar = new int[number_of_boxes+1];
    
    //Clear the nparticles_in_box arrays.
    for(j = 0;j<number_of_boxes;j++){
        nsources_in_box[j] = 0;
        ntargets_in_box[j] = 0;
    }
    
    //Assign the sources to boxes.
    for(j = 0;j<nsrc;j++) {
        int box_x = (int)floor(nside_x*(psrc[2*j]/Lx+0.5));
        int box_y = (int)floor(nside_y*(psrc[2*j+1]/Ly+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_src[j] =  box_y*nside_x + box_x;
        nsources_in_box[in_box_src[j]]++;
    }
    
    //Assign the targets to boxes.
    for(j = 0;j<ntar;j++) {
        int box_x = (int)floor(nside_x*(ptar[2*j]/Lx+0.5));
        int box_y = (int)floor(nside_y*(ptar[2*j+1]/Ly+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_tar[j] =  box_y*nside_x + box_x;
        ntargets_in_box[in_box_tar[j]]++;
    }
    
    //box_offsets is the offsets of the boxes in particle_offsets to be
    //computed below. particlenum_in_box is a temporary array.
    box_offsets_src[0] = particlenum_in_box_src[0] = 0;
    box_offsets_tar[0] = particlenum_in_box_tar[0] = 0;
    int isrc = 0;
    int itar = 0;
    
    for(j=1;j<number_of_boxes+1;j++) {
        isrc += nsources_in_box[j-1];
        
        particlenum_in_box_src[j] = isrc;
        
        box_offsets_src[j] = isrc;
        
        itar += ntargets_in_box[j-1];
        particlenum_in_box_tar[j] = itar;
        
        box_offsets_tar[j] = itar;
    }
    
    //The offsets of the particles in each box in the z and q arrays.
    for(j=0;j<nsrc;j++)
        particle_offsets_src[particlenum_in_box_src[in_box_src[j]]++] = j;
    for(j=0;j<ntar;j++)
        particle_offsets_tar[particlenum_in_box_tar[in_box_tar[j]]++] = j;
    
    delete in_box_src;
    delete particlenum_in_box_src;
    delete in_box_tar;
    delete particlenum_in_box_tar;
}

/*------------------------------------------------------------------------
 *This function finds the node to begin the Gaussian blur
 *------------------------------------------------------------------------
 */
void FindClosestNode(double x, double y, double Lx, double Ly, double h, int P,
        int* mx, int* my, double* px, double* py){
    
    double TOL = 1e-13;
    
    *px = x - h*floor(x/h);
    *py = y - h*floor(y/h);
    
    double Lhalf_x = static_cast<double>(Lx/2.0);
    double Lhalf_y = static_cast<double>(Ly/2.0);
    
    //(mx,my) is the center grid point in H1 and H2.
    *mx = static_cast<int>(floor((x+Lhalf_x)/h) - P/2);
    *my = static_cast<int>(floor((y+Lhalf_y)/h) - P/2);
    
    // correct for cases where target is very close to a grid node
    
    if (fabs(*px) < TOL || fabs(*px - h) < TOL)
        *px = 0;
    
    if (fabs(*py) < TOL || fabs(*py - h) < TOL)
        *py = 0;
    
    if (*px == 0 && remainder((x+Lhalf_x)/h - P/2, 1) < 0)
        (*mx)++;
    
    if (*py == 0 && remainder((y+Lhalf_y)/h - P/2, 1) < 0)
        (*my)++;
}

/*------------------------------------------------------------------------
 *This function speads a vector field to a uniform grid
 *------------------------------------------------------------------------
 */
void Spread(double* H1, double* H2, double* e1, double* psrc, double* f,
        int Nsrc, double Lx, double Ly, double xi, double w,
        double eta, int P, int Mx, int My, double h){
    
    //This is the precomputable part of the fast Gaussian gridding.
    double tmp = -2*xi*xi/eta*h*h;
    for(int j = -P/2;j<=P/2;j++)
        e1[j+P/2] = exp(tmp*j*j);
    
    //Spreading the sources to the grid is not a completely parallel
    //operation. We use the simple approach of locking the column of the
    //matrix we are working on currently. This might not be optimal but it
    //is simple to implement.
    omp_lock_t* locks = new omp_lock_t[Mx];
    for(int j = 0;j<Mx;j++)
        omp_init_lock(&locks[j]);
    
    //We use OpenMP for simple parallelization.
#pragma omp parallel for
    for(int k = 0;k<Nsrc;k++) {
        
        //The Gaussian bells are translation invariant. We exploit this
        //fact to avoid blow-up of the terms and the numerical instability
        //that follows. (px,py) is the center of the bell with the original
        //grid-alignment but close to the origin.
        
        double xsrc = psrc[2*k];
        double ysrc = psrc[2*k+1];
        
        int mx, my;
        double px, py;
        
        FindClosestNode(xsrc, ysrc, Lx, Ly, h, P, &mx, &my, &px, &py);
        
        //Some auxillary quantities for the fast Gaussian gridding.
        double tmp = -2*xi*xi/eta;
        double ex = exp(tmp*(px*px+py*py + 2*w*px));
        double e4y = exp(2*tmp*w*py);
        double e3x = exp(-2*tmp*h*px);
        double e3y = exp(-2*tmp*h*py);
        
        //We add the Gaussians column by column, and lock the one we are
        //working on to avoid race conditions.
        for(int x = 0;x<P+1;x++) {
            double ey = ex*e4y*e1[x];
            int xidx = ((x+mx+Mx)%Mx)*My;
            
            omp_set_lock(&locks[(x+mx+Mx)%Mx]);
            if(my >= 0 && my < My-P-1) {
                int idx = my+xidx;
                for(int y = 0;y<P+1;y++,idx++) {
                    double tmp = ey*e1[y];
                    
                    H1[idx] += tmp*f[2*k];
                    H2[idx] += tmp*f[2*k+1];
                    ey *= e3y;
                }
            }
            else{
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    int idx = ((y+my+My)%My)+xidx;
                    
                    H1[idx] += tmp*f[2*k];
                    H2[idx] += tmp*f[2*k+1];
                    ey *= e3y;
                }
            }
            omp_unset_lock(&locks[(x+mx+Mx)%Mx]);
            ex *= e3x;
        }
    }
    
    //Spreading is the only part of the k-space sum that needs locks, so
    //get rid of them.
    for(int j = 0;j<Mx;j++)
        omp_destroy_lock(&locks[j]);
    
    delete locks;
}

/*------------------------------------------------------------------------
 *This function performs the evaluation step, gathering the data at the 
 *target points
 *------------------------------------------------------------------------
 */
void Gather(double* H, int total_components, int component_number, 
        double* e1, double* ptar, double* output,
        int Ntar, double Lx, double Ly, double xi, double w,
        double eta, int P, int Mx, int My, double h){
    
#pragma omp parallel for
    for(int k = 0;k<Ntar;k++) {
        
        int target_index = total_components*k + component_number-1;
        //The Gaussian bells are translation invariant. We exploit this
        //fact to avoid blow-up of the terms and the numerical instability
        //that follows. (px,py) is the center of the bell with the original
        //grid-alignment but close to the origin.
        
        double xtar = ptar[2*k];
        double ytar = ptar[2*k+1];
        
        int mx, my;
        double px, py;
        
        FindClosestNode(xtar, ytar, Lx, Ly, h, P, &mx, &my, &px, &py);
        
        double tmp = -2*xi*xi/eta;
        double ex = exp(tmp*(px*px+py*py + 2*w*px));
        double e4y = exp(2*tmp*w*py);
        double e3x = exp(-2*tmp*h*px);
        double e3y = exp(-2*tmp*h*py);
        
        //If there is no wrap-around due to periodicity for this gaussian,
        //we use a faster loop. We go column by column, but as the
        //matrices Ht1 and Ht2 are only read from we have no need for
        //locks and such.
        if(mx >= 0 && my >= 0 && mx < Mx-P-1 && my < My-P-1) {
            int idx = mx*My+my;
            for(int x = 0;x<P+1;x++) {
                double ey = ex*e4y*e1[x];
                
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    output[target_index] += tmp*H[idx];
                    idx++;
                    ey *= e3y;
                }
                ex *= e3x;
                idx += My-P-1;
            }
        }
        else{
            for(int x = 0;x<P+1;x++) {
                double ey = ex*e4y*e1[x];
                int xidx = ((x+mx+Mx)%Mx)*My;
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    int idx = ((y+my+My)%My)+xidx;
                    output[target_index] += tmp*H[idx];
                    ey *= e3y;
                }
                ex *= e3x;
            }
        }
        
        tmp = 4*xi*xi/eta;
        tmp = tmp*tmp*h*h/pi;
        output[target_index] *= tmp / (4*pi);
    }
}

/*------------------------------------------------------------------------
 *This function extracts the real and imaginary parts of an FFT
 *------------------------------------------------------------------------
 */
void ExtractRealIm(mxArray *fftvector, double *Hhat_re, double *Hhat_im,
        int Mx, int My)
{
    
    mexPrintf("1\n");
    
    //The output of the FFT is complex. Get pointers to the real and
    //imaginary parts of Hhat1 and Hhat2.
    Hhat_re = mxGetPr(fftvector);
    Hhat_im = mxGetPi(fftvector);
    
    mwSize cs = Mx*My;
    
    mexPrintf("2\n");
    
    //We cannot assume that both the real and imaginary part of the
    //Fourier transforms are non-zero.
    if(Hhat_im == NULL) {
        Hhat_im = (double*) mxCalloc(cs,sizeof(double));
        mxSetPi(fftvector,Hhat_im);
    }
    if(Hhat_re == NULL) {
        Hhat_re = (double*) mxCalloc(cs,sizeof(double));
        mxSetPr(fftvector,Hhat_re);
    }
}
