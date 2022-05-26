#include "mex.h"
#include <math.h>
#include <omp.h>
#include <pmmintrin.h>

#include "mm_mxmalloc.h"
#include "ewald_tools.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*List used for translating sources. FF */
    static const int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    static const int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};
    
    if(nrhs != 9)
        mexErrMsgTxt("Incorrect number of input parameters");
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[2]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetM(prhs[3]) != 2)
        mexErrMsgTxt("n must be a 2xn matrix.");
    if(mxGetN(prhs[2]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    if(mxGetN(prhs[3]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and n must be the same size.");
    
    /*(x,y)-coordinates of the source and target points psrc,ptar. SP*/
    double *psrc = mxGetPr(prhs[0]);
    double *ptar = mxGetPr(prhs[1]);
    /*Number of sources and targets. FF*/
    int Nsrc = mxGetN(prhs[0]);
    int Ntar = mxGetN(prhs[1]);
    /*(x,y)-components of the forces at points p. FF*/
    double *f = mxGetPr(prhs[2]);
    double *n = mxGetPr(prhs[3]);
    /*The Ewald parameter xi. FF*/
    double xi = mxGetScalar(prhs[4]);
    /*Paramter s: #boxes along the periodic box. FF*/
    int nside_x = static_cast<int>(mxGetScalar(prhs[5]));
    int nside_y = static_cast<int>(mxGetScalar(prhs[6]));
    /*Length L of the periodic domain. FF*/
    double len_x = mxGetScalar(prhs[7]);
    double len_y = mxGetScalar(prhs[8]);
    
    /*#Boxes contained in L^2, given by s*s. FF*/
    int num_boxes = nside_x*nside_y;
    
    /*offsets refer to values that are added to a base pointer
     * in order to access a single element in a sequential list
     * of elements. FF */
    int* particle_offsets_src = new int[Nsrc];
    int* box_offsets_src = new int[num_boxes+1]; //NB added +1 SP
    int* nsources_in_box = new int[num_boxes];
    
    int* particle_offsets_tar = new int[Ntar];
    int* box_offsets_tar = new int[num_boxes+1]; //NB added +1 SP
    int* ntargets_in_box = new int[num_boxes];
    
    /*Assigns particles to boxes on the current grid. FF*/
    Assign(psrc,ptar,len_x,len_y,Nsrc,Ntar,nside_x,nside_y,
            particle_offsets_src,box_offsets_src,nsources_in_box,
            particle_offsets_tar,box_offsets_tar,ntargets_in_box);
    
    /*16-byte aligned arrays. Required to be compatible with SSE commands. FF*/
    double* psrc_a = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* ptar_a = static_cast<double*>(_mm_mxMalloc (2*Ntar*sizeof(double), 16));
    double* f_a  = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* n_a  = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* Ts = static_cast<double*>(_mm_mxCalloc (4*Ntar, sizeof(double),16));
    
    /*Write from p to ps and from f to fs, which are compatible with SSE commands. FF*/
    for(int j = 0;j<Nsrc;j++) {
        
        psrc_a[2*j] = psrc[2*particle_offsets_src[j]];
        psrc_a[2*j+1] = psrc[2*particle_offsets_src[j]+1];
        
        f_a[2*j] = f[2*particle_offsets_src[j]];
        f_a[2*j+1] = f[2*particle_offsets_src[j]+1];
        
        n_a[2*j] = n[2*particle_offsets_src[j]];
        n_a[2*j+1] = n[2*particle_offsets_src[j]+1];
        
    }
    
    for(int j = 0;j<Ntar;j++) {
        ptar_a[2*j] = ptar[2*particle_offsets_tar[j]];
        ptar_a[2*j+1] = ptar[2*particle_offsets_tar[j]+1];
    }
    /*Cut off radius squared. nside such that all points in a
     * box are within the distance sqrt(cutoffsq) of each other. FF*/
    double cutoffsq = len_x*len_y/nside_x/nside_y;
    double xi2 = xi*xi;

    /*Loop through boxes*/
#pragma omp parallel for    
    for(int current_box = 0;current_box<num_boxes;current_box++) {
        if(ntargets_in_box[current_box] == 0)
            continue;
        
        //Temporary pointers to the particles of current box.
        int tidx = box_offsets_tar[current_box];
        int sidx = box_offsets_src[current_box];
                
        //Compute the box self-interactions.
        for(int j=tidx;j<tidx+ntargets_in_box[current_box];j++) {
            
            for(int k=sidx;k<sidx+nsources_in_box[current_box];k++) {
                                
                double r1 = -(psrc_a[2*k] - ptar_a[2*j]);
                double r2 = -(psrc_a[2*k+1] - ptar_a[2*j+1]);
                
                double rSq = r1*r1+r2*r2;
                
                if(rSq == 0)
                    continue;
          
                double e2 = exp(-xi2*rSq);
                double f1 = f_a[2*k];
                double f2 = f_a[2*k+1];
                double n1 = n_a[2*k];
                double n2 = n_a[2*k+1];
                
                double rdotf = r1*f1 + r2*f2;
                double rdotn = r1*n1 + r2*n2;
                double fdotn = f1*n1 + f2*n2;
                
                //j = 1, p = 1
                Ts[4*j] += e2*(r1*r1*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(rdotf*rdotn + r1*f1*rdotn + r1*n1*rdotf)/rSq/rSq
                             +2*xi2*(2*f1*n1+fdotn-2*xi2*r1*(f1*rdotn+n1*rdotf+r1*fdotn)));
                
                //j = 2, p = 1
                Ts[4*j+1] += e2*(r1*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(r2*f1*rdotn + r2*n1*rdotf)/rSq/rSq
                             +2*xi2*(f1*n2+n1*f2-2*xi2*r1*(f2*rdotn+n2*rdotf+r2*fdotn)));
                
                //j = 1, p = 2
                Ts[4*j+2] += e2*(r1*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(r1*f2*rdotn + r1*n2*rdotf)/rSq/rSq
                             +2*xi2*(f2*n1+f1*n2-2*xi2*r2*(f1*rdotn+n1*rdotf+r1*fdotn)));
                
                //j = 2, p = 2
                Ts[4*j+3] += e2*(r2*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(rdotf*rdotn + r2*f2*rdotn + r2*n2*rdotf)/rSq/rSq
                             +2*xi2*(2*f2*n2+fdotn-2*xi2*r2*(f2*rdotn+n2*rdotf+r2*fdotn)));
            }
        }
        
        //Compute interactions from the nearest neighbors.
        for(int j=0;j<8;j++) {
            
            int per_source_x = current_box%nside_x+ilist_x[j];
            int per_source_y = current_box/nside_x+ilist_y[j];
            
            int t_x = (per_source_x+nside_x)%nside_x;
            int t_y = (per_source_y+nside_y)%nside_y;
            //The number of the source nearest neighbor box.
            int source_box = t_y*nside_x + t_x;
            
            if(nsources_in_box[source_box] > 0) {
                //z-offset of the source box corrected for periodicity.
                double zoff_re = (len_x*(per_source_x-t_x))/nside_x;
                double zoff_im = (len_y*(per_source_y-t_y))/nside_y;
                
                for(int k=0;k<ntargets_in_box[current_box];k++) {
                    
                    int idx = box_offsets_src[source_box];
                    for(int l=0;l<nsources_in_box[source_box];l++,idx++) {
                        double r1 = -(psrc_a[2*idx]-ptar_a[2*(tidx+k)]+zoff_re);
                        double r2 = -(psrc_a[2*idx+1]-ptar_a[2*(tidx+k)+1]+zoff_im);
                        
                        double rSq = r1*r1+r2*r2;
                        
                        if(rSq < cutoffsq) {

                            double e2 = exp(-xi2*rSq);
                            double f1 = f_a[2*idx];
                            double f2 = f_a[2*idx+1];
                            double n1 = n_a[2*idx];
                            double n2 = n_a[2*idx+1];
                            
                            double rdotf = r1*f1 + r2*f2;
                            double rdotn = r1*n1 + r2*n2;
                            double fdotn = f1*n1 + f2*n2;
                            
                            //j = 1, p = 1
                            Ts[4*(tidx+k)] += e2*(r1*r1*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(rdotf*rdotn + r1*f1*rdotn + r1*n1*rdotf)/rSq/rSq
                             +2*xi2*(2*f1*n1+fdotn-2*xi2*r1*(f1*rdotn+n1*rdotf+r1*fdotn)));
                            
                            //j = 2, p = 1
                            Ts[4*(tidx+k)+1] += e2*(r1*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(r2*f1*rdotn + r2*n1*rdotf)/rSq/rSq
                             +2*xi2*(f1*n2+n1*f2-2*xi2*r1*(f2*rdotn+n2*rdotf+r2*fdotn)));
                            
                            //j = 1, p = 2
                            Ts[4*(tidx+k)+2] += e2*(r1*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(r1*f2*rdotn + r1*n2*rdotf)/rSq/rSq
                             +2*xi2*(f2*n1+f1*n2-2*xi2*r2*(f1*rdotn+n1*rdotf+r1*fdotn)));
                            
                            //j = 2, p = 2
                            Ts[4*(tidx+k)+3] += e2*(r2*r2*rdotf*rdotn*(8*xi2*xi2/rSq + 16*xi2/rSq/rSq + 16/rSq/rSq/rSq)
                             -4*(1+xi2*rSq)*(rdotf*rdotn + r2*f2*rdotn + r2*n2*rdotf)/rSq/rSq
                             +2*xi2*(2*f2*n2+fdotn-2*xi2*r2*(f2*rdotn+n2*rdotf+r2*fdotn)));
                        }
                    }
                }
            }
        }
    }
    
    _mm_mxFree(psrc_a);
    _mm_mxFree(ptar_a);
    _mm_mxFree(f_a);
    _mm_mxFree(n_a);
    
    plhs[0] = mxCreateDoubleMatrix(4, Ntar, mxREAL);
    
    double* T = mxGetPr(plhs[0]);

    for(int j = 0;j<Ntar;j++) {
        T[4*particle_offsets_tar[j]] = Ts[4*j]/4/pi;
        T[4*particle_offsets_tar[j]+1] = Ts[4*j+1]/4/pi;
        T[4*particle_offsets_tar[j]+2] = Ts[4*j+2]/4/pi;
        T[4*particle_offsets_tar[j]+3] = Ts[4*j+3]/4/pi;
    }
    
    _mm_mxFree(Ts);
    
    delete particle_offsets_src;
    delete box_offsets_src;
    delete nsources_in_box;
    delete particle_offsets_tar;
    delete box_offsets_tar;
    delete ntargets_in_box;
    
}
