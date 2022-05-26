#include "mex.h"
#include <math.h>
#include <omp.h>
#include <pmmintrin.h>

#include "mm_mxmalloc.h"
#include "ewald_tools.h"

#define pi 3.1415926535897932385

/*Comments by Fredrik Fryklund denoted by FF*/
/*Comments by Sara Pålsson denoted by SP*/

void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc,
        int ntar, int nside_x, int nside_y,
        int* particle_offsets_src,int* box_offsets_src,int* nsources_in_box,
        int* particle_offsets_tar,int* box_offsets_tar,int* ntargets_in_box);

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
    
    /*Boxes along the periodic box. FF*/
    int nside_x = static_cast<int>(mxGetScalar(prhs[5]));
    int nside_y = static_cast<int>(mxGetScalar(prhs[6]));
    int num_boxes = nside_x*nside_y;
    
    /*Length L of the periodic domain. FF*/
    double len_x = mxGetScalar(prhs[7]);
    double len_y = mxGetScalar(prhs[8]);
    
    /*offsets refer to values that are added to a base pointer in order to 
     * access a single element in a sequential list of elements. FF */
    /* For source particles: */
    int* particle_offsets_src = new int[Nsrc];
    int* box_offsets_src = new int[num_boxes+1];
    int* nsources_in_box = new int[num_boxes];
    /* For target particles: */
    int* particle_offsets_tar = new int[Ntar];
    int* box_offsets_tar = new int[num_boxes+1];
    int* ntargets_in_box = new int[num_boxes];
    
    /*Assigns particles to boxes on the current grid. FF*/
    Assign(psrc,ptar,len_x,len_y,Nsrc,Ntar,nside_x,nside_y,
            particle_offsets_src,box_offsets_src,nsources_in_box,
            particle_offsets_tar,box_offsets_tar,ntargets_in_box);
    
    double* psrc_a = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* ptar_a = static_cast<double*>(_mm_mxMalloc (2*Ntar*sizeof(double), 16));
    double* f_a = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* n_a  = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* omega = static_cast<double*>(_mm_mxCalloc (Ntar,sizeof(double),16));

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
        
    /*Cut off radius squared. nside such that all points in a box are 
     * within the distance sqrt(cutoffsq) of each other. FF*/
    double cutoffsq = len_x*len_y/nside_x/nside_y;
    double xi2 = xi*xi;
    
#pragma omp parallel for
    for(int current_box = 0;current_box<num_boxes;current_box++) {
        if(ntargets_in_box[current_box] == 0)
            continue;
        
        /*Temporary pointers to the particles of current box.*/
        int tidx = box_offsets_tar[current_box];
        int sidx = box_offsets_src[current_box];
        
        /*Compute the box self-interactions.*/
        for(int j=tidx;j<tidx+ntargets_in_box[current_box];j++) {
            for(int k=sidx;k<sidx+nsources_in_box[current_box];k++) {
                
                double r1 = ptar_a[2*j] - psrc_a[2*k];
                double r2 = ptar_a[2*j+1] - psrc_a[2*k+1];
                
                double rSq = r1*r1+r2*r2;
                
                if(rSq < 1e-13) {
                    continue;
                }
          
                double e2 = exp(-xi2*rSq);
                double f1 = f_a[2*k];
                double f2 = f_a[2*k+1];
                double n1 = n_a[2*k];
                double n2 = n_a[2*k+1];
                
                double rdotn = r1*n1 + r2*n2;
                double rdotf = r1*f1 + r2*f2;

                double rdot = r1*(n2*rdotf+f2*rdotn) - r2*(n1*rdotf+f1*rdotn);
                double ndot = n1*(r2*rdotf) - n2*(r1*rdotf);
                double fdot = f1*(r2*rdotn) - f2*(r1*rdotn);
                
                omega[j] += e2*((1+xi2*rSq)*rdot/(rSq*rSq) + xi2*xi2*(ndot+fdot));
            }
        }
        
        /*Compute interactions from the nearest neighbors. 
         * On a uniform periodic grid, each box has eight neighbors. FF*/
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
                        
                        double r1 = ptar_a[2*(tidx+k)]- (psrc_a[2*idx]+zoff_re);
                        double r2 = ptar_a[2*(tidx+k)+1] -(psrc_a[2*idx+1]+zoff_im);
                        
                        double rSq = r1*r1+r2*r2;
                        
                        if(rSq < cutoffsq) {
                            
                            double e2 = exp(-xi2*rSq);
                            double f1 = f_a[2*idx];
                            double f2 = f_a[2*idx+1];
                            double n1 = n_a[2*idx];
                            double n2 = n_a[2*idx+1];

                            double rdotn = r1*n1 + r2*n2;
                            double rdotf = r1*f1 + r2*f2;

                            double rdot = r1*(n2*rdotf+f2*rdotn) - r2*(n1*rdotf+f1*rdotn);
                            double ndot = n1*(r2*rdotf) - n2*(r1*rdotf);
                            double fdot = f1*(r2*rdotn) - f2*(r1*rdotn);

                            omega[tidx+k] += e2*((1+xi2*rSq)*rdot/(rSq*rSq) + xi2*xi2*(ndot+fdot));
                        }
                    }
                }
            }
        }
    }
    
    /*Clean up. FF*/
    _mm_mxFree(psrc_a);
    _mm_mxFree(ptar_a);
    _mm_mxFree(f_a);
    _mm_mxFree(n_a);
    
    plhs[0] = mxCreateDoubleMatrix(1, Ntar, mxREAL);
    double* omega_out = mxGetPr(plhs[0]);

    for(int j = 0;j<Ntar;j++) {
        omega_out[particle_offsets_tar[j]] = omega[j]/pi;
    }
    
    /*Clean up. FF*/
    _mm_mxFree(omega);
    
    delete particle_offsets_src;
    delete box_offsets_src;
    delete nsources_in_box;
    delete particle_offsets_tar;
    delete box_offsets_tar;
    delete ntargets_in_box;
}
