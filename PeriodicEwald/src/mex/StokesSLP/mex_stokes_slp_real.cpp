#include "mex.h"
#include <math.h>
#include <omp.h>
#include <pmmintrin.h>

#include "mm_mxmalloc.h"
#include "ewald_tools.h"

#define pi 3.1415926535897932385

/*Comments by Fredrik Fryklund denoted by FF*/
/*Comments by Sara Pålsson denoted by SP*/


#define _mm_shuf2_pd(__A) (static_cast<__m128d>(__builtin_ia32_shufpd (static_cast<__v2df>(__A), static_cast<__v2df>(__A), 1)))

//Various series expansion coefficients for the computation of E1

    static const double __attribute__ ((aligned (16))) PQ0[8] = {
        -17.70313744792479226930481672752649,-2.03222164752550948918496942496859,
        -9.56230625623594754358691716333851,24.24775264986217493401454703416675,
        -0.99999982131078080094255255971802,34.82867640350680460414878325536847,
        -0.00000000083503088648841284312372,11.56228921223568129050818242831156
    };

    static const double __attribute__ ((aligned (16))) PQ1[22] = {
        -1185.45720315201027667L,-0.776491285282330997549L,
        -14751.4895786128450662L,1229.20784182403048905L,
        -54844.4587226402067411L,18455.4124737722049515L,
        -86273.1567711649528784L,86722.3403467334749201L,
        -66598.2652345418633509L,180329.498380501819718L,
        -27182.6254466733970467L,192104.047790227984431L,
        -6046.8250112711035463L,113057.05869159631492L,
        -724.581482791462469795L,38129.5594484818471461L,
        -43.3058660811817946037L,7417.37624454689546708L,
        -0.999999999999998811143L,809.193214954550328455L,
        -0.121013190657725568138e-18L,45.3058660811801465927L
    };
    static const double Y = 0.66373538970947265625L;

    static const double __attribute__ ((aligned (16))) PQ2[12] = {
        -0.000111507792921197858394L,-0.528611029520217142048e-6L,
        -0.00399167106081113256961L,0.000131049900798434683324L,
        -0.0368031736257943745142L,0.00427347600017103698101L,
        -0.245088216639761496153L,0.056770677104207528384L,
        0.0320913665303559189999L,0.37091387659397013215L,
        0.0865197248079397976498L,1L,
    };

/*
void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x, int nside_y,
        int* particle_offsets_src,int* box_offsets_src,int* nsources_in_box,
        int* particle_offsets_tar,int* box_offsets_tar,int* ntargets_in_box);

inline double expint(double x);*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /*List used for translating sources. FF */
    static const int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    static const int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};

    /* (x,y)-coordinates of source points SP */
    double *psrc = mxGetPr(prhs[0]);
    /* (x,y)-coordinates of target points SP */
    double *ptar = mxGetPr(prhs[1]);
    /* nbr of source points SP */
    int Nsrc = mxGetN(prhs[0]);
    /* nbr of target points */
    int Ntar = mxGetN(prhs[1]);
    /* (x,y)-components of force at points psrc SP */
    double *f = mxGetPr(prhs[2]);
    /*The Ewald parameter xi. FF*/
    double xi = mxGetScalar(prhs[3]);
    /* Number of boxes along the periodic box SP */
    /* Allow different number of boxes in x and y directions LB*/
    int nside_x = static_cast<int>(mxGetScalar(prhs[4]));
    int nside_y = static_cast<int>(mxGetScalar(prhs[5]));
    /*Length L of the periodic domain. FF*/
    /*Allow different lengths in x and y directions LB*/
    double len_x = mxGetScalar(prhs[6]);
    double len_y = mxGetScalar(prhs[7]);
    /*#Boxes contained in L^2, given by nside_x*nside_y. FF*/
    int num_boxes = nside_x*nside_y;

    /*offsets refer to values that are added to a base pointer in order to access a single element in a sequential list of elements. FF */
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


    /*16-byte aligned arrays. Required to be compatible with SSE commands. FF*/
    double* psrc_a = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* ptar_a = static_cast<double*>(_mm_mxMalloc (2*Ntar*sizeof(double), 16));
    double* fs = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* us = static_cast<double*>(_mm_mxCalloc (2*Ntar,sizeof(double),16));

/*Write from psrc/ptar to psrc_a/ptar_a, f to fs and u to us, which are compatible with SSE commands. FF + SP */
    for(int j = 0;j<Nsrc;j++) {
        psrc_a[2*j] = psrc[2*particle_offsets_src[j]];
        psrc_a[2*j+1] = psrc[2*particle_offsets_src[j]+1];

        fs[2*j] = f[2*particle_offsets_src[j]+1];
        fs[2*j+1] = f[2*particle_offsets_src[j]];

    }

    for(int j = 0;j<Ntar;j++) {
        ptar_a[2*j] = ptar[2*particle_offsets_tar[j]];
        ptar_a[2*j+1] = ptar[2*particle_offsets_tar[j]+1];
    }

    double self = -1.288607832450766155 - log(xi);

    /*Cut off radius squared. nside such that all points in a box are within the distance sqrt(cutoffsq) of each other. FF*/
    double cutoffsq = len_x*len_y/nside_x/nside_y;
   // mexPrintf("xi : %f, rc : %e, numboxes = %d\n", xi, sqrt(cutoffsq), num_boxes);
//     mexPrintf("R1 : %f, R2 : %f\n",cutoffsq,36.0/xi/xi);

#pragma omp parallel for
/*Loop through boxes*/
/*------------------------------------------------------------------------
 *Commands starting with _m or __m are SSE commands. SSE = Streaming SIMD Extensions.   *It is a set of CPU instructions dedicated to applications like scientific computation.
SIMD is an acronym itself: Single Instruction, Multiple Data.
 *------------------------------------------------------------------------
 FF*/

for(int current_box = 0;current_box<num_boxes;current_box++) {
    
    //mexPrintf("Box: %d, sources :%d, targets: %d\n", current_box, nsources_in_box[current_box], ntargets_in_box[current_box]);
    if(ntargets_in_box[current_box] == 0)
        continue;

    /*Temporary pointers to the particles of current box.*/
    int tidx = box_offsets_tar[current_box];
    int sidx = box_offsets_src[current_box];

    /*Compute the box self-interactions.*/
    for(int j=tidx;j<tidx+ntargets_in_box[current_box];j++) {
        for(int k=sidx;k<sidx+nsources_in_box[current_box];k++) {
            /*r2 = abs(x-y)*abs(x-y). FF*/
            double r2;
            /*loads x =(x1,y1) and y = (x2,y2), then xy = x-y = (x1-x2,y1-y2). FF */
            __m128d xy = _mm_sub_pd(_mm_load_pd(&ptar_a[2*j]),_mm_load_pd(&psrc_a[2*k]));
            /*tt= ((x1-x2)^2,(y1-y2)^2). FF*/
            __m128d tt = _mm_mul_pd(xy,xy);
            /*Store (x1-x2)^2 + (y1-y2)^2 as r2. FF*/
            _mm_storeh_pd(&r2,_mm_hadd_pd(tt,tt));
            
            if(fabs(r2) < 1e-15) {   
                us[2*j] += self*fs[2*k];
                us[2*j+1] += self*fs[2*k+1];
            }else{
              /*Compute E1, the exponential integral, by Padé approximant. FF*/
                double x = r2*xi*xi,y,e=exp(-x);
                if(x >= 16) {
                    /*Far field: 8-term Padé-approximant in 1/x*/
                    double recip = 1/x;
                    __m128d yy = _mm_load_pd(PQ0);
                    __m128d xx = _mm_set1_pd(recip);
                    double __attribute__ ((aligned (16))) res[2];

                    for(int k = 2;k<8;k+=2)
                        yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ0[k]));

                    _mm_store_pd(res,yy);

                    y = e*(recip+res[0]/(res[1] + x));
                }else if(x > 1) {
                    /*Mid field: 22-term Padé-approximant in 1/x*/
                    double recip = 1/x;
                    __m128d yy = _mm_load_pd(PQ1);
                    __m128d xx = _mm_set1_pd(recip);
                    double __attribute__ ((aligned (16))) res[2];

                    for(int k = 2;k<22;k+=2)
                        yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ1[k]));

                    _mm_store_pd(res,yy);

                    y = e*(recip+res[0]/(res[1] + x));

                }else{
                    /*Near field: 12-term Padé-approximant in in x*/
                    __m128d yy = _mm_load_pd(PQ2);
                    __m128d xx = _mm_set1_pd(x);
                    double __attribute__ ((aligned (16))) res[2];

                    for(int k = 2;k<12;k+=2)
                        yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ2[k]));

                    _mm_store_pd(res,yy);
                    y = res[0]/res[1] + x - log(x) - Y;
                }

                /*In order to avoid condfusion: let y  be the result from 
                 * computing E1 above, call this result E1. Then let y in 
                 * xy correspond to a point neighboring x, that is: 
                 * xy = x-y = (x1-x2,y1-y2), where y =(x2,y2) != E1. FF*/

                /* t1 = (0.5*E1,0.5*E1). FF*/
                __m128d t1 = _mm_set1_pd(0.5*y);
                /* t2 = (1/(r2))*e^(-xi*xi*r2) * (x1-x2,y1-y2). FF*/
                __m128d t2 = _mm_mul_pd(xy,_mm_set1_pd(e/r2));
                /*ff = (fk2,fk1). See the mapping of indices above. FF*/
                __m128d ff = _mm_load_pd(&fs[2*k]);
                /*tmp = (fk2*(x1-x2),fk1*(y1-y2)). FF*/
                __m128d tmp = _mm_mul_pd(xy,ff);
                /*tmp = (1/(r2))*e^(-xi*xi*r2) *((x1-x2)*(fk2*(x1-x2)-fk1*(y1-y2)),(y1-y2)*(fk2*(x1-x2)-fk1*(y1-y2))). FF*/
                tmp = _mm_mul_pd(t2,_mm_hsub_pd(tmp,tmp));
                /*tmp = (0.5*fk2*E1, 0.5*fk1*E1) + (1/(r2))*e^(-xi*xi*r2) * (-(x1-x2)*(fk2*(x1-x2)-fk1*(y1-y2)),(y1-y2)*(fk2*(x1-x2)-fk1*(y1-y2))). FF*/
                tmp = _mm_addsub_pd(_mm_mul_pd(ff,t1),tmp);
                /*us(j) = us(j) + (0.5*fk2*E1 + (1/(r2))*e^(-xi*xi*r2)*(x1-x2)*(-fk2*(x1-x2)+fk1*(y1-y2)), 0.5*fk1*E1+(1/(r2))*e^(-xi*xi*r2)*(y1-y2)*(fk2*(x1-x2)-fk1*(y1-y2))). FF*/
                _mm_store_pd(&us[2*j],_mm_add_pd(_mm_load_pd(&us[2*j]),tmp));
                
                //mexPrintf("r2 = %3.3g, x = %3.3g, y=%3.3g, us[%d] = %3.3g\n", r2, x, y, 2*j, us[2*j]);
            }
        }
    }

    /*Compute interactions from the nearest neighbors. On a uniform periodic grid, each box has eight neighbors. FF*/
    for(int j=0;j<8;j++) {
        /*Algorithm to loop over neighboring boxes .FF*/
        int per_source_x = current_box%nside_x+ilist_x[j];
        int per_source_y = current_box/nside_x+ilist_y[j];
        
        int t_x = (per_source_x+nside_x)%nside_x;
        int t_y = (per_source_y+nside_y)%nside_y;

        /*The number of the source nearest neighbor box.*/
        int source_box = t_y*nside_x + t_x;
       // mexPrintf("Box %d, looking in box %d\n", current_box, source_box);
        
        /*If the box is nonempty evaluate contribution. Look up "Verlet list" for more details. FF*/
        if(nsources_in_box[source_box] > 0) {
            /*z-offset of the source box corrected for periodicity.*/
            double zoff_re = (len_x*(per_source_x-t_x))/nside_x;
            double zoff_im = (len_y*(per_source_y-t_y))/nside_y;
            /*Loop over targets in the current box to evaluate contributions from source box. FF + SP */
            
            
            for(int k=0;k<ntargets_in_box[current_box];k++) {
                /*tu = (0,0). FF*/
                __m128d tu = _mm_setzero_pd();
                __m128d pk = _mm_sub_pd(_mm_load_pd(&ptar_a[2*(tidx+k)]),_mm_setr_pd(zoff_re,zoff_im));
                int idx = box_offsets_src[source_box];
                for(int l=0;l<nsources_in_box[source_box];l++,idx++) {
                    double r2;
                    __m128d xy = _mm_sub_pd(_mm_load_pd(&psrc_a[2*idx]),pk);
                    __m128d tt = _mm_mul_pd(xy,xy);
                    _mm_storeh_pd(&r2,_mm_hadd_pd(tt,tt));
                    /* Check if the points are within the cutoff. FF*/
                    if(r2< cutoffsq) {
                        /*Compute E1, the exponential integral*/
                        double x = r2*xi*xi,y,e=exp(-x);
                        if(x >= 16) {
                            /*Far field: 8-term Padé-approximant in 1/x*/
                            double recip = 1/x;
                            __m128d yy = _mm_load_pd(PQ0);
                            __m128d xx = _mm_set1_pd(recip);
                            double __attribute__ ((aligned (16))) res[2];

                            for(int k = 2;k<8;k+=2)
                                yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ0[k]));

                            _mm_store_pd(res,yy);

                            y = e*(recip+res[0]/(res[1] + x));
                        }else if(x > 1) {
                            /*Mid field: 22-term Padé-approximant in 1/x*/
                            double recip = 1/x;
                            __m128d yy = _mm_load_pd(PQ1);
                            __m128d xx = _mm_set1_pd(recip);
                            double __attribute__ ((aligned (16))) res[2];

                            for(int k = 2;k<22;k+=2)
                                yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ1[k]));

                            _mm_store_pd(res,yy);

                            y = e*(recip+res[0]/(res[1] + x));

                        }else{
                            /*Near field: 12-term Padé-approximant in in x*/
                            __m128d yy = _mm_load_pd(PQ2);
                            __m128d xx = _mm_set1_pd(x);
                            double __attribute__ ((aligned (16))) res[2];

                            for(int k = 2;k<12;k+=2)
                                yy = _mm_add_pd(_mm_mul_pd(xx,yy),_mm_load_pd(&PQ2[k]));

                            _mm_store_pd(res,yy);
                            y = res[0]/res[1] + x - log(x) - Y;
                        }

                        /*In order to avoid condfusion: let y  
                         * be the result from computing E1 above, 
                         * call this result E1. Then let y in xy correspond
                         * to a point neighboring x, that is: 
                         * xy = x-y = (x1-x2,y1-y2), where y =(x2,y2) != E1. FF*/

                        /* t1 = (0.5*E1,0.5*E1). FF*/
                        __m128d t1 = _mm_set1_pd(0.5*y);
                        /* t2 = ((1/(r2))*e^(-xi*xi*r2), (1/(r2))*e^(-xi*xi*r2)). FF*/
                        __m128d t2 = _mm_set1_pd(e/r2);
                        /*ff = (fidx f). FF*/
                        __m128d ff = _mm_load_pd(&fs[2*idx]);
                        __m128d tmp = _mm_mul_pd(xy,ff);
                        tmp = _mm_mul_pd(xy,_mm_hsub_pd(tmp,tmp));
                        tmp = _mm_mul_pd(tmp,t2);
                        tmp = _mm_addsub_pd(_mm_mul_pd(ff,t1),tmp);
                        tu = _mm_add_pd(tu,tmp);
                    }

                }
                _mm_store_pd(&us[2*(tidx+k)],_mm_add_pd(_mm_load_pd(&us[2*(tidx+k)]),tu));
                
               // mexPrintf("us[%d] = %3.3g\n", 2*(tidx+k), us[2*(tidx+k)]);
            }
        }
    }
}

/*Clean up. FF*/
_mm_mxFree(psrc_a);
_mm_mxFree(ptar_a);
_mm_mxFree(fs);

plhs[0] = mxCreateDoubleMatrix(2, Ntar, mxREAL);
double* u = mxGetPr(plhs[0]);
/*Write computed values to output vector u. Note the mapping of indices: u[i] = us[i+1]. FF*/
for(int j = 0;j<Ntar;j++) {
    u[2*particle_offsets_tar[j]] = us[2*j+1] / (4*pi);
    u[2*particle_offsets_tar[j]+1] = us[2*j] / (4*pi);
}

/*Clean up. FF*/
_mm_mxFree(us);

delete particle_offsets_src;
delete box_offsets_src;
delete nsources_in_box;
delete particle_offsets_tar;
delete box_offsets_tar;
delete ntargets_in_box;
}
