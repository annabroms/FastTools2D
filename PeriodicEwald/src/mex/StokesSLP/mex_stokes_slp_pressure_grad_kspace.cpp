#include <math.h>
#include <omp.h>
#include <string.h>
#include "ewald_tools.h"
#define pi 3.1415926535897932385

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[2]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetN(prhs[2]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    
    //Source and target points.
    double* psrc = mxGetPr(prhs[0]);
    int Nsrc = mxGetN(prhs[0]);
    double* ptar = mxGetPr(prhs[1]);
    int Ntar = mxGetN(prhs[1]);
    
    //Source strengths
    double* f = mxGetPr(prhs[2]);
    
    //Ewald parameter xi
    double xi = mxGetScalar(prhs[3]);
    
    //Splitting parameter eta
    double eta = mxGetScalar(prhs[4]);    
    
    //Number of grid intervals in each direction
    int Mx = static_cast<int>(mxGetScalar(prhs[5]));
    int My = static_cast<int>(mxGetScalar(prhs[6]));
    
    //The length of the domain
    double Lx = mxGetScalar(prhs[7]);
    double Ly = mxGetScalar(prhs[8]);
    
    //Width of the Gaussian bell curves
    double w = mxGetScalar(prhs[9]);
    
    //Number of support points for Gaussians
    int P = static_cast<int>(mxGetScalar(prhs[10]));
    
    //Grid spacing, here we assume hx = hy = h
    double h = Lx/Mx;
    
    //---------------------------------------------------------------------
    //Step 1 : Spreading to the grid
    //---------------------------------------------------------------------
    //Here we spread f1 and f2 to the grids H1 and H2
    mxArray *fft2rhs[2],*fft2lhs[2];
    fft2rhs[0] = mxCreateDoubleMatrix(My, Mx, mxREAL);    
    fft2rhs[1] = mxCreateDoubleMatrix(My, Mx, mxREAL);    
    
    double* H1 = mxGetPr(fft2rhs[0]);
    double* H2 = mxGetPr(fft2rhs[1]);
    
    //e1 contains some Gaussian gridding information that will be used in 
    //the gathering step later
    double* e1 = new double[P+1];    
    Spread(H1, H2, e1, psrc, f, Nsrc, Lx, Ly, xi, w, eta, P, Mx, My, h);
    
    //---------------------------------------------------------------------
    //Step 2 : Frequency space filter
    //---------------------------------------------------------------------
    //We apply the appropriate k-space filter associated with the
    //Stokeslet. This involves Fast Fourier Transforms.
    
    //Call Matlab's routines to compute the 2D FFT. Other choices,
    //such as fftw could possibly be better, mainly because of the poor
    //complex data structure Matlab uses which we now have to deal with.
    mexCallMATLAB(1,&fft2lhs[0],1,&fft2rhs[0],"fft2");
    mexCallMATLAB(1,&fft2lhs[1],1,&fft2rhs[1],"fft2");     

    //The output of the FFT is complex. Get pointers to the real and
    //imaginary parts of Hhat1 and Hhat2.
    double* Hhat1_re = mxGetPr(fft2lhs[0]);
    double* Hhat1_im = mxGetPi(fft2lhs[0]);
    
    double* Hhat2_re = mxGetPr(fft2lhs[1]);
    double* Hhat2_im = mxGetPi(fft2lhs[1]);
    
    mwSize cs = Mx*My;
    
    //We cannot assume that both the real and imaginary part of the
    //Fourier transforms are non-zero.
    if(Hhat1_im == NULL) {
        Hhat1_im = (double*) mxCalloc(cs,sizeof(double));        
        mxSetPi(fft2lhs[0],Hhat1_im);
    }
    if(Hhat2_im == NULL) {        
        Hhat2_im = (double*) mxCalloc(cs,sizeof(double));        
        mxSetPi(fft2lhs[1],Hhat2_im);
    }
       
    //Apply filter in the frequency domain. This is a completely
    //parallel operation. The FFT gives the frequency components in
    //non-sequential order, so we have to split the loops. One could
    //possibly use fast gaussian gridding here too, to get rid of the
    //exponentials, but as this loop accounts for a few percent of the
    //total runtime, this hardly seems worth the extra work.
#pragma omp parallel for
    for(int j = 0;j<Mx;j++) {
        int ptr = j*My;
        double k1;
        if(j <= Mx/2)
            k1 = 2.0*pi/Lx*j;
        else
            k1 = 2.0*pi/Lx*(j-Mx);
        for(int k = 0;k<=My/2;k++,ptr++) {
            double k2 = 2.0*pi/Ly*k;
            double Ksq = k1*k1+k2*k2;
            double e = exp(-Ksq*(1-eta)/(4*xi*xi));
            
            //Hhat1 contains density component 1 convolved with Gaussians
            //Hhat2 contains density component 2 convolved with Gaussians
            double kdotq_re = k1 * Hhat1_re[ptr] + k2 * Hhat2_re[ptr];
            double kdotq_im = k1 * Hhat1_im[ptr] + k2 * Hhat2_im[ptr];
            
            Hhat1_re[ptr] = -k1*kdotq_re * e / Ksq;
            Hhat1_im[ptr] = -k1*kdotq_im * e / Ksq; 
            
            Hhat2_re[ptr] = -k2*kdotq_re * e / Ksq;
            Hhat2_im[ptr] = -k2*kdotq_im * e / Ksq;   
        }
        for(int k = 0;k<My/2-1;k++,ptr++) {
            double k2 = 2.0*pi/Ly*(k-My/2+1);
            double Ksq = k1*k1+k2*k2;
            double e = exp(-Ksq*(1-eta)/(4*xi*xi));
            
            //Hhat1 contains density component 1 convolved with Gaussians
            //Hhat2 contains density component 2 convolved with Gaussians                
            double kdotq_re = k1 * Hhat1_re[ptr] + k2 * Hhat2_re[ptr];
            double kdotq_im = k1 * Hhat1_im[ptr] + k2 * Hhat2_im[ptr];
            
            Hhat1_re[ptr] = -k1*kdotq_re * e / Ksq;
            Hhat1_im[ptr] = -k1*kdotq_im * e / Ksq; 
            
            Hhat2_re[ptr] = -k2*kdotq_re * e / Ksq;
            Hhat2_im[ptr] = -k2*kdotq_im * e / Ksq; 
        }
    }
    
    //Remove the zero frequency term.
    Hhat1_re[0] = 0;
    Hhat1_im[0] = 0;
    Hhat2_re[0] = 0;
    Hhat2_im[0] = 0;
    
    
    //Get rid of the old H1 and H2 arrays. They are no longer needed.
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    
    //Do the inverse 2D FFT. We use Matlab's inbuilt functions again.
    mexCallMATLAB(1,&fft2rhs[0],1,&fft2lhs[0],"ifft2");
    mexCallMATLAB(1,&fft2rhs[1],1,&fft2lhs[1],"ifft2");
    
    //The pointer to the real part. We don't need the imaginary part.
    double* Ht1 = mxGetPr(fft2rhs[0]);
    double* Ht2 = mxGetPr(fft2rhs[1]);
    
    if(Ht1 == NULL)
        Ht1 = new double[Mx*My];    
    if(Ht2 == NULL)
        Ht2 = new double[Mx*My];
    
    

    //Get rid of the Hhat1 and Hhat2 arrays. They are no longer needed.
    mxDestroyArray(fft2lhs[0]);
    mxDestroyArray(fft2lhs[1]);
    
    //---------------------------------------------------------------------
    //Step 3 : Evaluating the pressure
    //---------------------------------------------------------------------
    //Here we compute the output pressure gradient as a convolution of the Ht1
    //and Ht2 functions with a properly scaled gaussian. This is basically
    //gaussian blur, and we again use fast gaussian gridding. This
    //procedure is fully parallel.
    
    //Create the output matrix.
    plhs[0] = mxCreateDoubleMatrix(2, Ntar, mxREAL);
    double* pressure_grad = mxGetPr(plhs[0]);
    
    Gather(Ht1, 2, 1, e1, ptar, pressure_grad, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    Gather(Ht2, 2, 2, e1, ptar, pressure_grad, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    
    //Clean up
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    delete e1;
}