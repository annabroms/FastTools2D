#include "ewald_tools.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 12)
        mexErrMsgTxt("Incorrect number of input parameters");
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[4]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetM(prhs[5]) != 2)
        mexErrMsgTxt("n must be a 2xn matrix.");
    if(mxGetN(prhs[4]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    if(mxGetN(prhs[5]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and n must be the same size.");
    
    //The points.
    double* psrc = mxGetPr(prhs[0]);
    int Nsrc = mxGetN(prhs[0]);
    double* ptar = mxGetPr(prhs[1]);
    int Ntar = mxGetN(prhs[1]);
    
    //The Ewald parameter xi
    double xi = mxGetScalar(prhs[2]);
    //The splitting parameter eta
    double eta = mxGetScalar(prhs[3]);
    //The Stresslet vectors.
    double* f = mxGetPr(prhs[4]);
    double* n = mxGetPr(prhs[5]);
    
    int Mx = static_cast<int>(mxGetScalar(prhs[6]));
    int My = static_cast<int>(mxGetScalar(prhs[7]));
    
    //The length of the domain
    double Lx = mxGetScalar(prhs[8]);
    double Ly = mxGetScalar(prhs[9]);
    
    //The width of the Gaussian bell curves. (unnecessary?)
    double w = mxGetScalar(prhs[10]);
    
    //The width of the Gaussian bell curves on the grid.
    int P = static_cast<int>(mxGetScalar(prhs[11]));
    
    //Grid spacing, assuming hx = hy = h
    double h = Lx/Mx;
    
    double mu = 1.0;
    
    //---------------------------------------------------------------------
    //Step 1 : Spreading to the grid
    //---------------------------------------------------------------------
    //We begin by spreading the sources onto the grid. This basically means
    //that we superposition properly scaled Gaussian bells, one for
    //each source. We use fast Gaussian gridding to reduce the number of
    //exps we need to evaluate.
    
    //The function H on the grid. We need to have these as Matlab arrays
    //since we call Matlab's in-built fft2 to compute the 2D FFT.
    mxArray *fft2rhs[4],*fft2lhs[4];
    fft2rhs[0] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    fft2rhs[1] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    fft2rhs[2] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    fft2rhs[3] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    
    double* H1 = mxGetPr(fft2rhs[0]);    
    double* H2 = mxGetPr(fft2rhs[1]);    
    double* H3 = mxGetPr(fft2rhs[2]);    
    double* H4 = mxGetPr(fft2rhs[3]);
    
    //Have to multiply components of f and n before speading
    double* v1 = new double[2*Nsrc];
    double* v2 = new double[2*Nsrc];
    for (int i = 0; i < Nsrc; i++)
    {
        v1[2*i] = f[2*i]*n[2*i];          //f1 * n1
        v1[2*i + 1] = f[2*i+1]*n[2*i];    //f2 * n1
        v2[2*i] = f[2*i]*n[2*i+1];        //f1 * n2
        v2[2*i + 1] = f[2*i+1]*n[2*i+1];  //n2 * n2
    }
    
    //This is the precomputable part of the fast Gaussian gridding.
    double* e1 = new double[P+1];
    Spread(H1, H2, e1, psrc, v1, Nsrc, Lx, Ly, xi, w, eta, P, Mx, My, h);
    Spread(H3, H4, e1, psrc, v2, Nsrc, Lx, Ly, xi, w, eta, P, Mx, My, h);
    
    //---------------------------------------------------------------------
    //Step 2 : Frequency space filter
    //---------------------------------------------------------------------
    //We apply the appropriate k-space filter associated with the
    //Stresslet. This involves Fast Fourier Transforms.
    
    //Call Matlab's routines to compute the 2D FFT. Other choices,
    //such as fftw could possibly be better, mainly because of the poor
    //complex data structure Matlab uses which we now have to deal with.
    mexCallMATLAB(1,&fft2lhs[0],1,&fft2rhs[0],"fft2");
    mexCallMATLAB(1,&fft2lhs[1],1,&fft2rhs[1],"fft2");
    mexCallMATLAB(1,&fft2lhs[2],1,&fft2rhs[2],"fft2");
    mexCallMATLAB(1,&fft2lhs[3],1,&fft2rhs[3],"fft2");
    
    //The output of the FFT is complex. Get pointers to the real and
    //imaginary parts of the Hhats.
    double* Hhat1_re = mxGetPr(fft2lhs[0]);
    double* Hhat1_im = mxGetPi(fft2lhs[0]);
    
    double* Hhat2_re = mxGetPr(fft2lhs[1]);
    double* Hhat2_im = mxGetPi(fft2lhs[1]);
    
    double* Hhat3_re = mxGetPr(fft2lhs[2]);
    double* Hhat3_im = mxGetPi(fft2lhs[2]);
    
    double* Hhat4_re = mxGetPr(fft2lhs[3]);
    double* Hhat4_im = mxGetPi(fft2lhs[3]);
    
    //We cannot assume that the imaginary parts of the
    //Fourier transforms are non-zero. Let Matlab take care of the
    //memory management.
    mwSize cs = Mx*My;
    
    if(Hhat1_im == NULL) {
        Hhat1_im = (double*) mxCalloc(cs,sizeof(double));
        mxSetPi(fft2lhs[0],Hhat1_im);
    }
    if(Hhat2_im == NULL) {
        Hhat2_im = (double*) mxCalloc(cs,sizeof(double));
        mxSetPi(fft2lhs[1],Hhat2_im);
    }
    
    if(Hhat3_im == NULL) {
        Hhat3_im = new double[Mx*My];
       mxSetPi(fft2lhs[2],Hhat3_im);
    }
    
    if(Hhat4_im == NULL) {
        Hhat4_im = new double[Mx*My];
        mxSetPi(fft2lhs[3],Hhat4_im);
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
            
            double f1n1_re = Hhat1_re[ptr];
            double f1n1_im = Hhat1_im[ptr];
            double f1n2_re = Hhat2_re[ptr];
            double f1n2_im = Hhat2_im[ptr];
            double f2n1_re = Hhat3_re[ptr];
            double f2n1_im = Hhat3_im[ptr];
            double f2n2_re = Hhat4_re[ptr];
            double f2n2_im = Hhat4_im[ptr];
            
            //j = 1, l = 1
            Hhat1_re[ptr] = e*((k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k1*(f1n1_re+f2n2_re)+k1*(k1*f1n1_re+k2*f1n2_re)
                +k1*(k1*f1n1_re+k2*f2n1_re)+k1*(k1*f1n1_re+k2*f1n2_re)+k1*(k1*f1n1_re+k2*f2n1_re)
                -4*k1*k1*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq));
            Hhat1_im[ptr] = e*((k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_im+f2n2_im)+k1*(k1*f1n1_im+k2*f1n2_im)
                +k1*(k1*f1n1_im+k2*f2n1_im)+k1*(k1*f1n1_im+k2*f1n2_im)+k1*(k1*f1n1_im+k2*f2n1_im)
                -4*k1*k1*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq));
            
            //j = 2, l = 1
            Hhat2_re[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k1*(f1n1_re+f2n2_re)+k2*(k1*f1n1_re+k2*f1n2_re)
                +k2*(k1*f1n1_re+k2*f2n1_re)+k1*(k1*f2n1_re+k2*f2n2_re)+k1*(k1*f1n2_re+k2*f2n2_re)
                -4*k2*k1*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq);
            Hhat2_im[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k1*(f1n1_im+f2n2_im)+k2*(k1*f1n1_im+k2*f1n2_im)
                +k2*(k1*f1n1_im+k2*f2n1_im)+k1*(k1*f2n1_im+k2*f2n2_im)+k1*(k1*f1n2_im+k2*f2n2_im)
                -4*k2*k1*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq);
            
            //j = 1, l = 2
            Hhat3_re[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_re+f2n2_re)+k1*(k1*f2n1_re+k2*f2n2_re)
                +k1*(k1*f1n2_re+k2*f2n2_re)+k2*(k1*f1n1_re+k2*f1n2_re)+k2*(k1*f1n1_re+k2*f2n1_re)
                -4*k1*k2*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq);
            Hhat3_im[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_im+f2n2_im)+k1*(k1*f2n1_im+k2*f2n2_im)
                +k1*(k1*f1n2_im+k2*f2n2_im)+k2*(k1*f1n1_im+k2*f1n2_im)+k2*(k1*f1n1_im+k2*f2n1_im)
                -4*k1*k2*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq);
            
            //j = 2, l = 2
            Hhat4_re[ptr] = e*((k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k2*(f1n1_re+f2n2_re)+k2*(k1*f2n1_re+k2*f2n2_re)
                +k2*(k1*f1n2_re+k2*f2n2_re)+k2*(k1*f2n1_re+k2*f2n2_re)+k2*(k1*f1n2_re+k2*f2n2_re)
                -4*k2*k2*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq));
            Hhat4_im[ptr] = e*((k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k2*(f1n1_im+f2n2_im)+k2*(k1*f2n1_im+k2*f2n2_im)
                +k2*(k1*f1n2_im+k2*f2n2_im)+k2*(k1*f2n1_im+k2*f2n2_im)+k2*(k1*f1n2_im+k2*f2n2_im)
                -4*k2*k2*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq));
            
        }
        for(int k = 0;k<My/2-1;k++,ptr++) {
            double k2 = 2.0*pi/Ly*(k-My/2+1);
            double Ksq = k1*k1+k2*k2;
            
            double e = exp(-Ksq*(1-eta)/(4*xi*xi));
            
            double f1n1_re = Hhat1_re[ptr];
            double f1n1_im = Hhat1_im[ptr];
            double f1n2_re = Hhat2_re[ptr];
            double f1n2_im = Hhat2_im[ptr];
            double f2n1_re = Hhat3_re[ptr];
            double f2n1_im = Hhat3_im[ptr];
            double f2n2_re = Hhat4_re[ptr];
            double f2n2_im = Hhat4_im[ptr];
            
            //j = 1, l = 1
            Hhat1_re[ptr] = e*((k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k1*(f1n1_re+f2n2_re)+k1*(k1*f1n1_re+k2*f1n2_re)
                +k1*(k1*f1n1_re+k2*f2n1_re)+k1*(k1*f1n1_re+k2*f1n2_re)+k1*(k1*f1n1_re+k2*f2n1_re)
                -4*k1*k1*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq));
            Hhat1_im[ptr] = e*((k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_im+f2n2_im)+k1*(k1*f1n1_im+k2*f1n2_im)
                +k1*(k1*f1n1_im+k2*f2n1_im)+k1*(k1*f1n1_im+k2*f1n2_im)+k1*(k1*f1n1_im+k2*f2n1_im)
                -4*k1*k1*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq));
            
            //j = 2, l = 1
            Hhat2_re[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k1*(f1n1_re+f2n2_re)+k2*(k1*f1n1_re+k2*f1n2_re)
                +k2*(k1*f1n1_re+k2*f2n1_re)+k1*(k1*f2n1_re+k2*f2n2_re)+k1*(k1*f1n2_re+k2*f2n2_re)
                -4*k2*k1*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq);
            Hhat2_im[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k1*(f1n1_im+f2n2_im)+k2*(k1*f1n1_im+k2*f1n2_im)
                +k2*(k1*f1n1_im+k2*f2n1_im)+k1*(k1*f2n1_im+k2*f2n2_im)+k1*(k1*f1n2_im+k2*f2n2_im)
                -4*k2*k1*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq);
            
            //j = 1, l = 2
            Hhat3_re[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_re+f2n2_re)+k1*(k1*f2n1_re+k2*f2n2_re)
                +k1*(k1*f1n2_re+k2*f2n2_re)+k2*(k1*f1n1_re+k2*f1n2_re)+k2*(k1*f1n1_re+k2*f2n1_re)
                -4*k1*k2*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq);
            Hhat3_im[ptr] = -e*mu*(1/Ksq+0.25/(xi*xi))*(2*k1*k2*(f1n1_im+f2n2_im)+k1*(k1*f2n1_im+k2*f2n2_im)
                +k1*(k1*f1n2_im+k2*f2n2_im)+k2*(k1*f1n1_im+k2*f1n2_im)+k2*(k1*f1n1_im+k2*f2n1_im)
                -4*k1*k2*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq);

            //j = 2, l = 2
            Hhat4_re[ptr] = e*((k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k2*(f1n1_re+f2n2_re)+k2*(k1*f2n1_re+k2*f2n2_re)
                +k2*(k1*f1n2_re+k2*f2n2_re)+k2*(k1*f2n1_re+k2*f2n2_re)+k2*(k1*f1n2_re+k2*f2n2_re)
                -4*k2*k2*(k1*k1*f1n1_re+k1*k2*f1n2_re+k2*k1*f2n1_re+k2*k2*f2n2_re)/Ksq));
            Hhat4_im[ptr] = e*((k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq
                -mu*(1/Ksq+0.25/(xi*xi))*(2*k2*k2*(f1n1_im+f2n2_im)+k2*(k1*f2n1_im+k2*f2n2_im)
                +k2*(k1*f1n2_im+k2*f2n2_im)+k2*(k1*f2n1_im+k2*f2n2_im)+k2*(k1*f1n2_im+k2*f2n2_im)
                -4*k2*k2*(k1*k1*f1n1_im+k1*k2*f1n2_im+k2*k1*f2n1_im+k2*k2*f2n2_im)/Ksq));
            
        }
    }
    
    //Remove the zero frequency term.
    Hhat1_re[0] = 0;
    Hhat1_im[0] = 0;
    Hhat2_re[0] = 0;
    Hhat2_im[0] = 0;
    Hhat3_re[0] = 0;
    Hhat3_im[0] = 0;
    Hhat4_re[0] = 0;
    Hhat4_im[0] = 0;
    
    //Get rid of the old H arrays. They are no longer needed.
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    mxDestroyArray(fft2rhs[2]);
    mxDestroyArray(fft2rhs[3]);
    
    //Do the inverse 2D FFTs. We use Matlab's inbuilt functions again.
    mexCallMATLAB(1,&fft2rhs[0],1,&fft2lhs[0],"ifft2");
    mexCallMATLAB(1,&fft2rhs[1],1,&fft2lhs[1],"ifft2");
    mexCallMATLAB(1,&fft2rhs[2],1,&fft2lhs[2],"ifft2");
    mexCallMATLAB(1,&fft2rhs[3],1,&fft2lhs[3],"ifft2");
    
    //The pointer to the real part. We don't need the imaginary part.
    double* Ht1 = mxGetPr(fft2rhs[0]);
    double* Ht2 = mxGetPr(fft2rhs[1]);
    double* Ht3 = mxGetPr(fft2rhs[2]);
    double* Ht4 = mxGetPr(fft2rhs[3]);
    
    //In case the inverse FFTs are purely imaginary. Not likely to happen, 
    //but the program would crash without guarding for this.
    if(Ht1 == NULL) {
        Ht1 = new double[Mx*My];
        memset(Ht1,0,Mx*My*sizeof(double));
    }
    if(Ht2 == NULL) {
        Ht2 = new double[Mx*My];
        memset(Ht2,0,Mx*My*sizeof(double));
    }
    if(Ht3 == NULL) {
        Ht3 = new double[Mx*My];
        memset(Ht3,0,Mx*My*sizeof(double));
    }
    if(Ht4 == NULL) {
        Ht4 = new double[Mx*My];
        memset(Ht3,0,Mx*My*sizeof(double));
    }
    
    //Get rid of the Hhat arrays. They are no longer needed.
    mxDestroyArray(fft2lhs[0]);
    mxDestroyArray(fft2lhs[1]);
    mxDestroyArray(fft2lhs[2]);
    mxDestroyArray(fft2lhs[3]);
    
    //---------------------------------------------------------------------
    //Step 3 : Evaluating the velocity
    //---------------------------------------------------------------------
    //Here we compute the output velocity as a convolution of the Ht1
    //and Ht2 functions with a properly scaled gaussian. This is basically
    //gaussian blur, and we again use fast gaussian gridding. This
    //procedure is fully parallel.
    
    //Create the output matrix.
    plhs[0] = mxCreateDoubleMatrix(4, Ntar, mxREAL);
    double* Tk = mxGetPr(plhs[0]);
    
    Gather(Ht1, 4, 1, e1, ptar, Tk, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    Gather(Ht2, 4, 2, e1, ptar, Tk, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    Gather(Ht3, 4, 3, e1, ptar, Tk, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    Gather(Ht4, 4, 4, e1, ptar, Tk, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
    
    //Clean up
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    mxDestroyArray(fft2rhs[2]);
    mxDestroyArray(fft2rhs[3]);
    
    delete e1;
}
