#MEX=/usr/local/MATLAB/R2021a/bin/mex
MEX=/opt/matlab/r2021a_u3/bin/mex
FORTRAN = gfortran
FFLAGS = -c -fpic -fopenmp -Wall -O3 -fdefault-real-8 -fdefault-double-8
#MEXLIBDIR = /usr/local/MATLAB/R2021a/sys/os/glnxa64
MEXLIBDIR = /opt/matlab/r2021a_u3/sys/os/glnxa64
MEXLIB = -lgfortran -liomp5 -lirc -lsvml -limf 

OBJS=\
lfmm2drouts.o\
cfmm2dpart.o\
rfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
prini.o\
stokesSLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts prini rfmm2dpart stokesSLP mex

mex: stokesSLPfmm.F $(OBJS) 
	$(MEX) -v -O stokesSLPfmm.F $(OBJS) -L$(MEXLIBDIR) $(MEXLIB)

cfmm2dpart: cfmm2dpart.f
	$(FORTRAN) $(FFLAGS) cfmm2dpart.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

stokesSLP: stokesSLP.f
	$(FORTRAN) $(FFLAGS) stokesSLP.f

clean: 
	rm *.o
	cp *.mex* ../src
	rm *.mex*
