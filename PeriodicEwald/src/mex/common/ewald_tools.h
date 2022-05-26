#ifndef EWALD_TOOLS  
#define EWALD_TOOLS

#include <math.h>
#include <string.h>
#include <omp.h>
#include "mex.h"

#define pi 3.1415926535897932385

void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x, int nside_y, int* particle_offsets_src,
        int* box_offsets_src,int* nsources_in_box, int* particle_offsets_tar,
        int* box_offsets_tar,int* ntargets_in_box);

void FindClosestNode(double x, double y, double Lx, double Ly, double h, int P, 
			int* mx, int* my, double* px, double* py);

void Spread(double* H1, double* H2, double* e1, double* psrc, double* f, 
                int Nsrc, double Lx, double Ly, double xi, double w,
                double eta, int P, int Mx, int My, double h);

void Gather(double* H, int total_components, int component_number, 
        double* e1, double* ptar, double* output,
        int Ntar, double Lx, double Ly, double xi, double w,
        double eta, int P, int Mx, int My, double h);
        
void ExtractRealIm(mxArray *fftvector, double *Hhat_re, double *Hhat_im, 
        int Mx, int My);
#endif
