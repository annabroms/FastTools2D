# Fast Summation Methods for 2D Potentials

This library provides tools to compute sums of various potentials. For periodic potentials the spectral Ewald method is provided, while for non-periodic potentials the fast mulitpole method (FMM) is provided. 

Currently includes
* FMM for Stokes single- and double-layer potentials
* Spectral Ewald for doubly-periodic Stokes single- and double-layer potentials

## Compilation Instructions

### FMM

To compile the FMM, makefiles are provided. On Linux, all you should need to do is to enter the correct Matlab path and run `make`, i.e.

	cd FMM/src/StokesDLP
	make
	cd ../StokesSLP
	make

The mex files will be placed in the same directory as the source code. **This was tested by Anna using Matlab R2019a (gfortran 6.x?). Compiled files contained in each source folder**


On Apple, Matlab no longer supports gfortran so you will have to use the Intel compilers (however see possible workaround [here](https://se.mathworks.com/matlabcentral/answers/338303-how-to-set-up-mex-with-gfortran-on-mac)).

### Spectral Ewald
In principle, compiling should be as simple as running the cmake script to create the make files, and then running make. In the `src` directory, you will have to change the directories of your Matlab installation in `CMakeLists.txt`. By running everything in the `build` directory all the cmake files will be created in one place. The mex files will be placed in a new `bin` directory. The cmake script has been tested on the following architectures:

#### Ubuntu 16.04 LTS / MATLAB 2017a

	cd PeriodicEwald/src/build
	cmake ..
	make

#### MAC OS Mojave / MATLAB 2019a

Here you have to make sure that gcc and g++ are called instead of clang. I also had to manually add the include and library directories to help the compiler find the GNU scientific library.

	cd src/build
	export CPATH=/usr/local/include/
	export LIBRARY_PATH=/usr/local/lib/
	export LD_LIBRARY_PATH=/usr/local/lib/
	cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 .. 
	make

## Testing

### FMM

In the `FMM/tests` directory there are several tests:
* stokes_DLP_FMM_tests.m: compares the DLP FMM to a direct sum, and checks timings, which should scale as O(N)
* stokes_SLP_FMM_tests.m: performs the same test for the SLP
* stresslet_identity_test.m: verifies the stresslet identity for points inside and outside a circle

### Spectral Ewald
In the `tests` directory, there are several tests that can be used to verify the compilation worked correctly:
* consistency_test.m: checks that changing the Ewald parameters and enlarging the periodic box by adding replicates of the reference cell don't change the results
* direct_sums_test.m: compares the spectral Ewald implementation to matlab direct sums of the real and Fourier parts. The Matlab direct sum does not truncate in real space, and in Fourier space it does not spread the data to a uniform grid and thus does not use FFTs
* timings_test.m: checks the timings of the code for increasing numbers of source and target points. The timing should scale as O(N log N), where N is the total number of points
* stresslet_indentity_test.m: verifies the stresslet identity for points inside and outside a circle


## To do

* add more kernels
* create CMake script for FMM
* check OMP timings

## Acknowledgements

All FMM code was written by Leslie Greengard and Zydrunas Gimbutas. The spectral Ewald code was written mainly by Sara Pålsson and Rikard Ojala, with some bug fixes and code structuring by Lukas Bystricky. 

## References

Ewald decomposition of two-dimensional Stokeslet and stresslet along with error estimates:

Pålsson and Tornberg, 2019. [*An integral equation method for closely interacting surfactant-covered droplets in wall-confined Stokes flow*](https://arxiv.org/abs/1909.12581).

For details on the spectral Ewald method in general:

Lindbo and Tornberg, 2010. [*Spectrally accurate fast summation for periodic Stokes potentials*](https://www.sciencedirect.com/science/article/pii/S0021999110004730). J. Comp. Phys. 229(23). 
