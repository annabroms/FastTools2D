% This is a test script to test the mex implementation of the spectral
% Ewald method against direct matlab sums. 

close all
clearvars
clc

initewald

%% Set up data

% Ewald tolerance
tol = 1e-16;
% k trunctation for direct summation
kinf = 70;

Nsrc = 20;
Ntar = 20;

Lx = 1;
Ly = 1;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Two components of normal vector
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target locations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

% Make sure the sources and targets are all inside the box. Needed for
% correct comparison between direct and Ewald Fourier sum due to the
% non-zero zero-mode
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
xtar = mod(xtar+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;
ytar = mod(ytar+Ly/2,Ly)-Ly/2;

%% Check single-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES SINGLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
[~,~, ur_ewald, uk_ewald, xi] = StokesSLP_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, f1, f2, Lx, Ly, 'verbose', 1, 'tol', tol);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
ur_direct = stokes_slp_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(ur_direct - ur_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

uk_direct = stokes_slp_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi, kinf);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(uk_direct - uk_ewald))));

%% Check double-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES DOUBLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
[~,~, ur_ewald, uk_ewald, xi] = StokesDLP_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, n1, n2, f1, f2, Lx, Ly, 'verbose', 1);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
ur_direct = stokes_dlp_real_ds(xsrc, ysrc, xtar, ytar, n1, n2,...
                        f1, f2, Lx, Ly, xi);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(ur_direct - ur_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

uk_direct = stokes_dlp_kspace_ds(xsrc, ysrc, xtar, ytar, n1, n2,...
                        f1, f2, Lx, Ly, xi, kinf);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(uk_direct - uk_ewald))));
