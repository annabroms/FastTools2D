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

Nsrc = 10;
Ntar = 10;

Lx = 1;
Ly = 2;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Two components of normal vector
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

% xtar = xsrc;
% ytar = ysrc;
xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

% %% Check single-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES SINGLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
[~, pr_ewald, pk_ewald, xi] = StokesSLP_pressure_grad_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, f1, f2, Lx, Ly, 'verbose', 1, 'tol', tol);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
pr_direct = stokes_slp_pressure_grad_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(pr_direct - pr_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

pk_direct = stokes_slp_pressure_grad_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi, kinf);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(pk_direct - pk_ewald))));

% %% Check double-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES SINGLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
[~, pr_ewald, pk_ewald, xi] = StokesDLP_pressure_grad_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, n1, n2, f1, f2,  Lx, Ly, 'verbose', 1, 'tol', tol);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
pr_direct = stokes_dlp_pressure_grad_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2,  Lx, Ly, xi);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(pr_direct - pr_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

pk_direct = stokes_dlp_pressure_grad_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, Lx, Ly, xi, kinf);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(abs(pk_direct - pk_ewald))));
