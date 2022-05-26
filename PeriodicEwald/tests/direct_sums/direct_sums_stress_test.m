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
Ly = 1;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Two components of normal vector
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

% Check single-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES SINGLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
% row 1
[~, ~, sigma1_r_ewald, sigma1_k_ewald, xi] = StokesSLP_stress_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, 'verbose', 0, 'tol', tol);

% row 2
[~, ~, sigma2_r_ewald, sigma2_k_ewald, ~] = StokesSLP_stress_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, 'verbose', 0, 'tol', tol);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
sigma1_r_direct = stokes_slp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, xi);

sigma2_r_direct = stokes_slp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, xi);
                    
fprintf('MAXIMUM ERROR IN ROW1: %.5e\n',max(max(abs(sigma1_r_direct - sigma1_r_ewald))));
fprintf('MAXIMUM ERROR IN ROW2: %.5e\n',max(max(abs(sigma2_r_direct - sigma2_r_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

sigma1_k_direct = stokes_slp_stress_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, xi, kinf);
                    
sigma2_k_direct = stokes_slp_stress_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, xi, kinf);
                    
fprintf('MAXIMUM ERROR IN ROW 1: %.5e\n',max(max(abs(sigma1_k_direct - sigma1_k_ewald))));
fprintf('MAXIMUM ERROR IN ROW 2: %.5e\n',max(max(abs(sigma2_k_direct - sigma2_k_ewald))));

%% Check double-layer potential
fprintf("*********************************************************\n");
fprintf('TESTING DIRECT SUMS FOR STOKES DOUBLE-LAYER POTENTIAL\n');
fprintf("*********************************************************\n");

% Compute solution with Spectral Ewald
% row 1
[~, ~, sigma1_r_ewald, sigma1_k_ewald, xi] = StokesDLP_stress_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, n1, n2, f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, 'verbose', 0, 'tol', tol);

% row 2
[~, ~, sigma2_r_ewald, sigma2_k_ewald, ~] = StokesDLP_stress_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, n1, n2, f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, 'verbose', 0, 'tol', tol);

fprintf('CHECKING REAL SUM...\n');
% Compute direct sums
sigma1_r_direct = stokes_dlp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, xi);

sigma2_r_direct = stokes_dlp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, xi);

fprintf('MAXIMUM ERROR IN ROW1: %.5e\n',max(max(abs(sigma1_r_direct - sigma1_r_ewald))));
fprintf('MAXIMUM ERROR IN ROW2: %.5e\n',max(max(abs(sigma2_r_direct - sigma2_r_ewald))));

fprintf("*********************************************************\n");
fprintf('CHECKING FOURIER SUM...\n');

sigma1_k_direct = stokes_dlp_stress_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, ones(Ntar,1), zeros(Ntar,1), Lx, Ly, xi, kinf);

sigma2_k_direct = stokes_dlp_stress_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, zeros(Ntar,1), ones(Ntar,1), Lx, Ly, xi, kinf);
       
fprintf('MAXIMUM ERROR IN ROW 1: %.5e\n',max(max(abs(sigma1_k_direct - sigma1_k_ewald))));
fprintf('MAXIMUM ERROR IN ROW 2: %.5e\n',max(max(abs(sigma2_k_direct - sigma2_k_ewald))));
