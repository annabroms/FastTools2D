% This is a test script to check the consistency of the Ewald method. In
% particular we look at what happens when we change the Ewald parameters,
% and what happens when we replicate the reference cell. 

close all
clearvars
clc

initewald

% Test parameters
test_self = 1;
Nsrc = 1000;
Ntar = 1000;
Lx_value = 2;
Ly_value = 1;
rng(123);

%% Set up data
% Length of periodic box
Lx = Lx_value;
Ly = Ly_value;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target locations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

if test_self
    xtar = xsrc;
    ytar = ysrc;
else
    xtar = Lx*rand(Ntar,1);
    ytar = Ly*rand(Ntar,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of single-layer potential...\n');
fprintf("*********************************************************\n\n");
%% Compute solution with Spectral Ewald
% Here we change Nb, the number of points in each box in the real space
% sum. This controls the cutoff radius in real space, and that in turn
% changes the Ewald parameter xi, and the cutoff mode in Fourier space. 
% Changing Nb should not change the result, but it will have an impact on
% the relative times to compute the real space and Fourier sums. This
% parameter can be studied to find an optimum balance.

Nb = [3, 9, 27];

omega = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    omega_tmp = StokesSLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly,...
                'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    omega(:,j) = omega_tmp;
end

%% Check that using a different number of bins doesn't affect solution

E1 = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E1(i,j) = max(abs(omega(:,i) - omega(:,j))./abs(omega(:,i)));
    end
end
fprintf('\nMaximum error from changing number of bins for SLP: %.5e\n',...
    max(max(E1)));

%% Check that replicating boxes doesn't affect solution
% Note that we have to subtract off the zero mode, because it depends on
% the source locations

omega1 = StokesSLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];

Lx = 2*Lx;

omega2 = StokesSLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);

fprintf('\nMaximum error from creating periodic replicate for SLP: %.5e\n',...
    max(abs(omega1 - omega2)./abs(omega1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of double-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Set up data
% Length of periodic box
Lx = Lx_value;
Ly = Ly_value;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Normal vector at source points
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target locations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

if test_self
    xtar = xsrc;
    ytar = ysrc;
else
    xtar = Lx*rand(Ntar,1);
    ytar = Ly*rand(Ntar,1);
end

%% Compute solution with Spectral Ewald
Nb = [3, 9, 27];

omega = zeros(Ntar, length(Nb));
xi = zeros(length(Nb),1);

for j = 1:length(Nb)
    tic
    omega_tmp = StokesDLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly,...
                'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    omega(:,j) = omega_tmp;
end

%% Check that using a different number of bins doesn't affect solution

E1 = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E1(i,j) = max(abs(omega(:,i) - omega(:,j))./abs(omega(:,i)));
    end
end
fprintf('\nMaximum error from changing number of bins for DLP: %.5e\n',...
    max(max(E1)));

%% Check that replicating boxes doesn't affect solution

omega1 = StokesDLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];
n1 = [n1; n1];
n2 = [n2; n2];

Lx = 2*Lx;

omega2 = StokesDLP_vorticity_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

fprintf('\nMaximum error from creating periodic replicate for DLP: %.5e\n',...
    max(abs(omega1 - omega2)./abs(omega1)));
