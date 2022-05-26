function [sigma1, sigma2, sigmar, sigmak, xi] = StokesSLP_stress_ewald_2p(xsrc, ysrc,...
            xtar, ytar, f1, f2, b1, b2, Lx, Ly, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Spectral Ewald evaluation of the doubly-periodic Stokeslet.
%
% Input:
%       xsrc, x component of source points
%       ysrc, y component of source points 
%       xtar, x component of target points
%       ytar, y component of target points 
%       f1, x component of density function
%       f2, y component of density function
%       b1, x component of target direction vector
%       b2, y component of target_direction vector
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%         'P', integer giving support points in each direction (default 24)
%         'Nb', average number of points per box (default P*log2(#pts))
%         'tol', error tolerance for truncation of sums (default 1e-16)
%         'verbose', flag to write out parameter information
% Output:
%       sigma1, x component of stress
%       sigma2, y component of stress
%       sigmar, real component of Ewald decomposition (as a 2xN matrix)
%       sigmar, Fourier component of Ewald decomposition (as a 2xN matrix)
%       xi, Ewald parameter
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

npts = length(xsrc)+length(xtar);

%% set default parameter values 
% support points in each direction
P = 24;                    
% average number of points per box for real space sum
Nb = min(P*round(log2(npts)), npts/4);
% tolerance, used to get parameters from estimates
tol = 1e-16;  
% print diagnostic information
verbose = 0;

%% read in optional input parameters
if nargin > 8
    % Go through all other input arguments and assign parameters
    jv = 1;
    while jv <= length(varargin)-1
       switch varargin{jv}
  
           case 'P'
               P = varargin{jv+1};
               
           case 'Nb'
               Nb = varargin{jv+1};
               
           case 'tol'
               tol = varargin{jv+1};
               
           case 'verbose'
               verbose = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

% TO DO: ADD CHECKS ON INPUT DATA HERE
%% Fix for matlab 2018/2019, not sure why this is necessary, but it seems 
% to work. 

% Check to see first if it's necessary, for Matlab 2017a at least it isn't.
v=ver('MATLAB');
if v.Release~="(R2017a)"
    offset = 1e-60;
    b1 = b1 - offset;
    if b1 < -1
        b1 = b1 + 2*offset;
    end
    b2 = sqrt(1 - b1.^2).*sign(b2);
    
    f1 = f1 + offset;
    f2 = f2 + offset;
end

if verbose
    fprintf("*********************************************************\n");
    fprintf("SPECTRAL EWALD FOR THE STOKES SINGLE-LAYER POTENTIAL\n\n")
    fprintf("NUMBER OF SOURCES: %d\n", length(xsrc));
    fprintf("NUMBER OF TARGETS: %d\n", length(xtar));
    fprintf("TOLERANCE: %3.3e\n", tol);
    fprintf("P: %d\n", P);
    fprintf("Points per box: %d\n", Nb);
end

%  Make sure the sources and targets are all inside the box.
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
xtar = mod(xtar+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;
ytar = mod(ytar+Ly/2,Ly)-Ly/2;

psrc = [xsrc';ysrc'];
ptar = [xtar';ytar'];
f = [f1';f2'];
b = [b1';b2'];

% compute parameters, rc, xi and kinf
[A,B] = rat(Lx/Ly);

Q = sum(sum(f.^2))+1;
a = ceil(sqrt(npts/(Nb*A*B)));

m = 0.95*sqrt(pi*P);
nside_x = a*A;
nside_y = a*B;
rc = Lx/nside_x;

xi = find_xi(Q,Lx,Ly,rc,tol);
kinfx = find_kinfb(Q,Lx,Lx,xi,tol);

Mx = min(2*kinfx,10000);

My = B * Mx;
Mx = A * Mx;

w = P*Lx/Mx/2;
eta = (2*xi*w/m)^2;

if verbose
    fprintf("\nPARAMETER INFORMATION:\n")
    fprintf("\txi: %3.3f\n", xi);
    fprintf("\trc: %3.3f\n", rc);
    fprintf("\tkinf: %d\n", kinfx);
    fprintf("\tMx: %d\n", Mx);
    fprintf("\tMy: %d\n", My);
    fprintf("\tw: %3.3f\n", w);
    fprintf("\teta: %3.3f\n", eta);
    fprintf("*********************************************************\n");
    tic
end

sigmar_tmp = mex_stokes_slp_stress_real(psrc,ptar,f,xi,nside_x,nside_y,Lx,Ly);

sigmar(1,:) = sigmar_tmp(1,:).*b1' + sigmar_tmp(3,:).*b2';
sigmar(2,:) = sigmar_tmp(2,:).*b1' + sigmar_tmp(4,:).*b2';

if verbose
    fprintf("TIME FOR REAL SUM: %3.3g s\n", toc);
    tic
end

if verbose
    fprintf("TIME FOR FOURIER SUM: %3.3g s\n", toc);
    fprintf("*********************************************************\n\n");
end

% old mex that does not work.
% don't know why row 2 and 3 are different, they should be equal. If I use
% row 2 in both places I get good result in pipe flow when running with eta
% = inf, i.e. having only the SLP.
%sigmak_tmp = mex_stokes_slp_stress_kspace_old(psrc,ptar,xi,eta,f,Mx,My,Lx,Ly,w,P);

sigmak_tmp = mex_stokes_slp_stress_kspace(psrc,ptar,xi,eta,f,Mx,My,Lx,Ly,w,P);

sigmak = zeros(2,length(xtar));
sigmak(1,:) = sigmak_tmp(1,:).*b1' + sigmak_tmp(3,:).*b2';
sigmak(2,:) = sigmak_tmp(2,:).*b1' + sigmak_tmp(4,:).*b2';

sigma = sigmar + sigmak;

% add on zero mode (from pressure)
zero_mode = sum((f1.*xsrc + f2.*ysrc)) / (2*Lx*Ly);
sigma = sigma + zero_mode;

sigma1 = sigma(1,:)';
sigma2 = sigma(2,:)';

end

%% Computing error estimates. Estimates come from Pålsson and Tornberg 2019 
% https://arxiv.org/pdf/1909.12581.pdf

% -------------------------------------------------------------------------
% Given xi and tol, finds kinfbar according to estimate
% -------------------------------------------------------------------------
function k = find_kinfb(Q,Lx,Ly,xi,tol)

% Find k using a Newton iteration (OBS kinf)
k = round(5*xi);       % Initial guess
maxit = 1e2; it = 0;

f = @(k) sqrt(4*Q*pi*max(Lx,Ly)/(Lx^3*Ly^3*k))*exp(-k^2/(4*xi^2)) - tol;
fp = @(k) -sqrt(4*Q*pi*max(Lx,Ly)/(Lx^3*Ly^3))*exp(-k^2/(4*xi^2))*...
                (0.5*k^(-1.5) + sqrt(1/k)*2*k/(4*xi^2));

%kdiff = 1;
while abs(f(k)) > tol
    if it > maxit
        disp("Couldn't find suitable k! Break")
        break;
    end
    
    % limit the step size in the negative direction
    kdiff = min(f(k)/fp(k), 0.9*k);
    k = k - kdiff;
    
    it = it + 1;
end

k = round(Lx*k/(2*pi));

end

% -------------------------------------------------------------------------
% Given rc and tol, finds xi according to estimate
% -------------------------------------------------------------------------
function x = find_xi(Q,Lx,Ly,rc,tol)

% Find xi using a Newton iteration
x = 4/rc;       % Initial guess
maxit = 1e2; it = 0;

% estimates from paper
f = @(x) sqrt(Q * pi/(4*Lx*Ly*x))*exp(-x^2*rc^2) - tol;
fp = @(x) -sqrt(Q*pi/(4*Lx*Ly))*exp(-x^2*rc^2)*(0.5*x^(-1.5)+...
                2*rc^2*x*sqrt(1/x));

%xdiff = 1;
while abs(f(x)) > tol
    if it > maxit
        disp("Couldn't find suitable xi! Break")
        break;
    end
    
    % limit the step size in the negative direction
    xdiff = min(f(x)/fp(x), 0.9*x); 
    x = x - xdiff;
    
    it = it + 1;
end

end