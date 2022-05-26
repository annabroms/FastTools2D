function sigmak = stokes_dlp_stress_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, b1, b2, Lx, Ly, xi, kinf)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the Stokeslet directly. Considers one
% periodic replicate box in each direction.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points
%       xsrc, x component of target points
%       ytar, y component of target points
%       f1, x component of density function
%       f2, y component of density function
%       n1, x component of normal vector at source points
%       n2, y component of normal vector at source points
%       b1, x component of target direction vector
%       b2, y component of target direction vector
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
%       kinf, truncation length (integer)
% Output:
%       sigmak, stress as an 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Nsrc = length(xsrc);
Ntar = size(xtar,1);
sigmak = zeros(2,Ntar);
 
%Source points
for n = 1:Nsrc  
    
    sigma_tmp = zeros(2,Ntar);
    
    %Target points   
    for m = 1:Ntar      
        
        r1 = xtar(m) - xsrc(n);
        r2 = ytar(m) - ysrc(n);
        
        % Loop through Fourier modes in x and y directions
        for jk1 = -kinf+1:kinf 
            
            k1 = 2*pi/Lx*jk1;
            
            for jk2 = -kinf+1:kinf
                k2 = 2*pi/Ly*jk2;
                
                if k1 == 0 && k2 == 0
                    continue
                end
                
                tmp = stokes_dlp_stress_k_sum(k1, k2, ...
                            n1(n), n2(n), f1(n), f2(n), b1(m), b2(m), xi);
                sigma_tmp(:,m) = sigma_tmp(:,m) + tmp*exp(1i*(k1*r1+k2*r2));
                
            end
        end
    end
    
    sigmak = sigmak + sigma_tmp;
    
end

sigmak = real(sigmak) / (Lx*Ly);

% add on zero mode
sigmak = sigmak + sum((n1.*f1 + n2.*f2))/(2*Lx*Ly);

end

function sigmak = stokes_dlp_stress_k_sum(k1, k2, n1, n2, f1, f2, b1, b2, xi)
sigmak = [0;0];

kdotf = k1*f1 + k2*f2;
kdotn = k1*n1 + k2*n2;
fdotn = f1*n1 + f2*n2;
kdotb = k1*b1 + k2*b2;
ndotb = n1*b1 + n2*b2;
fdotb = f1*b1 + f2*b2;

kSq = k1^2 + k2^2;
xi2 = xi*xi;
e2 = exp(-kSq/(4*xi2));
mu = 1;

sigmak(1) = e2*(b1*kdotf*kdotn/kSq-mu*(1/kSq+1/(4*xi2))*(2*k1*kdotb*fdotn + ...
    k1*kdotn*fdotb+k1*kdotf*ndotb+f1*kdotb*kdotn+n1*kdotb*kdotf - ...
    4*k1*kdotb*kdotf*kdotn/kSq));

sigmak(2) = e2*(b2*kdotf*kdotn/kSq-mu*(1/kSq+1/(4*xi2))*(2*k2*kdotb*fdotn + ...
    k2*kdotn*fdotb+k2*kdotf*ndotb+f2*kdotb*kdotn+n2*kdotb*kdotf - ...
    4*k2*kdotb*kdotf*kdotn/kSq));

end



