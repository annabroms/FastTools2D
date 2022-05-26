function omega_k = stokes_slp_vorticity_kspace_ds(xsrc,ysrc,xtar,ytar,f1,f2,Lx,Ly,xi,kinf)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the Fourier space sum of the vorticity directly. 
% Evaluates between -kinf and kinf, and does not use FFTs or fast Gaussian 
% gridding.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points
%       xsrc, x component of target points
%       ytar, y component of target points
%       f1, x component of density function
%       f2, y component of density function
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
%       kinf, truncation length (integer)
% Output:
%       omega_k, vorticity as 1xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
omega_k = zeros(1, Ntar);
 
%Source points
for n = 1:Nsrc  
    
    omega_k_tmp = zeros(1,Ntar);
    
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
                
                tmp = vorticity_slp_k_sum(k1, k2, f1(n), f2(n), xi);
                omega_k_tmp(m) = omega_k_tmp(m) + tmp*exp(1i*(k1*r1+k2*r2));
 
            end
        end
    end
    
    omega_k = omega_k + omega_k_tmp;
    
end

omega_k = real(omega_k)/(Lx*Ly);

end

function omega_k = vorticity_slp_k_sum(k1, k2, f1, f2, xi)
fdotkperp = f1*k2 - f2*k1;

kSq = k1^2+k2^2;

omega_k = -1i*fdotkperp*(1/kSq+1/(4*xi*xi))*exp(-kSq/(4*xi*xi)); 
end
