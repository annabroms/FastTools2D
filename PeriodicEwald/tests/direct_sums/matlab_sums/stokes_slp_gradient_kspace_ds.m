function ugrad_k = stokes_slp_gradient_kspace_ds(xsrc,ysrc,xtar,ytar,f1,f2,b1,b2,Lx,Ly,xi,kinf)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the Fourier space sum of the Stokeslet pressure directly. 
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
%       b1, x component of target direction vector
%       b2, y component of target direction vector
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
%       kinf, truncation length (integer)
% Output:
%       ugrad_k, velocity gradient as 1xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
ugrad_k = zeros(2, Ntar);
 
%Source points
for n = 1:Nsrc  
    
    ugrad_k_tmp = zeros(2,Ntar);
    
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
                
                tmp = stokeslet_gradient_k_sum(k1, k2, f1(n), f2(n), b1(m), b2(m), xi);
                ugrad_k_tmp(:,m) = ugrad_k_tmp(:,m) + tmp*exp(1i*(k1*r1+k2*r2));
 
            end
        end
    end
    
    ugrad_k = ugrad_k + ugrad_k_tmp;
    
end

ugrad_k = real(ugrad_k)/(Lx*Ly);

end

function ugrad_k = stokeslet_gradient_k_sum(k1, k2, f1, f2, b1, b2, xi)

    ugrad_k = [0;0];
    
    kdotf = k1*f1 + k2*f2;
    kdotb = k1*b1 + k2*b2;
    
    k = sqrt(k1^2 + k2^2);
    
    ugrad_k(1) = 1i*(f1*kdotb/k^2 - k1*kdotb*kdotf/k^4)*(1 + k^2/(4*xi*xi))...
                    *exp(-k^2/(4*xi*xi));
    ugrad_k(2) = 1i*(f2*kdotb/k^2 - k2*kdotb*kdotf/k^4)*(1 + k^2/(4*xi*xi))...
                    *exp(-k^2/(4*xi*xi));            
end


