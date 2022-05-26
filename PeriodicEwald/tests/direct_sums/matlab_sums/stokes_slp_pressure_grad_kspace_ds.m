function pgradk = stokes_slp_pressure_grad_kspace_ds(xsrc,ysrc,xtar,ytar,f1,f2,Lx,Ly,xi,kinf)
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
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
%       kinf, truncation length (integer)
% Output:
%       pk, pressure as 1xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
pgradk = zeros(2, Ntar);
 
%Source points
for n = 1:Nsrc  
    
    pk_tmp = zeros(2,Ntar);
    
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
                
                tmp = stokeslet_pressure_k_sum(k1, k2, f1(n), f2(n), xi);
                pk_tmp(:,m) = pk_tmp(:,m) + tmp*exp(1i*(k1*r1+k2*r2));
 
            end
        end
    end
    
    pgradk = pgradk + pk_tmp;
    
end

pgradk = real(pgradk)/(Lx*Ly);

end

function pk = stokeslet_pressure_k_sum(k1, k2, f1, f2, xi)

    kdotf = k1*f1 + k2*f2;
    k = sqrt(k1^2 + k2^2);
    
    pk = -[k1;k2]*kdotf*exp(-k^2/(4*xi*xi))/k^2;
end


