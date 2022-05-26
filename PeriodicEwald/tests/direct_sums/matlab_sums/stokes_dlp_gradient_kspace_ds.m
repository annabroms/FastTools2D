function uk = stokes_dlp_gradient_kspace_ds(xsrc, ysrc, xtar, ytar,...
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
%       uk, velocity as an 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Nsrc = length(xsrc);
Ntar = size(xtar,1);
uk1 = zeros(1, Ntar);
uk2 = zeros(1, Ntar);
 
%Source points
for n = 1:Nsrc  
    
    uk1_tmp = zeros(1,Ntar);
    uk2_tmp = zeros(1,Ntar);
    
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
                
                [tmp_x, tmp_y] = stokes_dlp_gradient_k_sum(k1, k2, ...
                            n1(n), n2(n), f1(n), f2(n), b1(m), b2(m), xi);
                uk1_tmp(m) = uk1_tmp(m) + tmp_x*exp(1i*(k1*r1+k2*r2));
                uk2_tmp(m) = uk2_tmp(m) + tmp_y*exp(1i*(k1*r1+k2*r2));
                
            end
        end
    end
    
    uk1 = uk1 + uk1_tmp;
    uk2 = uk2 + uk2_tmp;
    
end

uk = real([uk1; uk2]) / (Lx*Ly);

% Add on zero mode
% uk(1,:) = uk(1,:) + sum((f1.*n1 + f2.*n2).*xsrc) / (Lx*Ly);
% uk(2,:) = uk(2,:) + sum((f1.*n1 + f2.*n2).*ysrc) / (Lx*Ly);

end

function [uk1, uk2] = stokes_dlp_gradient_k_sum(k1, k2, n1, n2, f1, f2, b1, b2, xi)

    kdotf = k1*f1 + k2*f2;
    kdotn = k1*n1 + k2*n2;
    fdotn = f1*n1 + f2*n2;
    kdotb = k1*b1 + k2*b2;
    
    k = sqrt(k1^2 + k2^2);
    
    uk1 = -kdotb*(f1*kdotn + n1*kdotf + k1*fdotn -...
            2*k1*kdotn*kdotf/k^2)*(1+k^2/(4*xi^2))*exp(-k^2/(4*xi^2))/k^2;
    uk2 = -kdotb*(f2*kdotn + n2*kdotf + k2*fdotn -...
            2*k2*kdotn*kdotf/k^2)*(1+k^2/(4*xi^2))*exp(-k^2/(4*xi^2))/k^2;
end



