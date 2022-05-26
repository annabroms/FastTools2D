function omega_r = stokes_slp_vorticity_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the vorticity directly. Considers one
% periodic replicate box in each direction.
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
% Output:
%       omega_r, vorticity as a 1xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
omega_r = zeros(1, Ntar);

for n=1:Nsrc
    
    for m=1:Ntar 
        
        for jpx = -1:1 %Go through a layer of boxes in the x direction
            
            for jpy = -1:1 %Go through a layer of boxes in the y direction
                
                %Compute periodic source point
                xsrc_p = xsrc(n) - jpx*Lx;
                ysrc_p = ysrc(n) - jpy*Ly;
                
                r1 = xtar(m) - xsrc_p;
                r2 = ytar(m) - ysrc_p;
                r = sqrt(r1^2 + r2^2);
                               
                if abs(r) < 1e-13
                    continue
                else
                    omega_tmp = vorticity_slp_real_sum(r1,r2,f1(n),f2(n),xi);
                    omega_r(m) = omega_r(m) + omega_tmp;
                end
            end
            
        end
    end
end

end

function omega_real = vorticity_slp_real_sum(r1, r2, f1, f2, xi)
fdotrperp = f1*r2 - f2*r1;

rSq = r1^2 + r2^2;
xi2 = xi*xi;

omega_real = exp(-xi2*rSq)*(1/rSq-xi2)*fdotrperp/(2*pi);

end
