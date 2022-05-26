function preal = stokes_slp_pressure_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi)
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
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
% Output:
%       pressure, velocity as a 1xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
preal = zeros(1, Ntar);

for n=1:Nsrc
    
    for m=1:Ntar 
        
        for jpx = -1:1 %Go through a layer of boxes in the x direction
            
            for jpy = -1:1 %Go through a layer of boxes in the y direction
                
                %Compute periodic source point
                xsrc_p = xsrc(n) + jpx*Lx;
                ysrc_p = ysrc(n) + jpy*Ly;
                
                r1 = -(xsrc_p - xtar(m));
                r2 = -(ysrc_p - ytar(m));
                r = sqrt(r1^2 + r2^2);
                               
                if abs(r) < 1e-13
                    continue
                else
                    ptmp =  stokeslet_pressure_real_sum(r1,r2,f1(n),f2(n),xi);
                    preal(m) = preal(m) + ptmp;
                end
            end
            
        end
    end
end

preal = preal / (2*pi);

end

function pressure = stokeslet_pressure_real_sum(r1, r2, f1, f2, xi)

rdotf = r1*f1 + r2*f2;
r = sqrt(r1^2 + r2^2);

pressure = rdotf*exp(-xi*xi*r^2) / r^2;

end
