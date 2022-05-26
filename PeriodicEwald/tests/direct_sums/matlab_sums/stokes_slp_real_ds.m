function ureal = stokes_slp_real_ds(xsrc, ysrc, xtar, ytar,...
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
%       ureal, velocity as an 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
ureal1 = zeros(1, Ntar);
ureal2 = zeros(1, Ntar);

for n=1:Nsrc
    
    for m=1:Ntar 
        
        for jpx = -1:1 %Go through a layer of boxes in the x direction
            
            for jpy = -1:1 %Go through a layer of boxes in the y direction
                
                %Compute periodic source point
                xsrc_p = xsrc(n) + jpx*Lx;
                ysrc_p = ysrc(n) + jpy*Ly;
                
                r1 = (xsrc_p - xtar(m));
                r2 = (ysrc_p - ytar(m));
                r = sqrt(r1^2 + r2^2);
                               
                % if r == 0, add self interaction term
                if abs(r) < 1e-13
                    ureal1(m)  = ureal1(m) - (double(eulergamma)/2 + log(xi) + 1)*f1(n);
                    ureal2(m)  = ureal2(m) - (double(eulergamma)/2 + log(xi) + 1)*f2(n);
                else
                    [utmp1, utmp2] =  stokeslet_real_sum(r1,r2,f1(n),f2(n),xi);
                    ureal1(m) = ureal1(m) + utmp1;
                    ureal2(m) = ureal2(m) + utmp2;
                end
            end
            
        end
    end
end

ureal = [ureal1; ureal2] / (4*pi);

end

function [u1, u2] = stokeslet_real_sum(r1, r2, f1, f2, xi)

rdotf = r1*f1 + r2*f2;
r = sqrt(r1^2 + r2^2);

u1 = exp(-xi^2*r^2)*(r1 * rdotf / r^2 - f1) + f1*expint(xi^2*r^2)/2;
u2 = exp(-xi^2*r^2)*(r2 * rdotf / r^2 - f2) + f2*expint(xi^2*r^2)/2;

end
