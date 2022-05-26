function ureal = stokes_dlp_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, Lx, Ly, xi)
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
%       n2, y component of normal vecotr at source points
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
                
                r1 = -(xsrc_p - xtar(m));
                r2 = -(ysrc_p - ytar(m));
                r = sqrt(r1^2 + r2^2);
                               
                % if r == 0, skip
                if abs(r) < 1e-13
                    continue
                else
                    [utmp1, utmp2] =  stokeslet_real_sum(r1,r2,...
                                n1(n), n2(n), f1(n),f2(n),xi);
                    ureal1(m) = ureal1(m) + utmp1;
                    ureal2(m) = ureal2(m) + utmp2;
                end
            end
            
        end
    end
end

ureal = [ureal1; ureal2] / (4*pi);

end

function [u1, u2] = stokeslet_real_sum(r1, r2, n1, n2, f1, f2, xi)

rdotf = r1*f1 + r2*f2;
rdotn = r1*n1 + r2*n2;
ndotf = n1*f1 + n2*f2;

r = sqrt(r1^2 + r2^2);

u1 = exp(-xi^2*r^2)*(-4*r1*rdotn*rdotf/r^4*(1 + xi^2*r^2) +...
            2*xi^2*(f1*rdotn + n1*rdotf + ndotf*r1));
u2 = exp(-xi^2*r^2)*(-4*r2*rdotn*rdotf/r^4*(1 + xi^2*r^2) +...
            2*xi^2*(f2*rdotn + n2*rdotf + ndotf*r2));

end
