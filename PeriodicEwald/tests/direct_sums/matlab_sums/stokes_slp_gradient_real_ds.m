function ugrad_r = stokes_slp_gradient_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, b1, b2, Lx, Ly, xi)
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
%       b1, x component of target direction vector
%       b2, y component of target direction vector
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       xi, ewald parameter
% Output:
%       gradient of u, velocity as a 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
ugrad_r = zeros(2, Ntar);

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
                    ugrad_tmp =  stokeslet_gradient_real_sum(r1,r2,f1(n),f2(n),b1(m),b2(m),xi);
                    ugrad_r(:,m) = ugrad_r(:,m) + ugrad_tmp;
                end
            end
            
        end
    end
end

ugrad_r = ugrad_r / (4*pi);

end

function ugrad_real = stokeslet_gradient_real_sum(r1, r2, f1, f2, b1, b2, xi)

ugrad_real = [0;0];

rdotf = r1*f1 + r2*f2;
rdotb = r1*b1 + r2*b2;
bdotf = f1*b1 + f2*b2;

r = sqrt(r1^2 + r2^2);

ugrad_real(1) = exp(-xi*xi*r^2)*(2*xi*xi*rdotb*f1 + ...
        (b1*rdotf + r1*bdotf - f1*rdotb )/r^2 - 2*r1*rdotb*rdotf*(xi*xi + 1/r^2)/r^2);

ugrad_real(2) = exp(-xi*xi*r^2)*(2*xi*xi*rdotb*f2 + ...
        ( b2*rdotf + r2*bdotf - f2*rdotb)/r^2 - 2*r2*rdotb*rdotf*(xi*xi + 1/r^2)/r^2);    

end
