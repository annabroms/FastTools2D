function ugrad_r = stokes_dlp_gradient_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, b1, b2, Lx, Ly, xi)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the Stokeslet directly. Considers one
% periodic replicate box in each direction.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points
%       xsrc, x component of target points
%       ytar, y component of target points
%       n1, x component of normal vector
%       n2, y component of normal vector
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
                    ugrad_tmp =  stresslet_gradient_real_sum(r1,r2,n1(n),n2(n),f1(n),f2(n),b1(m),b2(m),xi);
                    ugrad_r(:,m) = ugrad_r(:,m) + ugrad_tmp;
                end
            end
            
        end
    end
end

ugrad_r = ugrad_r / (4*pi);

end

function ugrad_real = stresslet_gradient_real_sum(r1, r2, n1, n2, f1, f2, b1, b2, a)

ugrad_real = [0;0];

rdotf = r1*f1 + r2*f2;
rdotn = r1*n1 + r2*n2;
rdotb = r1*b1 + r2*b2;
bdotf = f1*b1 + f2*b2;
bdotn = n1*b1 + n2*b2;
fdotn = n1*f1 + n2*f2;

r = sqrt(r1^2 + r2^2);

ugrad_real(1) = exp(-a^2*r^2)*(r1*rdotf*rdotn*rdotb*(8*a^4/r^2+16*a^2/r^4+16/r^6)...
             -4*(1+a^2*r^2)*(b1*rdotf*rdotn+r1*bdotf*rdotn+r1*bdotn*rdotf)/r^4 ...
            +2*a^2*(f1*bdotn+n1*bdotf+b1*fdotn-2*a^2*rdotb*(f1*rdotn+n1*rdotf+r1*fdotn)));

ugrad_real(2) = exp(-a^2*r^2)*(r2*rdotf*rdotn*rdotb*(8*a^4/r^2+16*a^2/r^4+16/r^6)...
             -4*(1+a^2*r^2)*(b2*rdotf*rdotn+r2*bdotf*rdotn+r2*bdotn*rdotf)/r^4 ...
             +2*a^2*(f2*bdotn+n2*bdotf+b2*fdotn-2*a^2*rdotb*(f2*rdotn+n2*rdotf+r2*fdotn)));

end
