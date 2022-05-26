function sigmar = stokes_dlp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, b1, b2, Lx, Ly, xi)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the stress directly. Considers one
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
%       sigmar, stress as a 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
sigmar = zeros(2, Ntar);

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
                    sigma_tmp =  stress_dlp_real_sum(r1,r2,n1(n),n2(n),f1(n),f2(n),b1(m),b2(m),xi);
                    sigmar(:,m) = sigmar(:,m) + sigma_tmp;
                end
            end
            
        end
    end
end
end

function sigma_real = stress_dlp_real_sum(r1, r2, n1, n2, f1, f2, b1, b2, xi)

sigma_real = [0;0];

rdotf = r1*f1 + r2*f2;
rdotn = r1*n1 + r2*n2;
rdotb = r1*b1 + r2*b2;
fdotb = f1*b1 + f2*b2;
bdotn = n1*b1 + n2*b2;
fdotn = n1*f1 + n2*f2;

rSq = r1^2 + r2^2;
xi2 = xi*xi;
e = exp(-xi2*rSq);
mu = 1;

sigma_real(1) = e*(b1/(2*pi)*((fdotn-2*xi2*rdotf*rdotn)/rSq-2*rdotf*rdotn/(rSq*rSq)) + ...
    1/(4*pi)*2*mu*(r1*rdotb*rdotf*rdotn*((8*xi2*xi2)/rSq+(16*xi2)/(rSq*rSq)+16/(rSq*rSq*rSq)) - ...
    2*(1+xi2*rSq)*(f1*rdotb*rdotn+n1*rdotb*rdotf+r1*fdotb*rdotn+r1*bdotn*rdotf+2*b1*rdotf*rdotn)/(rSq*rSq) + ...
    2*xi2*(b1*fdotn+n1*fdotb+f1*bdotn-xi2*(2*r1*rdotb*fdotn+r1*rdotn*fdotb+r1*rdotf*bdotn+f1*rdotb*rdotn+n1*rdotb*rdotf))));

sigma_real(2) = e*(b2/(2*pi)*((fdotn-2*xi2*rdotf*rdotn)/rSq-2*rdotf*rdotn/(rSq*rSq)) + ...
    1/(4*pi)*2*mu*(r2*rdotb*rdotf*rdotn*((8*xi2*xi2)/rSq+(16*xi2)/(rSq*rSq)+16/(rSq*rSq*rSq)) - ...
    2*(1+xi2*rSq)*(f2*rdotb*rdotn+n2*rdotb*rdotf+r2*fdotb*rdotn+r2*bdotn*rdotf+2*b2*rdotf*rdotn)/(rSq*rSq) + ...
    2*xi2*(b2*fdotn+n2*fdotb+f2*bdotn-xi2*(2*r2*rdotb*fdotn+r2*rdotn*fdotb+r2*rdotf*bdotn+f2*rdotb*rdotn+n2*rdotb*rdotf))));

end
