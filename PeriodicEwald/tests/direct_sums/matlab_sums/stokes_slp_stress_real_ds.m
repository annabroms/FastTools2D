function sigmar = stokes_slp_stress_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, b1, b2, Lx, Ly, xi)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the stresslet. Considers one
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
%       sigmar, stress as an 2xN matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Nsrc = length(xsrc);
Ntar = size(xtar,1);
sigmar = zeros(2, Ntar);

for n=1:Nsrc
    
    for m=1:Ntar 
        
        for jpx = -1:1 %Go through a layer of boxes in the x direction
            
            for jpy = -1:1 %Go through a layer of boxes in the y direction
                
                %Compute periodic source point
                xsrc_p = xsrc(n) + jpx*Lx;
                ysrc_p = ysrc(n) + jpy*Ly;
                
                r1 = xtar(m) - xsrc_p;
                r2 = ytar(m) - ysrc_p;
                r = sqrt(r1^2 + r2^2);
                               
                % if r == 0, skip
                if abs(r) < 1e-13
                    continue
                else
                    sigma_tmp =  stresslet_real_sum(r1,r2,...
                                f1(n),f2(n),b1(m),b2(m),xi);
                    sigmar(:,m) = sigmar(:,m) + sigma_tmp;
                end
            end
            
        end
    end
end

sigmar = sigmar / (4*pi);

end

function sigma_real = stresslet_real_sum(r1, r2, f1, f2, b1, b2, xi)

sigma_real = [0;0];

rdotf = r1*f1 + r2*f2;
rdotb = r1*b1 + r2*b2;
bdotf = f1*b1 + f2*b2;

rSq = r1*r1 + r2*r2;
xi2 = xi*xi;

sigma_real(1) = exp(-xi2*rSq)*(2*xi2*(b1*rdotf+f1*rdotb+r1*bdotf) - ...
    (4*r1*rdotb*rdotf/(rSq*rSq))*(1+xi2*rSq));

sigma_real(2) = exp(-xi2*rSq)*(2*xi2*(b2*rdotf+f2*rdotb+r2*bdotf) - ...
    (4*r2*rdotb*rdotf/(rSq*rSq))*(1+xi2*rSq));

end
