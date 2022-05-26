function omega_r = stokes_dlp_vorticity_real_ds(xsrc, ysrc, xtar, ytar,...
                        n1, n2, f1, f2, Lx, Ly, xi)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Evaluates the real space sum of the vorticity directly. Considers one
% periodic replicate box in each direction.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points
%       xsrc, x component of target points
%       ytar, y component of target points
%       n1, x component of normal vector at source points
%       n2, y component of normal vecotr at source points
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
                    continue;
                else
                    omega_tmp = vorticity_dlp_real_sum(r1,r2,n1(n),n2(n),f1(n),f2(n),xi);
                    omega_r(m) = omega_r(m) + omega_tmp;
                end
            end
            
        end
    end
end

end

function omega_real = vorticity_dlp_real_sum(r1, r2, n1, n2, f1, f2, xi)
rdotf = r1*f1 + r2*f2;
rdotn = r1*n1 + r2*n2;

rdot = r1*(n2*rdotf+f2*rdotn) - r2*(n1*rdotf+f1*rdotn);
ndot = n1*(r2*rdotf) - n2*(r1*rdotf);
fdot = f1*(r2*rdotn) - f2*(r1*rdotn);

rSq = r1^2 + r2^2;
xi2 = xi*xi;

omega_real = exp(-xi2*rSq)*((1+xi2*rSq)*rdot/(rSq*rSq) + xi2*xi2*(ndot+fdot))/pi;

% April 2022 (gives exactly the same as above)
% rdotf = r1*f1 + r2*f2;
% rdotn = r1*n1 + r2*n2;
% 
% ndotrperp = n1*r2*rdotf + n2*(-r1)*rdotf;
% fdotrperp = f1*r2*rdotn + f2*(-r1)*rdotn;
% 
% rdotnperp = r1*n2*rdotf + r2*(-n1)*rdotf;
% rdotfperp = r1*f2*rdotn + r2*(-f1)*rdotn;
% 
% omega_real_new = -exp(-xi2*rSq)*((1+xi2*rSq)*(ndotrperp+fdotrperp)/(rSq*rSq) + ...
%     xi2*xi2*(rdotnperp+rdotfperp))/pi;

end



