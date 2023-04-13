function [xi,vacf] = vacf_fbm_theo(alpha)
%-------------------------------------------------
% normalized analytical VACF for FBM 
%-------------------------------------------------
% alpha  anomaly exponent 
%-------------------------------------------------

if ((alpha <= 0) | (alpha>= 2))
    error('choose 0<alpha<2');
end
xi   = 0.01*(1:500); 
vacf = 0.5*((xi+1).^alpha+abs(xi-1).^alpha-2*xi.^alpha);





