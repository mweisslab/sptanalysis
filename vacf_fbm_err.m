function [xi,vacf] = vacf_fbm_err(dn,theta,alpha)
%-------------------------------------------------
% normalized analytical VACF for FBM with loc. errors
%-------------------------------------------------
% dn      time steps used for velocity (\delta t= dn*\Delta t)
% anz     numer of lagtime steps in frame times dt 
% theta   rel. static loc. error, sigma^2/(K \Delta t)^alpha 
% alpha   anomaly exponent
%------------------------------------------------------

N    = 5*dn;
A    = (1:2*N+dn)*0;
xi   = (1:N)/dn;
vacf = 0*xi;
a1   = alpha+1;
a2   = alpha+2;
for m=1:2*N+dn
    A(m) = (m+1)^a2 + (m-1)^a2 - 2*m^a2;
end
for m=1:N
    nenn = 2*A(dn)-4+2*theta*a1*a2;
    if (m == dn)
        zahl = A(2*dn)-2*A(dn)+2-theta*a1*a2;
    else
        zahl = A(m+dn)-2*A(m)+A(abs(m-dn));
    end    
    vacf(m) = zahl/nenn;
end   

