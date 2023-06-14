function [tau,mom] = ta_mom(dt,xx,yy,q,dim,dis)
%------------------------------------------------------
% calculate TA k-th moment of single trajectory 
%------------------------------------------------------
% dt    time increment / frame time
% xx    x coordinates 
% yy    y coordinates 
% q     order of the moment to be calculated
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
%------------------------------------------------------


[x,y] = dimcheck(xx,yy);
N     = numel(x);  %--> # positions of trajectory
nmax  = N/2;       %--> max. lag time in units dt

if (dis == 'lin')
    m     = nmax;
    s     = 1:m;
elseif (dis == 'log')
    df    = 1.2; 
    s     = unique(round(df.^[1:round(log(nmax)/log(df))]));
    m     = max(size(s));
else
    error('> TA_MSD stopped, wrong parameter 5 <')
end

tau   = s*dt;
quatx = 0*(1:m);
quaty = 0*(1:m);

for j=1:m   
    quatx(j) = mean(abs(x(1:N-s(j))-x(1+s(j):N)).^q);
    quaty(j) = mean(abs(y(1:N-s(j))-y(1+s(j):N)).^q);
end
    
if (dim == 1)
    mom = quatx;
elseif (dim == 2)
    mom = quaty;
else
    mom = quatx+quaty;
end
