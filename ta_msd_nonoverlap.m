function [tau,msdt] = ta_msd_nonoverlap(dt,xx,yy,dim,dis)
%------------------------------------------------------
% calculate TA-MSD of single trajectory w/o overlappings
%------------------------------------------------------
% dt    time increment / frame time
% xx    x coordinates 
% yy    y coordinates 
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
msdtx = 0*(1:m);
msdty = 0*(1:m);

for j=1:m   
    msdtx(j) = mean((x(1:s(j):N-s(j))-x(1+s(j):s(j):N)).^2);
    msdty(j) = mean((y(1:s(j):N-s(j))-y(1+s(j):s(j):N)).^2);
end
    
if (dim == 1)
    msdt = msdtx;
elseif (dim == 2)
    msdt = msdty;
else
    msdt = msdtx+msdty;
end