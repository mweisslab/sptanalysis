function [tau,g] = ta_gaussianity(dt,xx,yy,dim,dis)
%------------------------------------------------------
% calculate TA-Gaussianity of single trajectory 
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
    error('> ta_gaussianity stopped, wrong parameter 5 <')
end

[tau,msdt]  = ta_msd(dt,x,y,dim,dis);
[tau,quat]  = ta_quad(dt,x,y,dim,dis);

if ((dim == 1) | (dim == 2))
    g  =   quat./(3*msdt.^2)-1; 
else
    g  = 2*quat./(3*msdt.^2)-1; 
end
