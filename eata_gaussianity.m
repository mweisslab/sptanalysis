function [tau,g] = eata_gaussianity(dt,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate EA-TA-Gaussianity of trajectory ensemble
%------------------------------------------------------
% dt    time increment / frame time
% N     length of trajectories 
% M     ensemble size (number of trajectories
% xx    array of x coordinates
% yy    array of y coordinates
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
%------------------------------------------------------

[x,y] = enscheck(N,M,xx,yy);
nmax  = N/2;      %--> max. lag time in units dt

if (dis == 'lin')
    m     = nmax;
    s     = 1:m;
elseif (dis == 'log')
    df    = 1.2; 
    s     = unique(round(df.^[1:round(log(nmax)/log(df))]));
    m     = max(size(s));
else
    error('> EATA_gaussianity stopped, wrong parameter 5 <')
end

tau    = s*dt;
msdte  = 0*(1:m);
quate  = 0*(1:m);

for j=1:M
    [tau,msd]  = ta_msd(dt,x(j,1:N),y(j,1:N),dim,dis);
    [tau,qua]  = ta_quad(dt,x(j,1:N),y(j,1:N),dim,dis);
    msdte      = msdte+msd;
    quate      = quate+qua;
end
msdte = msdte/M;
quate = quate/M;

if ((dim == 1) | (dim == 2))
    g  = quate./(3*msdte.^2)-1; 
else
    g  = 2*quate./(3*msdte.^2)-1; 
end
