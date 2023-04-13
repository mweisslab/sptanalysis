function [tau,msde] = ea_msd(dt,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate EA-MSD of trajectory ensemble
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
    error('> EA_MSD stopped, wrong parameter 5 <')
end

tau   = s*dt;
msdex = 0*(1:m);
msdey = 0*(1:m);

for j=1:m   
    msdex(j) = mean((x(1:M,1+s(j))-x(1:M,1)).^2);
    msdey(j) = mean((y(1:M,1+s(j))-y(1:M,1)).^2);
end
    
if (dim == 1)
    msde = msdex;
elseif (dim == 2)
    msde = msdey;
else
    msde = msdex+msdey;
end