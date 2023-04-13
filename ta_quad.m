function [tau,quat] = ta_quad(dt,xx,yy,dim,dis)
%------------------------------------------------------
% calculate TA 4th moment of single trajectory 
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
quatx = 0*(1:m);
quaty = 0*(1:m);

for j=1:m   
    quatx(j) = mean((x(1:N-s(j))-x(1+s(j):N)).^4);
    quaty(j) = mean((y(1:N-s(j))-y(1+s(j):N)).^4);
end
    
if (dim == 1)
    quat = quatx;
elseif (dim == 2)
    quat = quaty;
else
    quat = quatx+quaty;
end