function [xi,vacf] = vacf(dt,dn,xx,yy,dim,dis)
%------------------------------------------------------
% calculate normalized VACF of single trajectory 
%------------------------------------------------------
% dt    time increment / frame time
% dn    frame lag (\delta t=dn*\Delta t, frames i->i+dn)
% xx    x coordinates 
% yy    y coordinates 
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
%------------------------------------------------------
[x,y]    = dimcheck(xx,yy);
N        = numel(x);  %--> # positions of trajectory
dnt      = dn*dt;     %--> \delta t
zz       = N-dn;
vx(1:zz) = x(1+dn:N)-x(1:N-dn);
vy(1:zz) = y(1+dn:N)-y(1:N-dn);

if (dis == 'lin')
    m     = N/2;
    s     = 1:m;
elseif (dis == 'log')
    df    = 1.2; 
    s     = unique(round(df.^[1:round(log(zz/2)/log(df))]));
    m     = max(size(s));
else
    error('>>>> TA_VACF stopped, wrong input parameter 6 <<<')
end

tau   = s*dt;
xi    = tau/dnt;
vacfx = 0*(1:m);
vacfy = 0*(1:m);
vacf2 = 0*(1:m);

for j=1:m   
    vacfx(j) = mean(vx(1:zz-s(j)).*vx(1+s(j):zz));
    vacfy(j) = mean(vy(1:zz-s(j)).*vy(1+s(j):zz));
    vacf2(j) = mean(vx(1:zz-s(j)).*vx(1+s(j):zz)+vy(1:zz-s(j)).*vy(1+s(j):zz));
end

if (dim == 1)
    nox  = mean(vx.^2);
    vacf = vacfx/nox;
elseif (dim == 2)
    noy  = mean(vy.^2);
    vacf = vacfy/noy;
else
    noxy = mean(vx.^2+vy.^2);
    vacf = vacf2/noxy;
end