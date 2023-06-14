function [tau] = passage_time(dt,xx,yy,r,dim)
%------------------------------------------------------
% get first passage time (-1 if distance r not reached)
%------------------------------------------------------
% dt    time increment / frame time
% xx    x coordinates 
% yy    y coordinates
% r     escape radius  
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
%------------------------------------------------------

[x,y] = dimcheck(xx,yy);
r2    = r^2;
tau   = -1;

if (dim == 1)
    dr2 = (x(2:end)-x(1)).^2;
elseif (dim == 2)
    dr2 = (y(2:end)-y(1)).^2;
else
    dr2 = (x(2:end)-x(1)).^2 + (y(2:N)-y(1)).^2;
end
jj  = min(find(dr2 > r2));

if (jj > 0)
    tau=dt*jj;
end

