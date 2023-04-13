function [dx,dy] = get_increm(dn,xx,yy,chi)
%-------------------------------------------------
% calculate step increments in single trajectory 
%-------------------------------------------------
% dn    frame lag (steps between frames i, i+dn)
% xx    x coordinates 
% yy    y coordinates 
% chi   true/false: normalize y/n
%-------------------------------------------------

[x,y] = dimcheck(xx,yy);
N     = numel(x);  %--> # positions of trajectory
dx    = (x(1+dn:N)-x(1:N-dn));
dy    = (y(1+dn:N)-y(1:N-dn));
if (chi)
    dx=dx/sqrt(mean(dx.^2));
    dy=dy/sqrt(mean(dy.^2));    
end
