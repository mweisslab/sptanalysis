function [tau,Sd] = lch(dt,dn,xx,yy)
%------------------------------------------------------
% max. normalized diameter of LCH, single trajectory 
%------------------------------------------------------
% dt    time increment / frame time
% dn    integration window (# vertices for convex hull)
% xx    x coordinates 
% yy    y coordinates 
%------------------------------------------------------
[x,y]  = dimcheck(xx,yy);
N      = numel(x);     %--> # positions of trajectory
tau    = [];
Sd     = [];
for i=dn+1:N-dn
    sx = x(i-dn:i+dn);
    sy = y(i-dn:i+dn);
    Dt = convhull(sx,sy);
    nn = numel(Dt);
    ab = 0;
    for jj=1:nn
        pf = Dt(jj+1:nn);
        aa = max(sqrt((sx(Dt(jj))-sx(pf)).^2+(sy(Dt(jj))-sy(pf)).^2));    
        if (aa > ab) ab=aa; end;   
    end
    Sd  = [Sd,ab];
    tau = [tau,dt*i];
end
Sd=Sd/mean(Sd);