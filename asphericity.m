function [An,Ad] = asphericity(xx,yy)
%------------------------------------------------------
% calculates  nominator & denominator of asphericity for a single 2D trajectory
%------------------------------------------------------
% xx    x coordinates 
% yy    y coordinates 
%------------------------------------------------------

[x,y] = dimcheck(xx,yy);
N     = numel(x);  %--> # positions of trajectory

rx  = x-mean(x);
ry  = y-mean(y);
a   = mean(rx.^2);  %matrix element 1,1
b   = mean(ry.^2);  %matrix element 2,2
c   = mean(rx.*ry); %matrix element 1,2
R2  = 0.5*(a+b)+sqrt(0.25d0*(a-b)^2+c^2)*[-1,+1];
An  = (R2(1)-R2(2))^2;
Ad  = (R2(1)+R2(2))^2;
%A = mean(zahl)/mean(nenn);

