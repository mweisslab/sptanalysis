function S = straightness(xx,yy,lb,ub)
%------------------------------------------------------
% calculate trajectory straightness between points lb,ub
%------------------------------------------------------
% xx    x coordinates of whole trajectory 
% yy    y coordinates of whole trajectory 
% lb    first position to be included
% ub    last position  to be included
%------------------------------------------------------

[x,y] = dimcheck(xx,yy);
N     = numel(x);  %--> # positions of trajectory
S     = 0;

if ((ub <= lb) | (ub > N) | (lb > N) | (lb < 1) | (ub < 1))
    error('lb,ub: wrong/mising positions in trajectory');
end

za = sqrt((x(ub)-x(lb))^2+(y(ub)-y(lb))^2);
ne = sum(sqrt((x(lb+1:ub)-x(lb:ub-1)).^2+(y(lb+1:ub)-y(lb:ub-1)).^2));

if (ne > 0) 
    S=za/ne; 
end