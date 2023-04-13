function [xx,yy] = dimcheck(xx,yy)
%----------------------------------
% check proper array dimensions 
%----------------------------------
nx = size(xx); 
ny = size(yy);

%--> check array dimensions of xx, yy 
if ((ndims(xx) ~= 2) | (min(nx) ~= 1) | (max(nx) < 2) | ...
    (ndims(yy) ~= 2) | (min(ny) ~= 1) | (max(ny) < 2) | (max(nx) ~=max(ny))) 
    error('xx,yy: wrong array dimensions, both must be (1:N,1) or (1,1:N)'); 
end

%--> force xx and yy to be row vectors (1,1:N)
if (size(xx,1) > size(xx,2))
    xx=xx';
end
if (size(yy,1) > size(yy,2))
    yy=yy';
end
