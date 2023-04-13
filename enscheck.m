function [xx,yy] = enscheck(N,M,xx,yy)
%------------------------------------------------
% check proper dimensions for trajectory ensemble
%------------------------------------------------

pf = sort([N,M]);
nx = sort(size(xx));
ny = sort(size(yy));

%--> check array dimensions of xx, yy 
if ((ndims(xx) ~= 2) | (ndims(yy) ~= 2) | any(nx ~= pf) | any(ny ~= pf))
    error('xx,yy: wrong array dimensions, must be (1:N,1:M) or (1:M,1:N)'); 
end

if (N ~= M)
    %--> force xx and yy to have trajectories as row vectors (.,1:N)
    if (size(xx,1) == N)
        xx=xx';
    end
    if (size(yy,1) == N)
        yy=yy';
    end
else
    display(['CAUTION: M=N, no self-check that trajectories are row vectors!'])
end