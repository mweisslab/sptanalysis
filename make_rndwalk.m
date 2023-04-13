function pos = make_rndwalk(N,alpha,dx)
%------------------------------------------------
% create FBM random walk in 2D (independent x,y)
%------------------------------------------------
% N       number of positions
% alpha   anomaly exponent = 2*Hurst parameter
% dx      desired mean step length
%------------------------------------------------

inc          = zeros(N-1,2);
pos          = zeros(N,2);
%---> create FBM, get increments
xx           = wfbm(alpha/2,N); 
yy           = wfbm(alpha/2,N); 
inc(1:N-1,1) = xx(2:N)-xx(1:N-1);
inc(1:N-1,2) = yy(2:N)-yy(1:N-1);
%---> scale increments, create trajectory 
mwx          = mean(inc(1:N-1,1).^2);
mwy          = mean(inc(1:N-1,2).^2);
pos(2:N,1)   = cumsum(inc(1:N-1,1))*dx/sqrt(mwx);
pos(2:N,2)   = cumsum(inc(1:N-1,2))*dx/sqrt(mwy);

