function pos = make_blurwalk(N,alpha,dx,np)
%------------------------------------------------------
% create 2D FBM trajectories with localization errors
%------------------------------------------------------
% N       number of positions
% alpha   anomaly exponent = 2*Hurst parameter
% dx      desired mean step length
% np      number of photons
%------------------------------------------------

nst   = 10;             %--> # substeps per frame
sig   = 0.22;           %--> width of PSF
niter = N*nst;
tra   = make_rndwalk(niter,alpha,dx/sqrt(nst^alpha)); %--> adapt step size
pos   = zeros(N,2);

for t=1:N
    xs=[]; ys=[];
    for i=1:nst
        xs=[xs,tra((t-1)*nst+i,1)+random('Normal',0,sig,1,np)];
        ys=[ys,tra((t-1)*nst+i,2)+random('Normal',0,sig,1,np)];
    end
    pos(t,1) = mean(xs);
    pos(t,2) = mean(ys);
end
