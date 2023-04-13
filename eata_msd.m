function [tau,msdte] = eata_msd(dt,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate EA-TA-MSD of trajectory ensemble
%------------------------------------------------------
% dt    time increment / frame time
% N     length of trajectories 
% M     ensemble size (number of trajectories
% xx    array of x coordinates
% yy    array of y coordinates
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
%------------------------------------------------------

[x,y]       = enscheck(N,M,xx,yy);
[tau,msdte] = ta_msd(dt,x(1,1:N),y(1,1:N),dim,dis);
for j=2:M
    [tau,msd]  = ta_msd(dt,x(j,1:N),y(j,1:N),dim,dis);
    msdte      = msdte+msd;
end
msdte = msdte/M;    
