function [tau,acfe] = ea_acf_sqinc(dt,dn,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate EA-ACF of fluctuations in squared increments in dn frames)
%------------------------------------------------------
% dt    time increment / frame time
% dn    frame lag (differences between frames i, i+dn)
% N     length of trajectories 
% M     ensemble size (number of trajectories
% xx    array of x coordinates
% yy    array of y coordinates
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': times equi-distr. on lin/log  
%------------------------------------------------------

[x,y]      = enscheck(N,M,xx,yy);
[tau,acfe] = acf_sqinc(dt,dn,x(1,1:N),y(1,1:N),dim,dis);
for j=2:M
    [f,psd]  = acf_sqinc(dt,dn,x(j,1:N),y(j,1:N),dim,dis);
    acfe     = acfe+psd;
end
acfe = acfe/M;    
