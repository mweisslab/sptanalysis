function [xi,vacfte] = ea_vacf(dt,dn,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate EA-TA-VACF of trajectory ensemble
%------------------------------------------------------
% dt    time increment / frame time
% dn    frame lag (\delta t=dn*\Delta t, frames i->i+dn)
% N     length of trajectories 
% M     ensemble size (number of trajectories
% xx    array of x coordinates
% yy    array of y coordinates
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
%------------------------------------------------------

[x,y]       = enscheck(N,M,xx,yy);
[xi,vacfte] = vacf(dt,dn,x(1,1:N),y(1,1:N),dim,dis);
for j=2:M
    [xi,v]  = vacf(dt,dn,x(j,1:N),y(j,1:N),dim,dis);
    vacfte  = vacfte+v;
end
vacfte = vacfte/M;    
