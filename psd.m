function [f,psdt] = psd(dt,xx,yy,dim,dis)
%------------------------------------------------------
% calculate PSD of single trajectory 
%------------------------------------------------------
% dt    time increment / frame time
% xx    x coordinates 
% yy    y coordinates 
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': frequencies equi-distr. on lin/log  
%------------------------------------------------------
%x(t) is real => X(-f)=(X(f))^* => one-sided PSD 

[x,y]    = dimcheck(xx,yy);
N        = numel(x);  %--> # positions of trajectory
m        = fix(N/2);
f        = (0:m)/(N*dt);
no       = dt/N; 

xf       = no*abs(fft(x)).^2;
psx      = xf(1:m+1);
psx(2:m) = 2*psx(2:m);

yf       = no*abs(fft(y)).^2;
psy      = yf(1:m+1);
psy(2:m) = 2*psy(2:m);

if (dim == 1)
    psdt = psx;
elseif (dim == 2)
    psdt = psy;
else
    psdt  = psx+psy;
end

if (dis == 'lin')
    ;
elseif (dis == 'log')
    df    = 1.2; 
    pf    = unique(round(df.^[1:round(log(m)/log(df))]));
    pf(pf > m)=[]; f=f(pf); psdt=psdt(pf);
else
    error('> psd stopped, wrong input parameter 5 <')
end
