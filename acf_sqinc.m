function [tau,acf] = acf_sqinc(dt,dn,xx,yy,dim,dis)
%------------------------------------------------------
% ACF of fluctuations in squared steps taken in dn frames, single trajectory 
%------------------------------------------------------
% dt    time increment / frame time
% dn    frame lag (differences between frames i, i+dn)
% xx    x coordinates 
% yy    y coordinates 
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': lag times equi-distr. on lin/log  
% NB: code avoids correlated values for dn>1 by superimposing
%     individual realizations within interval 1:dn
%------------------------------------------------------

[x,y] = dimcheck(xx,yy);
b     = mod(length(x),dn);
N     = numel(x)-b;  %--> # usable positions of trajectory 

for loop=1:dn
    sx  = (x(loop+dn:dn:N)-x(loop:dn:N-dn)).^2;
    sy  = (y(loop+dn:dn:N)-y(loop:dn:N-dn)).^2;
    if (dim == 1)
        dr2 = sx;
    elseif (dim == 2)
        dr2 = sy;
    else
        dr2 = sx+sy;        
    end
    dr2 = dr2/mean(dr2)-1;
    if (loop == 1)
        m    = numel(dr2);
        tau  = dn*dt*(1:m-1);
        acf  = 0*(1:m-1);
    end
    for j=1:m-1   
        acf(j)=acf(j)+mean(dr2(1+j:m).*dr2(1:m-j));
    end
end
acf=acf/dn;

if (dis == 'lin')
    ;
elseif (dis == 'log')
    df    = 1.2; 
    pf    = unique(round(df.^[1:round(log(m-1)/log(df))]));
    pf(pf > m-1)=[];
    tau=tau(pf);acf=acf(pf);
else
    error('> acf_sq_inc stopped, wrong input parameter 6 <<<')
end




