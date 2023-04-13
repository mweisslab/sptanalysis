function [fT,gam] = cov_gamma(dt,N,M,xx,yy,dim,dis)
%------------------------------------------------------
% calculate PSD coefficient of variation (gamma)
%------------------------------------------------------
% dt    time increment / frame time
% N     length of trajectories 
% M     ensemble size (number of trajectories
% xx    array of x coordinates
% yy    array of y coordinates
% dim   1: 1D (x) | 2: 1D (y) | else: 2D (x,y) 
% dis  'lin'/'log': frequencies equi-distr. on lin/log  
%------------------------------------------------------

[x,y] = enscheck(N,M,xx,yy);
psdt  =[];
if ((dim == 1) | (dim == 2))
    [f,psde] = ea_psd(dt,N,M,x,y,dim,'lin');
    for j=1:M
        [f,ppd]  = psd(dt,x(j,1:N),y(j,1:N),dim,'lin');
        psdt     = [psdt;ppd];
    end
else
    [f,psde] = ea_psd(dt,N,M,x,y,1,'lin');
    [f,ppd ] = ea_psd(dt,N,M,x,y,2,'lin');    
    psde     = (psde+ppd)/2;
    for j=1:M
        [f,ppd]  = psd(dt,x(j,1:N),y(j,1:N),1,'lin');
        psdt     = [psdt;ppd];
        [f,ppd]  = psd(dt,x(j,1:N),y(j,1:N),2,'lin');
        psdt     = [psdt;ppd];
    end
end

fT  = f*N*dt;
gam = std(psdt)./psde;

if (dis == 'log')
    df    = 1.15; 
    m     = fix(N/2);
    pf    = unique(round(df.^[1:round(log(m)/log(df))]));
    pf(pf > m)=[]; 
    fT    = fT(pf);
    ppd   = gam;
    gam   = ppd(1)+fT*0;
    for i=2:numel(pf)
        lb = pf(i-1);
        ub = pf(i);
        gam(i)=mean(ppd(lb:ub-1));
    end
end


