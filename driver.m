clear; close all; 
%##########################################################
% test routine of the analysis package
%##########################################################

%%%%%%%%%%%  from here on: create data sets %%%%%%%%%%%%%
%---> uncomment those parts that you need to run
dt    = 0.1;    %--> frame time in [s]
dx    = 0.01;   %--> average step size in [um]
N     = 500;    %--> trajectory length
M     = 100;    %--> ensemble size
alpha = 0.6;    %--> anomaly exponent (=2*H) 
np    = 900;
ext   = num2str(alpha);
if (alpha == 1.0) ext='1.0'; end;

%%---> create & store  2D FBM trajectories
% $$$ y     = zeros(M,N,2);
% $$$ for loop=1:M
% $$$     y(loop,1:N,1:2) = make_rndwalk(N,alpha,dx); 
% $$$ end
% $$$ xx = y(:,:,1);
% $$$ yy = y(:,:,2);
% $$$ dlmwrite(['FBM_xx_a_',ext,'.dat'],xx,'delimiter','\t');
% $$$ dlmwrite(['FBM_yy_a_',ext,'.dat'],yy,'delimiter','\t');

%%---> create & store intermittent 2D FBM trajectories
% $$$ xx       = dlmread(['./data/FBM_xx_a_0.6.dat']);
% $$$ yy       = dlmread(['./data/FBM_yy_a_0.6.dat']);
% $$$ kon      = 0.4;
% $$$ koff     = 0.1;
% $$$ fra      = 2;
% $$$ [xx,yy]  = make_switchwalk(dt,N,M,xx,yy,kon,koff,fra);
% $$$ dlmwrite(['FBM_switch_xx_a_0.6.dat'],xx,'delimiter','\t');
% $$$ dlmwrite(['FBM_switch_yy_a_0.6.dat'],yy,'delimiter','\t');
% $$$ 
% $$$ xx       = dlmread(['./data/FBM_xx_a_1.0.dat']);
% $$$ yy       = dlmread(['./data/FBM_yy_a_1.0.dat']);
% $$$ kon      = 0.2;
% $$$ koff     = 0.3;
% $$$ fra      = 3;
% $$$ [xx,yy]  = make_switchwalk(dt,N,M,xx,yy,kon,koff,fra);
% $$$ dlmwrite(['FBM_switch_xx_a_1.0.dat'],xx,'delimiter','\t');
% $$$ dlmwrite(['FBM_switch_yy_a_1.0.dat'],yy,'delimiter','\t');

%%---> create & store 2D FBM trajectories with localization errors
% $$$ y     = zeros(M,N,2);
% $$$ for loop=1:M
% $$$     y(loop,1:N,1:2) = make_blurwalk(N,alpha,dx,np); 
% $$$ end
% $$$ xx = y(:,:,1);
% $$$ yy = y(:,:,2);
% $$$ dlmwrite(['FBM_xx_a_',ext,'_np_',num2str(np),'.dat'],xx,'delimiter','\t');
% $$$ dlmwrite(['FBM_yy_a_',ext,'_np_',num2str(np),'.dat'],yy,'delimiter','\t');
% $$$ 

%%%%%%%%%%%  from here on: use of tools, production of figures %%%%%%%%
%%---> Fig. 1A,B,C
cw = cos(2*pi*(0:500)/500);
sw = sin(2*pi*(0:500)/500);
for loop=1:3
    if (loop == 1) df=0.15; alpha=0.6; ext='0.6'; end
    if (loop == 2) df=0.45; alpha=1.0; ext='1.0'; end    
    if (loop == 3) df=1.2;  alpha=1.4; ext='1.4'; end    
    xx          = dlmread(['data/FBM_xx_a_',ext,'.dat']);
    yy          = dlmread(['data/FBM_yy_a_',ext,'.dat']);
    [tau,msdte] = eata_msd(1,N,M,xx,yy,0,'lin');
    plot(xx(1,:),yy(1,:),'k','LineWidth',2);
    hold on;
    plot(xx(3,:),yy(3,:),'r','LineWidth',2);
    plot(xx(5,:),yy(5,:),'b','LineWidth',2);
    plot(xx(6,:),yy(6,:),'g','LineWidth',2);
    plot(sqrt(msdte(50))*cw,sqrt(msdte(50))*sw,':','Color',0.5*[1,1,1],'LineWidth',4)    
    plot(sqrt(msdte(200))*cw,sqrt(msdte(200))*sw,':','Color',0.5*[1,1,1],'LineWidth',4)    
    plot(sqrt(msdte(250)*2^alpha)*cw,sqrt(msdte(250)*2^alpha)*sw,':', 'Color',0.5*[1,1,1],'LineWidth',4)    
    plot([df-0.15,df-0.05],-0.9*df*[1,1],'k','LineWidth',5);
    hold off
    pbaspect([1 1 1]);
    axis(df*[-1 1 -1 1]);
    %saveas(1,['fig01a_FBM_a_',num2str(alpha),'.svg'],'svg')
    close(1)
end

%%---> Fig. 1D
ax=plotter([7e-2 4e1 1e-4 2e1],'\tau [s]','\langler^2(\tau)\rangle_t [\mum^2]','xlog','ylog');
for loop=1:3
    if (loop == 1) alpha=0.6; ext='0.6'; sh=1; vf=0.5; end
    if (loop == 2) alpha=1.0; ext='1.0'; sh=5; vf=0.45; end    
    if (loop == 3) alpha=1.4; ext='1.4'; sh=25;vf=2; end    
    xx          = dlmread(['./data/FBM_xx_a_',ext,'.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'.dat']);
    [taut,msdt] = ta_msd(dt,xx(1,:),yy(1,:),0,'log');
    [tau,msdte] = eata_msd(dt,N,M,xx,yy,0,'log');
    tt          = [2*taut(1),taut(end)/2];
    oplot(tau,sh*msdte,'-rs','none',10,2);
    oplot(taut,sh*msdt,'ko','none',15,2);
    oplot(tt,sh*vf*msdt(1)*(tt/dt).^alpha,'--k','none',10,3);
end
yticks([1e-4 1e-3 1e-2 1e-1 1e0 1e1])
hold off 
%saveas(1,'fig01d.svg','svg')
close(1)

%%---> Fig. 1E
y1       = dlmread(['./data/FBM_xx_a_0.6.dat']);
y2       = dlmread(['./data/FBM_yy_a_1.0.dat']);
y3       = dlmread(['./data/FBM_yy_a_1.4.dat']);
[tau,E1] = eb(dt,N,M,y1,y1,1,'log');
[tau,E2] = eb(dt,N,M,y2,y2,1,'log');
[tau,E3] = eb(dt,N,M,y3,y3,1,'log');
ax=plotter([0.2 40 3e-3 2],'\tau [s]','E(\tau)','xlog','ylog');
oplot(tau,E1,'rd','none',12,2);
oplot(tau,E2,'ko','none',12,2);
oplot(tau,E3,'bs','none',12,2);
oplot([0.2,40],4*[0.2,40]/(3*N*dt),'--k','none',10,3);
%saveas(1,'fig01e.svg','svg')
close(1)

%%---> Fig. 2
ax=plotter([0 4 -0.4 1.5],'\xi','C(\xi)','xlin','ylin');
for loop=1:3
    if (loop == 1) alpha=0.6; ext='0.6'; sh=0;   s='r'; end
    if (loop == 2) alpha=1.0; ext='1.0'; sh=0.2; s='k'; end    
    if (loop == 3) alpha=1.4; ext='1.4'; sh=0.4; s='b'; end    
    xx          = dlmread(['./data/FBM_xx_a_',ext,'.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'.dat']);
    [xit,vacft] = vacf(dt,3,xx(1,:),yy(1,:),0,'lin');
    [xi3,vacf3] = ea_vacf(dt,3,N,M,xx,yy,0,'lin');
    [xi5,vacf5] = ea_vacf(dt,5,N,M,xx,yy,0,'lin');
    [xi7,vacf7] = ea_vacf(dt,7,N,M,xx,yy,0,'lin');
    [tt,vv]     = vacf_fbm_theo(alpha);
    oplot(tt,vv+sh,['-',s],'none',12,3);               
    oplot(xit,vacft+sh,['*',s],'none',12,2);
    oplot(xi3,vacf3+sh,['s',s],'none',12,2);
    oplot(xi5,vacf5+sh,['d',s],'none',12,2);
    oplot(xi7,vacf7+sh,['o',s],'none',12,2);
end
%saveas(1,'fig02.svg','svg')
close(1)

%%---> Fig. 3A
xx1 = dlmread(['./data/FBM_xx_a_1.4.dat']);
yy1 = dlmread(['./data/FBM_yy_a_1.4.dat']);
xx2 = dlmread(['./data/FBM_xx_a_1.0.dat']);
yy2 = dlmread(['./data/FBM_yy_a_1.0.dat']);
xx  = dlmread(['./data/FBM_xx_a_0.6.dat']);
yy  = dlmread(['./data/FBM_yy_a_0.6.dat']);

An=[]; Ad=[]; S1=[];
Bn=[]; Bd=[]; S2=[];
Cn=[]; Cd=[]; S3=[];

lb=1;ub=50;
for i=1:M
    S1=[S1,straightness(xx(i,lb:ub),yy(i,lb:ub),lb,ub)];
    S2=[S2,straightness(xx2(i,lb:ub),yy2(i,lb:ub),lb,ub)];        
    S3=[S3,straightness(xx1(i,lb:ub),yy1(i,lb:ub),lb,ub)];            
    [aa,bb]=asphericity(xx(i,:),yy(i,:));
    An=[An,aa]; Ad=[Ad,bb]; 
    [aa,bb]=asphericity(xx2(i,:),yy2(i,:));
    Bn=[Bn,aa]; Bd=[Bd,bb];
    [aa,bb]=asphericity(xx1(i,:),yy1(i,:));
    Cn=[Cn,aa]; Cd=[Cd,bb];
end
A1 = mean(An)/mean(Ad); sA1=[];
A2 = mean(Bn)/mean(Bd); sA2=[];
A3 = mean(Cn)/mean(Cd); sA3=[];
for i=1:5
    sA1=[sA1,mean(An((i-1)*20+1:i*20))/mean(Ad((i-1)*20+1:i*20))];
    sA2=[sA2,mean(Bn((i-1)*20+1:i*20))/mean(Bd((i-1)*20+1:i*20))];
    sA3=[sA3,mean(Cn((i-1)*20+1:i*20))/mean(Bd((i-1)*20+1:i*20))];
end
sx = [0.6,1,1.4];
sy = [A1,A2,A3];
sz = [mean(S1),mean(S2),mean(S3)];
ax = plotter([0.5,1.5,-0.1,0.9],'\alpha','A_e, \langleS\rangle','xlin','ylin');
oplot(sx,sy,'ko','none',15,2);
errorbar(sx,sy,[std(sA1),std(sA2),std(sA2)],'k.','LineWidth',3,'CapSize',18);
oplot([0.9,1.1],4*[1,1]/7,'--k','none',15,2);               
oplot(sx,sz,'rd','none',15,2);               
errorbar(sx,sz,[std(S1),std(S2),std(S3)],'r.','LineWidth',3,'CapSize',18);
hold off
xticks([0.6,1.0,1.4])
%saveas(1,'fig03a.svg','svg')
close(1)

%%---> Fig. 3B
xx1  = dlmread(['./data/FBM_xx_a_0.6.dat']);
yy1  = dlmread(['./data/FBM_yy_a_0.6.dat']);
xx2  = dlmread(['./data/FBM_switch_xx_a_0.6.dat']);
yy2  = dlmread(['./data/FBM_switch_yy_a_0.6.dat']);
mm   = 7;
who  = 3;
[tau,Sd1]=lch(dt,mm,xx1(who,:),yy1(who,:));
[tau,Sd2]=lch(dt,mm,xx2(who,:),yy2(who,:));
ax = plotter([0,50,0,5],'\tau [s]','S_d/\langleS_d\rangle','xlin','ylin');
oplot(tau,Sd1/mean(Sd1),'b','none',15,2)
oplot(tau,2+Sd2/mean(Sd2),'r','none',15,2);
oplot([0,50],[1,1],'--k','none',15,4);
oplot([0,50],[3,3],'--k','none',15,4);
%saveas(1,'fig03b.svg','svg')
close(1)    

%%---> Fig. 4A
xx  = dlmread(['./data/FBM_xx_a_0.6.dat']);
yy  = dlmread(['./data/FBM_yy_a_0.6.dat']);
xx1 = dlmread(['./data/FBM_switch_xx_a_0.6.dat']);
yy1 = dlmread(['./data/FBM_switch_yy_a_0.6.dat']);
x   = -5+10*(1:1000)/1000;
dx1 = [];  dy1 = [];
dx5 = [];  dy5 = [];
for j=1:M
    [ddx,ddy] = get_increm(1,xx(j,:),yy(j,:),true);
    dx1       = [dx1,ddx,ddy];
    [ddx,ddy] = get_increm(5,xx(j,:),yy(j,:),true);
    dx5       = [dx5,ddx,ddy];

    [ddx,ddy] = get_increm(1,xx1(j,:),yy1(j,:),true);
    dy1       = [dy1,ddx,ddy];
    [ddx,ddy] = get_increm(5,xx1(j,:),yy1(j,:),true);
    dy5       = [dy5,ddx,ddy];
end
ax = plotter([-5,5,1e-4,1],'\chi','p(\chi)','xlin','ylog');
xticks(-5:4);
oplot(x,exp(-x.^2/2)/sqrt(2*pi),'--k','none',10,3);
[yh,xh] = histcounts(dx1,'BinMethod','scott','Normalization','pdf');
xh      = xh+(xh(2)-xh(1))/2;
xh      = xh(1:end-1);
oplot(xh,yh,'bo','none',10,2);

[yh,xh] = histcounts(dx5,'BinMethod','scott','Normalization','pdf');
xh      = xh+(xh(2)-xh(1))/2;
xh      = xh(1:end-1);
oplot(xh,yh,'bs','none',10,2);

[yh,xh] = histcounts(dy1,'BinMethod','scott','Normalization','pdf');
xh      = xh+(xh(2)-xh(1))/2;
xh      = xh(1:end-1);
oplot(xh,yh,'ro','none',10,2);

[yh,xh] = histcounts(dy5,'BinMethod','scott','Normalization','pdf');
xh      = xh+(xh(2)-xh(1))/2;
xh      = xh(1:end-1);
oplot(xh,yh,'rs','none',10,2);
%saveas(1,'fig04a.svg','svg')
close(1)    

%%---> Fig. 4B 
[tau,ge]   = eata_gaussianity(dt,N,M,xx,yy,0,'log');
[tau1,ge1] = eata_gaussianity(dt,N,M,xx1,yy1,0,'log');
ax = plotter([0.05,40,-0.1,0.6],'\tau [s]','g(\tau)','xlog','ylin');
oplot(tau,ge,'-bo','none',10,2);
oplot(tau1,ge1,'-rs','none',10,2);
%saveas(1,'fig04b.svg','svg')
close(1)    

%%---> Fig. 4C
[tau,ge]  = ea_acf_sqinc(dt,1,N,M,xx,yy,0,'log');
[tau1,ge1] = ea_acf_sqinc(dt,1,N,M,xx1,yy1,0,'log');
ax = plotter([0.05,50,-0.1,0.6],'\tau [s]','G(\tau)','xlog','ylin');
oplot(tau,ge,'-bo','none',10,2);
oplot(tau1,ge1,'-rs','none',10,2);
%saveas(1,'fig04c.svg','svg')
close(1)

%%---> Fig. 5A
ax = plotter([1e-2 1e1 1e-5 2e0],'f [Hz]','S(f)','xlog','ylog');
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
oplot([3e-2,2],5e-4./[3e-2,2].^2,'-.k','none',10,3);
oplot([2e-2,1],3e-5./[2e-2,1].^1.6,'--k','none',10,3);
text(4e-2,5e-4,'\sim1/\tau^{0.6}','FontSize',20);
for loop=1:3
    if (loop == 1) alpha=0.6; ext='0.6'; s='rs'; ss='none'; end
    if (loop == 2) alpha=1.0; ext='1.0'; s='ko'; ss='none'; end    
    if (loop == 3) alpha=1.4; ext='1.4'; s='bo'; ss='b'; end    
    xx       = dlmread(['./data/FBM_xx_a_',ext,'.dat']);
    yy       = dlmread(['./data/FBM_yy_a_',ext,'.dat']);
    [f,psde] = ea_psd(dt,N,M,xx,yy,0,'log');
    oplot(f,psde,s,ss,12,2);
end
yticks([1e-4 1e-3 1e-2 1e-1 1e0 1e1])
%saveas(1,'fig05a.svg','svg')
close(1)

%%---> Fig. 5B
ax = plotter([1e-2 1e1 1e-6 1e-1],'f [Hz]','S(f)','xlog','ylog');
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
oplot([2e-2,3],3e-5./[2e-2,3].^1.6,'--k','none',10,3);
text(4e-2,3e-4,'\sim1/\tau^{0.6}','FontSize',20)

alpha    = 0.6; 
xx       = dlmread('./data/FBM_xx_a_0.6.dat');
yy       = dlmread('./data/FBM_yy_a_0.6.dat');
for i=1:5
    [f,psdt] = psd(dt,xx(i,:),yy(i,:),0,'log');  
    oplot(f,psdt,'-k','none',12,1);
end
[f,psde] = ea_psd(dt,N,M,xx,yy,0,'log');
oplot(f,psde,'rs','none',12,2);
%saveas(1,'fig05b.svg','svg')
close(1)
%%---> Fig.5c
ax = plotter([1 300 0.8 1.5],'fT','\gamma','xlog','ylin');
oplot([1,3e2],[1,1],'--r','none',12,3);
oplot([1,3e2],sqrt(5)*0.5*[1,1],'--k','none',12,3);
oplot([1,3e2],sqrt(2)*[1,1],'--b','none',12,3);

for loop=1:3
    if (loop == 1) alpha=0.6; ext='0.6'; s='rs'; ss='none'; end
    if (loop == 2) alpha=1.0; ext='1.0'; s='ko'; ss='none'; end    
    if (loop == 3) alpha=1.4; ext='1.4'; s='bo'; ss='b'; end    
    xx       = dlmread(['./data/FBM_xx_a_',ext,'.dat']);
    yy       = dlmread(['./data/FBM_yy_a_',ext,'.dat']);
    [fT,gam] = cov_gamma(dt,N,M,xx,yy,0,'log');
    oplot(fT,gam,s,ss,12,2);                
end
%saveas(1,'fig05c.svg','svg')
close(1)

%%---> Fig.6 all
theta=2;
for loop=1:3
    if (loop == 1) alpha=0.6; ext='0.6';  end
    if (loop == 2) alpha=1.0; ext='1.0';  end    
    if (loop == 3) alpha=1.4; ext='1.4';  end    
    ax=plotter([7e-2 1e1 1e-5 5e-2],'\tau [s]','\langler^2(\tau)\rangle_t [\mum^2]','xlog','ylog');
    xx          = dlmread(['./data/FBM_xx_a_',ext,'.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'.dat']);
    [te,msdte]  = eata_msd(dt,N,M,xx,yy,0,'log');
    xx          = dlmread(['./data/FBM_xx_a_',ext,'_np_900.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'_np_900.dat']);
    [td,msdtd]  = eata_msd(dt,N,M,xx,yy,0,'log');
    xx          = dlmread(['./data/FBM_xx_a_',ext,'_np_50.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'_np_50.dat']);
    [ts,msdts]  = eata_msd(dt,N,M,xx,yy,0,'log');
    tt          = [5*taut(1),taut(end)/3];
    oplot(te,sh*msdte,'-ko','k',8,2);
    oplot(ts,sh*msdts,'rd','none',15,2);
    oplot(td,sh*msdtd,'bs','none',15,2);
    yticks([1e-5 1e-4 1e-3 1e-2 ])
    hold off 
    %saveas(1,['fig06a_',num2str(loop),'.svg'],'svg')
    close(1)    
    xx          = dlmread(['./data/FBM_xx_a_',ext,'_np_50.dat']);
    yy          = dlmread(['./data/FBM_yy_a_',ext,'_np_50.dat']);
    [xi3,vacf3] = ea_vacf(dt,3,N,M,xx,yy,0,'lin');
    [xi5,vacf5] = ea_vacf(dt,5,N,M,xx,yy,0,'lin');
    [xi7,vacf7] = ea_vacf(dt,7,N,M,xx,yy,0,'lin');
    [tt,vv]     = vacf_fbm_theo(alpha);
    [rr,qq]     = vacf_fbm_err(5,theta,alpha);
    ax          = plotter([0 3 -0.5 1],'\xi','C(\xi)','xlin','ylin');        
    oplot(tt,vv,'-k','none',5,3);
    oplot(rr,qq,'-b','none',5,3);
    oplot(xi3,vacf3,'ks','none',10,2);
    oplot(xi5,vacf5,'bd','none',10,2);
    oplot(xi7,vacf7,'ro','none',10,2);
    %saveas(1,['fig06b_',num2str(loop),'.svg'],'svg')
    close(1)
end

%%--> Fig. 7A
N  = 2000;
M  = 100;
dt = 0.125;
xx = dlmread('data/telomers_xx.dat');
yy = dlmread('data/telomers_yy.dat');
[tau,msdte] = eata_msd(dt,N,M,xx,yy,0,'log');
A           = [tau',msdte'];
msdtg       = msdte*0;
for i=1:M
    [taut,msdt] = ta_msd(dt,xx(i,1:N),yy(i,1:N),0,'log');
    msdtg       = msdtg+log(msdt);
    if (mod(i,M/5) == 0) A=[A,msdt']; end;
end
A = [A, exp(msdtg'/M)];
ax=plotter([1d-1 2d2 5e-4 2e-1],'\tau [s]','\langler^2(\tau)\rangle_t [\mum^2]','xlog','ylog');
oplot(A(:,1),A(:,2),'ro','none',12,2);
for i=3:7  
    oplot(A(:,1),A(:,i),'b','none',1,1);
end
oplot([0.15,5],6d-3*[0.15,5].^0.5,'--k','none',1,3);
oplot([7,50],3d-3*[7,50],'--k','none',1,3);    
text(0.63,1e-2,'\sim\tau^{0.5}','FontSize',20);
text(12.6,0.09,'\sim\tau','FontSize',20);
%saveas(1,'fig07a.svg','svg')
close(1)     

%%--> Fig. 7B
[xi3,vacf3] = ea_vacf(dt,3,N,M,xx,yy,0,'lin');
[xi5,vacf5] = ea_vacf(dt,5,N,M,xx,yy,0,'lin');
[xi7,vacf7] = ea_vacf(dt,7,N,M,xx,yy,0,'lin');
[tt,vv]     = vacf_fbm_theo(0.5);
ax=plotter([0 3.2 -0.35,1.0],'\xi','C(\xi)','xlin','ylin');
oplot(xi3,vacf3,'sr','none',12,2);
oplot(xi5,vacf5,'dk','none',12,2);
oplot(xi7,vacf7,'ob','none',12,2);
oplot(tt,vv,'-k','none',12,3);               
%saveas(1,'fig07b.svg','svg')
close(1)

%--> inset
chi = []; %symmetrized, normalized
for j=1:M
    [ddx,ddy] = get_increm(3,xx(j,:),yy(j,:),true);
    chi       = [chi,abs(ddx),abs(ddy)];
end
[yh,xh] = histcounts(chi,'BinWidth',0.15,'Normalization','pdf');
xh      = xh+(xh(2)-xh(1))/2;
xh      = xh(1:end-1);
xr      = (0:400)/100;         
ax      = plotter([0 4 5e-4 1],'|\chi|','p(|\chi|)','xlin','ylog');
oplot(xh,yh,'ro','none',12,2);
oplot(xr,2*exp(-xr.^2/2)/sqrt(2*pi),'--k','none',2,3);         
%saveas(1,'fig07bi.svg','svg')
close(1)

%%--> Fig. 7C
[f,psde] = ea_psd(dt,N,M,xx,yy,0,'log');
A        = [f',psde'];
psdg     = psde*0;
for i=1:M
    [ff,ddd] = psd(dt,xx(i,1:N),yy(i,1:N),0,'log');
    psdg       = psdg+log(ddd);
    if (mod(i,M/5) == 0) A=[A,ddd']; end;
end
A  = [A, exp(psdg'/M)];
ax = plotter([3e-3 8 2e-5 5e1],'f','S(f)','xlog','ylog');
oplot(A(:,1),A(:,2),'ro','none',12,2);
for i=3:7 
    oplot(A(:,1),A(:,i),'-b','none',1,1);
end
oplot([1e-2,1],3d-5./[1e-2,1].^1.5,'--k','none',1,3);
text(0.03,4e-4,'\sim1/f^{1.5}','FontSize',20);
yticks([1e-4 1e-3 1e-2 1e-1 1e0 1e1]);
%saveas(1,'fig07c.svg','svg')
close(1)

%--> inset
for j=1:M
    [ddx,ddy] = get_increm(3,xx(j,:),yy(j,:),false);
    xx(j,:)=xx(j,:)/sqrt(mean(ddx.^2));
    yy(j,:)=yy(j,:)/sqrt(mean(ddy.^2));    
end
[fT,gam] = cov_gamma(dt,N,M,xx,yy,0,'log');
ax       = plotter([1 1e3 0.8 2.5],'fT','\gamma','xlog','ylin');
oplot(fT,gam,'ro','none',12,2);
oplot([1,1d3],[1,1],'--k','none',2,3);         
%saveas(1,'fig07ci.svg','svg')
close(1)

%%--> from here: helper routines for plotting

function erg=plotter(rang,xla,yla,stx,sty)
f=figure('Color','white','Units','inches','Position',[0,0,7.5,5]);
set(f,'renderer','Painters')
plot(rang(1),rang(3),'k.');
hold on
axis(rang);
ax = gca; % current axes
if (stx == 'xlog') ax.XScale = 'log'; end;
if (sty == 'ylog') ax.YScale = 'log'; end;
ax.XLabel.String = xla; % x-label
ax.YLabel.String = yla; % y-label
ax.LineWidth     = 2;   % axes line thickness
ax.TickLength    = [0.02 0.035]; % tick length
ax.FontSize      = 18;  
ax.FontName      = 'Arial';
erg=ax;
end

function res=oplot(x,y,psym,fil,siz,lnw)
plot(x,y,psym,'MarkerSize',siz,'LineWidth',lnw,'MarkerFaceColor',fil);
res=1;
end


