function [sx,sy] = make_switchwalk(dt,N,M,xx,yy,k1,k2,fact)
%------------------------------------------------------
% add switching of diffusion coefficient in trajectory ensemble 
%------------------------------------------------------
% dt      time increment / frame time
% N       length of trajectories 
% M       ensemble size (number of trajectories
% xx      array of x coordinates
% yy      array of y coordinates
% k1      forward switch rate (to slower diffusion) 
% k2      backward switch rate (to faster diffusion) 
% fact    factor for faster diffusion
%------------------------------------------------------

[x,y] = enscheck(N,M,xx,yy);
pon   = dt*k1;     %-- switch probability 
poff  = dt*k2;     %-- switch back probability 
sx    = x;
sy    = y;

if (fact <1)
    error('> make_switchwalk stopped, wrong parameter 8 <');
end

for j=1:M
    dx = x(j,2:N)-x(j,1:N-1);
    dy = y(j,2:N)-y(j,1:N-1);    
    ask   = (random('Uniform',0,1) < k1/(k1+k2));  %-- initial state
    hop   = random('Uniform',0,1,N,1);
    for i=1:N-1 
        if (ask)  %-- "on" => slow
            sx(j,i+1) = sx(j,i)+dx(i);
            sy(j,i+1) = sy(j,i)+dy(i);        
            if (hop(i) < poff) ask=false; end; 
        else      %-- "off" => fast
            sx(j,i+1) = sx(j,i)+dx(i)*fact;
            sy(j,i+1) = sy(j,i)+dy(i)*fact;        
            if (hop(i) < pon) ask=true; end; 
        end       
    end
end
