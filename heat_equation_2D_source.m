%Heat Equation 2D

clearvars; close all; clc;

dt=0.01; % time step
dx=0.01; % x-direction space step
dy=dx; % y-direction space step same as dx
alpha=12.5e-6; % alpha=k/rho*c
k=28; % W/mK
e_dot=2000000; % W/m^3
t_steps=1000000;  %Number of time steps
L=1;  %Length of bar 100 cm
sp_steps=floor(L/dx);
s=alpha*dt/dx^2;
 
p=100; % plot every p time steps

T_old=zeros(sp_steps,sp_steps);
T_new=zeros(sp_steps,sp_steps);
% posx=zeros(1,sp_steps);
% posy=zeros(1,sp_steps);

% Set boundaries
for m=1:sp_steps
    T_old(m,1)=0;
    T_old(m,sp_steps)=0;
    T_new(m,1)=0; % for consistent ouput
    T_new(m,sp_steps)=0; % for consistent ouput
end    
for n=1:sp_steps
    T_old(1,n)=0;
    T_old(sp_steps,n)=0;
    T_new(1,n)=0; % for consistent ouput
    T_new(sp_steps,n)=0; % for consistent ouput
end   

% Initial condition
for m=2:sp_steps-1
    for n=2:sp_steps-1
        T_old(m,n)=0;
    end
end    

for i=1:t_steps
    
    for m=2:sp_steps-1
        for n=2:sp_steps-1
            if(m<=floor(sp_steps/2) && n<=floor(sp_steps/2))
                T_new(m,n)=s*(T_old(m+1,n)+T_old(m-1,n)+T_old(m,n+1)+T_old(m,n-1))+(1-4*s)*T_old(m,n)+e_dot*alpha*dt/k;
            else   
                T_new(m,n)=s*(T_old(m+1,n)+T_old(m-1,n)+T_old(m,n+1)+T_old(m,n-1))+(1-4*s)*T_old(m,n);
            end
%             if m<=55 && m>=45 && n<=55 && n>=45
%                 T_new(m,n)=s*(T_old(m+1,n)+T_old(m-1,n)+T_old(m,n+1)+T_old(m,n-1))+(1-4*s)*T_old(m,n)+e_dot*alpha*dt/k;
%             else   
%                 T_new(m,n)=s*(T_old(m+1,n)+T_old(m-1,n)+T_old(m,n+1)+T_old(m,n-1))+(1-4*s)*T_old(m,n);
%             end
            
%            T_new(m,n)=s*(T_old(m+1,n)+T_old(m-1,n)+T_old(m,n+1)+T_old(m,n-1))+(1-4*s)*T_old(m,n)+e_dot*alpha*dt/k;
%             posx(m)=m*dx;
%             posy(n)=n*dy;
        end
    end
      if(mod(i,p)==0)
         %image(T_new);
         surf(T_new,'edgecolor','none');  
          view([0 90])
          view([33 44])
          axis([1 100 1 100 -10 1200]);
          xlabel('x-direction');
          ylabel('y-direction');
          zlabel('Temperature');
         pause(0.000000001);
      end  

%     for m=2:sp_steps-1
%         for n=2:sp_steps-1
%             T_old(m,n)=T_new(m,n);
%         end
%     end
T_old=T_new;
    
end    
    