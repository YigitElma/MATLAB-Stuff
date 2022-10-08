%Yiğit Günsür ELMACIOĞLU
%2017405120
%18 JANUARY 2021

%ME303 TERM PROJECT 2020 FALL
%QUESTION 2
%RESTIRCTED 3-BODY PROBLEM

clc
clear all
close all

%Given numerical values
m1 = 0.012277471;
m2 = 1-m1;
T  = 17.06521656015796;
disp_T=2;         %for much longer display times (for example disp_T*T = 2*T)
                  %I indicated this on the graphs
c_rk=6000;      %for much higher precisions I'll change this value
h_rk = T/c_rk;     %for Runge-Kutta
c_e=24000;       %for much higher precisions I'll change this value
h_e  = T/c_e;      %for Euler

t_rk = 0:h_rk:disp_T*T;  %for Runge-Kutta
t_e  = 0:h_e:disp_T*T;    %for Euler

N=length(t_rk);
M=length(t_e);

%I'll store Runge-Kutta and Euler solutions in different matrices
%for Runge-Kutta Solution
x_rk=zeros(1,N);
y_rk=zeros(1,N);
z_rk=zeros(1,N);
w_rk=zeros(1,N);
%for Euler Solution
x_e=zeros(1,M);
y_e=zeros(1,M);
z_e=zeros(1,M);
w_e=zeros(1,M);

%Differential equations after dividing them into 4 first order differential equation
d_X=@(z) z;
d_Y=@(w) w;
d_Z=@(x,y,w) x + 2*w - m2*(x+m1) / (((x+m1)^2 + y^2)^(3/2)) - m1*(x-m2) / (((x-m2)^2 + y^2)^(3/2));
d_W=@(x,y,z) y - 2*z - m2*y / (((x+m1)^2 + y^2)^(3/2)) - m1*y / (((x-m2)^2 + y^2)^(3/2));

%Initial conditions for Runge-Kutta solution
x_rk(1)=0.994;
y_rk(1)=0;
z_rk(1)=0;
w_rk(1)=-2.0015851063790825;

%Initial conditions for Euler solution
x_e(1)=0.994;
y_e(1)=0;
z_e(1)=0;
w_e(1)=-2.0015851063790825;

%Runge-Kutta Method solution
for i=1:N-1    
    f1_X=d_X( z_rk(i) );
    f1_Y=d_Y( w_rk(i) );
    f1_Z=d_Z( x_rk(i), y_rk(i), w_rk(i) );
    f1_W=d_W( x_rk(i), y_rk(i), z_rk(i) ); 
    
    f2_X=d_X( z_rk(i)+h_rk/2*f1_Z );
    f2_Y=d_Y( w_rk(i)+h_rk/2*f1_W );
    f2_Z=d_Z( x_rk(i)+h_rk/2*f1_X, y_rk(i)+h_rk/2*f1_Y, w_rk(i)+h_rk/2*f1_W );
    f2_W=d_W( x_rk(i)+h_rk/2*f1_X, y_rk(i)+h_rk/2*f1_Y, z_rk(i)+h_rk/2*f1_Z );
   
    f3_X=d_X( z_rk(i)+h_rk/2*f2_Z );
    f3_Y=d_Y( w_rk(i)+h_rk/2*f2_W );
    f3_Z=d_Z( x_rk(i)+h_rk/2*f2_X, y_rk(i)+h_rk/2*f2_Y, w_rk(i)+h_rk/2*f2_W );
    f3_W=d_W( x_rk(i)+h_rk/2*f2_X, y_rk(i)+h_rk/2*f2_Y, z_rk(i)+h_rk/2*f2_Z );
    
    f4_X=d_X( z_rk(i)+h_rk*f3_Z );
    f4_Y=d_Y( w_rk(i)+h_rk*f3_W );
    f4_Z=d_Z( x_rk(i)+h_rk*f3_X, y_rk(i)+h_rk*f3_Y, w_rk(i)+h_rk*f3_W );
    f4_W=d_W( x_rk(i)+h_rk*f3_X, y_rk(i)+h_rk*f3_Y, z_rk(i)+h_rk*f3_Z );
    
    x_rk(i+1)=x_rk(i) + h_rk/6 * (f1_X + 2*f2_X + 2*f3_X + f4_X);
    y_rk(i+1)=y_rk(i) + h_rk/6 * (f1_Y + 2*f2_Y + 2*f3_Y + f4_Y);
    z_rk(i+1)=z_rk(i) + h_rk/6 * (f1_Z + 2*f2_Z + 2*f3_Z + f4_Z);
    w_rk(i+1)=w_rk(i) + h_rk/6 * (f1_W + 2*f2_W + 2*f3_W + f4_W);
end

%Euler's method solution 
for j=1:M-1
    d_X_e=d_X( z_e(j) );
    d_Y_e=d_Y( w_e(j) );
    d_Z_e=d_Z( x_e(j), y_e(j), w_e(j) );
    d_W_e=d_W( x_e(j), y_e(j), z_e(j) );
    
    x_e(j+1)=x_e(j) + h_e*d_X_e;
    y_e(j+1)=y_e(j) + h_e*d_Y_e;
    z_e(j+1)=z_e(j) + h_e*d_Z_e;
    w_e(j+1)=w_e(j) + h_e*d_W_e;
end

%Graph of Runge-Kutta Solution
plot(x_rk, y_rk)
hold on
plot(-m1, 0, 'b*','linewidth',5)  %this represent the Earth
hold on
plot(m2, 0, 'k.')                 %this represents the Moon
%name_rk is only for visual purposes, to show each parameter of the graph on the title
%without typing for each case, for example when I change h it's automaticly
%adapts to new value
name_rk=sprintf('Runge-Kutta Solution with h=T/%d and for %d*T', c_rk, disp_T);
title(name_rk)
xlabel('x axis (y1)')
ylabel('y axis (y2)')

%Graph of Euler Solution
figure
plot(x_e,y_e)
hold on
plot(-m1, 0, 'b*','linewidth',5)
hold on
plot(m2, 0, 'k.')
%only for visual purposes, to show each parameter of the graph on the title
name_e=sprintf('Euler Solution with h=T/%d and for %d*T', c_e, disp_T);
title(name_e)
xlabel('x axis (y1)')
ylabel('y axis (y2)')

%Graph of ODE45 Solution
tspan=[0 disp_T*T];
[t,y]=ode45(@odefun,tspan,[0.994; 0; 0; -2.0015851063790825]);
figure 
plot(y(:,1),y(:,2))
hold on
plot(-m1, 0, 'b*','linewidth',5)
hold on
plot(m2, 0, 'k.')
%only for visual purposes, to show each parameter of the graph on the title
name_45=sprintf('ODE45 Solution for %d*T',disp_T);
title(name_45)
xlabel('x axis (y1)')
ylabel('y axis (y2)')

function dydt=odefun(t,y)
    m1 = 0.012277471;
    m2 = 1-m1;    
    
    dydt(1,1)=y(3);
    dydt(2,1)=y(4);
    dydt(3,1)=y(1) + 2*y(4) - m2*(y(1)+m1) / (((y(1)+m1)^2 + y(2)^2)^(3/2)) - m1*(y(1)-m2) / (((y(1)-m2)^2 + y(2)^2)^(3/2));
    dydt(4,1)=y(2) - 2*y(3) - m2*y(2) / (((y(1)+m1)^2 + y(2)^2)^(3/2)) - m1*y(2) / (((y(1)-m2)^2 + y(2)^2)^(3/2));
end
