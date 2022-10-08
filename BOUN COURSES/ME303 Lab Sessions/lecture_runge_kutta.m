% ME303 W, December 30, 2020
% Second order ODE solution using RK4 method

close all
clear all
clc

t0=0;
tf=2;
h=0.1;
M=(tf-t0)/h;

x=zeros(1,M+1);
y=zeros(1,M+1);
t=t0:h:tf;

f=@(t,x,y) y;
g=@(t,x,y) 0.5*(5*y+3*x+45*exp(2*t));

x(1)=2;
y(1)=1;

for i=1:M
    f1=f(t(i),x(i),y(i));
    g1=g(t(i),x(i),y(i));
    
    f2=f(t(i)+h/2, x(i)+h*0.5*f1, y(i)+h*0.5*g1);
    g2=g(t(i)+h/2, x(i)+h*0.5*f1, y(i)+h*0.5*g1);
    
    f3=f(t(i)+h/2, x(i)+h*0.5*f2, y(i)+h*0.5*g2);
    g3=g(t(i)+h/2, x(i)+h*0.5*f2, y(i)+h*0.5*g2);
    
    f4=f(t(i)+h, x(i)+h*f3, y(i)+h*g3);
    g4=g(t(i)+h, x(i)+h*f3, y(i)+h*g3);
    
    x(i+1)=x(i)+(h/6)*(f1 +2*f2 +2*f3 +f4);
    y(i+1)=y(i)+(h/6)*(g1 +2*g2 +2*g3 +g4);
end

f_real=@(k) 4*exp(-k/2)+7*exp(3*k)-9*exp(2*k);

figure(1)
plot(t,x,'r')
hold on
p=0:0.001:2;
plot(p,f_real(p),'b')
title('RK4')
legend('RK4 Approx.','analytic sol.')
grid on

tSpan = [0 2];
x0 = [2 1];
[t1 x1] = ode45(@(t,x) odefun(t,x), tSpan, x0);

figure(2)
plot(t1,x1(:,1));
title('ode45 Solution')
grid on

function dx_dt=odefun(t,x);
dx_dt=zeros(2,1);
dx_dt(1)=x(2);
dx_dt(2)= 0.5*(5*x(2)+3*x(1)+45*exp(2*t));
end



