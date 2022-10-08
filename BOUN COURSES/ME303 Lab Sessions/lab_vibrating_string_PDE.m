clear all
close all
clc

%ME303 LAB SESSION
%SOLUTION OF PDEs

% Vibration String
% PDE: u_tt(x,t)=(c^2)*u_xx(x,t) for a<=x<=b and t0<=t<=tf
% BCs: u(x1,t)=0 and u(xf,t)=0
% ICs: u(xi,0)=f(x) u_t(x,0)=v(x)

L=1.4;   %wire length in meter
T=1000; %tensin in newtons
m=0.11;  %linear mass density (kg/m)
c=sqrt(T/m);  %wave speed (m/s)

A=0.05;  %magnitude of the wave (m)
f=@(x) A*sin(pi*x/L).*x;
v=@(x) 0;

n=11;
dx=L/(n-1);
x=0:dx:L;

tf=0.1; %total time span
dt=0.0004;
t=0:dt:tf;
M=tf/dt+1;

c_disc=dx/dt;
r=c/c_disc;

%Displacement matrix u(x,t)
%BCs: u(0,t)=0, u_t(0,t)=0
u=zeros(n,M);

for i=2:n-1
    u_tt=(c^2)*(f(x(i+1))-2*f(x(i))+f(x(i-1)))/dx^2;
    u_t=v(x(i));
    u(i,1)=f(x(i));
    u(i,2)=u(i,1)+u_t*dt+u_tt/2*dt^2;
end

for j=3:M
    for i=2:n-1
        u(i,j)=(2-2*r^2)*u(i,j-1)+r^2*(u(i+1,j-1)+u(i-1,j-1)) - u(i,j-2);
    end
end

% Animating string

figure
axis([0 L -2*A 2*A])
xlabel('x(m)')
ylabel('u(m)')
grid on
ax=gca;
ax.NextPlot='replaceChildren';

for j=1:M
    plot(x,u(:,j),'linewidth',1.5)
    F(j)=getframe;
end
