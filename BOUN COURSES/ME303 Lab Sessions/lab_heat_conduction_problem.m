%ME303 LAB SESSION
%1 dimensional transient heat conductiton problem

clc
clear all 
close all

%k/rho/c=1
k=1;
rho=1;
c=1;

dt=0.001;
dx=0.05;
xl=0;
xr=1;
x=xl:dx:xr;
N=length(x);
T=zeros(N,1);
T(1)=20;
a=dt/dx^2;   %a = k*dt/(rho*c*dx^2)
Tp=T;
A=zeros(N,N);

for i=1:N
    if i==1
        A(i,i)=1;
        A(i,i+1)=0;
    elseif i==N
        A(i,i-1)=2*a;
        A(i,i)=-2*a+1;
    else 
        A(i,i-1)=a;
        A(i,i)=-2*a+1;
        A(i,i+1)=a;
    end
end

figure 
hold on
xlabel('x (meters)','Interpreter','tex','FontSize',16)
ylabel('Temperature (C)','Interpreter','tex','FontSize',16)

for i=1:10000
    T=A*Tp;
    if mod(i,250)==0
        plot(x,T)
        txt=text(x((N+1)/2),T((N+1)/2),['t=' num2str(i*dt)],...
            'VerticalAlignment','middle','HorizontalAlignment',...
            'center','Interpreter','latex','FontSize',10);
    end
    if max(abs(T-Tp))<10^-4
        tss=i*dt;
        break
    end
    Tp=T;
end
       




