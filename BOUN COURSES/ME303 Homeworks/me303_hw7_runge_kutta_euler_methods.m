clc
close all
clear all

a=1.3;
b=-0.25;
c=-0.0001;
h=0.5;
x=0:h:20;
N=length(x);

y_e=zeros(1,N);    %I will store Euler method results in this one
y_rk=zeros(1,N);   %also Runge Kutta results stored here
y_e(1)=1;          %initial condition
y_rk(1)=1;         %initial condition

T_e(1)=0;      %this trapezoidal differentiation will use Euler method y values
for i=1:N-1
    if i>1     %since we know T(1) no need to itterate
        T_e(i)=T_e(i-1)+h/2*(y_e(i-1)+y_e(i));
    end        %calculation of T(i) must be before y_e(i+1) calculation because it uses T(i)
    y_e(i+1)=y_e(i)+h*( a*y_e(i) + b*y_e(i)^2 + c*y_e(i)*T_e(i) );   %I implemented given y' function in the code
end
plot(x,y_e,'b','linewidth',2)
grid on 
hold on

T_rk(1)=0;     %this trapezoidal differentiation will use Runge Kutta method y values
for j=1:N-1
    if j>1
    T_rk(j) = T_rk(j-1)+h/2*(y_rk(j-1)+y_rk(j));
    end
    
    f = @(t,y) a*y+b*y^2+c*y.*T_rk(j);   %given first derivative
    
    f1_rk=f(x(j),y_rk(j));               %these are Runge Kutta variables
    f2_rk=f(x(j)+h/2,y_rk(j)+h/2*f1_rk);
    f3_rk=f(x(j)+h/2,y_rk(j)+h/2*f2_rk);
    f4_rk=f(x(j)+h,y_rk(j)+h*f3_rk);
    y_rk(j+1)=y_rk(j)+h*(f1_rk+2*f2_rk+2*f3_rk+f4_rk)/6;
end
plot(x,y_rk,'r','linewidth',2)
legend('y_e','y_rk')
xlabel('x')
ylabel('y')



