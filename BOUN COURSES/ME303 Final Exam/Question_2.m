%ME303 FINAL EXAM  04.01.2021
%QUESTION 2
%SOLVING DIFFERENTIAL EQUATIONS

%NOTE: Solutions are so close to each other, at first sight you may think
%there is only one graph but actually each one is there. If you really zoom
%in you can see them.

clear all
close all
clc

h=0.05;
t=0:h:2;

d_x=@(y)  y;  %after dividing equation into 2 first order DE's I get
d_y=@(t,x,y) 0.5*( 5*y+3*x+45*exp(2*t) );   %these functions

N=length(t);
x=zeros(1,N);
y=zeros(1,N);

x(1)=2;  %These are the initial conditions
y(1)=1;

for i=1:N-1  %this part classic Runge-Kutta formuation for solving DE'
    f1_x = d_x( y(i) );
    f1_y = d_y( t(i) , x(i) , y(i) );
    
    f2_x = d_x( y(i)+h/2*f1_y );
    f2_y = d_y( t(i)+h/2 , x(i)+h/2*f1_x , y(i)+h/2*f1_y );
    
    f3_x = d_x( y(i)+h/2*f2_y );
    f3_y = d_y( t(i)+h/2 , x(i)+h/2*f2_x , y(i)+h/2*f2_y );
    
    f4_x = d_x( y(i)+h*f3_y );
    f4_y = d_y( t(i)+h , x(i)+h*f3_x , y(i)+h*f3_y );
    
    x(i+1) = x(i) + h/6*( f1_x + 2*f2_x + 2*f3_x + f4_x );
    y(i+1) = y(i) + h/6*( f1_y + 2*f2_y + 2*f3_y + f4_y );
end

analytic=@(t) 4*exp(-t/2) +7*exp(3*t) -9*exp(2*t);  %this is the analytical solution
t_continuous=0:0.000001:2;            %of the equation. I defined t_continuous
                                      %to plot "analytic"
tspan=[0 2];  %for ODE45 we need this vector
initial=[2 1];  %for ODE45 initial conditions vector
[t_ode X_ode] = ode45( @odefun, tspan, initial );  %this part is for ODE45 solution

%This part is for plotting every solution I found for this equation
plot( t, x, 'r')
grid on 
hold on
plot( t_continuous, analytic(t_continuous), 'b', 'linewidth', 2)
hold on
plot(t_ode, X_ode(:,1), 'k')
legend('Runge-Kutta Solution','Analytic Solution','ODE45 Solution')
title('Solutions')
xlabel('t (time)')
ylabel('x(t)')

%I created odefun function for ODE45
function X_ode = odefun(t,x)
    X_ode(1,1)= x(2);   %these are the same DE's that I used earlier
    X_ode(2,1)= 0.5 * ( 5*x(2) + 3*x(1) + 45*exp(2*t) );
end

    