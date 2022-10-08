% ME 303 5 January 2021 Lab Session
% Systems of Differental Equations
% Predator-Prey Model

% x'(t) = A*x(t) - B*x(t)*y(t)
% y'(t) = C*x(t)*y(t) - D*y(t)
% Given A = 2, B = 0.02, C = 0.0002, D = 0.8
% h = 0.1 in weeks x(0)=3000 rabbits y(0)=120 foxes
% solve in the interval [0, 5]

clc, clear all, close all

A = 2; B = 0.02; C = 0.0002; D = 0.8;
f = @(t,x,y) A*x - B*x*y;
g = @(t,x,y) C*x*y - D*y;

t0 = 0; tf = 5; h = 0.1;
M = (tf-t0)/h + 1;
t = zeros(1,M);
x = zeros(1,M);
y = zeros(1,M);
t = t0:h:tf;
x(1) = 3000;
y(1) = 120;

for j=1:M-1
    f1 = f(t(j),x(j),y(j));
    g1 = g(t(j),x(j),y(j));
    f2 = f(t(j)+h/2 , x(j)+h/2*f1 , y(j)+h/2*g1);
    g2 = g(t(j)+h/2 , x(j)+h/2*f1 , y(j)+h/2*g1);
    f3 = f(t(j)+h/2 , x(j)+h/2*f2 , y(j)+h/2*g2);
    g3 = g(t(j)+h/2 , x(j)+h/2*f2 , y(j)+h/2*g2);
    f4 = f(t(j)+h , x(j)+h*f3 , y(j)+h*g3);
    g4 = g(t(j)+h , x(j)+h*f3 , y(j)+h*g3);
    
    x(j+1) = x(j) + (f1+2*f2+2*f3+f4)*h/6;
    y(j+1) = y(j) + (g1+2*g2+2*g3+g4)*h/6;
    
end
R = [t' x' y'];
x_end = round(x(end))
y_end = round(y(end))

figure(1)
subplot 211
plot(R(:,1),R(:,2),'r') % rabbits
title('Rabbits vs. Time')
grid on
legend('Rabbits')
xlabel('time')
ylabel('x(t) vs. y(t)')
hold on
subplot 212
plot(R(:,1),R(:,3),'g') % foxes
xlabel('time')
ylabel('x(t) vs. y(t)')
legend('Foxes')
title('Foxes vs. Time')
grid on
