%ME 303 LAB SESSION
%Euler Method
%Let y' = x^3 + y^2 initial values problem with y(0)=0.5, find y(2)=?
%We will use Taylor Series 

clc
clear all
close all

h=0.1;
x=0:h:2;
y=zeros(size(x));
y(1)=0.5;

for i=1:length(x)-1
    f=x(i)^3+y(i)^2;
    y(i+1)=y(i)+f*h;
end

y(end)
%%
% ME 303 Euler Method 22 December, 2020
% y' = 0.1*y initial value problem
% with y(0) = 1000, find y(5)=?
% Remember Taylar Series:
% f(x+h) = f(x) + f'(x)*h + f''(x)*h^2/2! + ...

clc
clear all
close all

h = 1/60; 
x = 0:h:5;
y = zeros(size(x));
y(1) = 1000;

for i=1:length(x)-1
    f = 0.1*y(i);
    y(i+1) = y(i) + f*h;
end

y(end)






