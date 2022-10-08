clear all; 
clc;
close all;

a=0.1;
b=0.9;
delta=10^(-10);
epsilon=delta;
max=100;

r=1;
V=4/3*pi*r^3;
V1=@(x) pi*(r-x).^2.*(2*r+x)/3;
V2=@(x) V-V1(x);
f=@(x) V2(x)-3*V1(x);

x=0:0.1:r;
plot(x,f(x))
grid on
