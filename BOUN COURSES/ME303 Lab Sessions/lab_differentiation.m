%ME 303 LAB SESSION
%

clc
clear all
close all

f=@(x) 60*x.^45-32*x.^33+233*x.^5-47*x.^2-77;
x=1/sqrt(3);

tol=10^-9;
eps=10^-9;
maxI=15;

h=1;
H(1)=h;
D(1)=(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
E(1)=0;
R(1)=0;

for n=1:2
    h=h/10;
    H(n+1)=h;
    D(n+1)=(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
    E(n+1)=abs(D(n+1)-D(n));
    R(n+1)=E(n+1)/(abs(D(n+1))+eps);
end

while E(n)>E(n+1) && R(n)>tol && n<maxI
    h=h/10;
    H(n+2)=h;
    D(n+2)=(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
    E(n+2)=abs(D(n+2)-D(n+1));
    R(n+2)=E(n+2)/(abs(D(n+2))+eps);
    n=n+1;
end

k=length(D)-1
D_approx=D(k)
L=[H' D' E' R']


