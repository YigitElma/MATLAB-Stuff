% ME303  Fall 2020
% Data linearization Example 2

%finding f(x)=1/(Ax+B) type of fit to a given set of data points
%main idea is same for all curve fits, try to express function as a
%polynomial of x and make necessary change of variables.

clc
clear all
close all


x=[-1 0 1 2 3];
y=[6.62 3.94 2.17 1.35 0.89];

N=length(x);
X=x;
Y=1./y;

m11=X*X';
m12=N*mean(X);
m21=m12;
m22=N;

l1=Y*X';
l2=N*mean(Y);

M=[m11 m12
    m21 m22];
L=[l1;l2];

U=M\L;
A=U(1);
B=U(2);

f=@(x) 1./(A*x+B);

a=-1.1:0.01:4;
plot(a,f(a),x,y,'*')
legend('f(x)=1/(Ax+B)','Data Points')
xlabel('x')
ylabel('y')
title('Least-Squares Curve with the Data Points')
grid on
