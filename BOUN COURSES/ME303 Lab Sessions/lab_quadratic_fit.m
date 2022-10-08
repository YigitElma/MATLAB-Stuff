%me303 fall 2020
%Quadratic Fit

%This code is for finding appropriate quadratic fit to given set of data
%points. Calculation of coefficients of matrix M is standart. Only
%difference will be size of the matrix 3*3. When making linear fit this 
%matrix will be 2*2.  

clc
clear all
close all

X=[-3 0 2 4];
Y=[3 1 1 3];

scatter(X,Y,'r*')

N=length(X);
m11=sum(X.^4);
m12=sum(X.^3);
m13=sum(X.^2);
m21=m12;
m22=m13;
m23=sum(X);
m31=m22;
m32=m23;
m33=N;

L1=Y*X'.^2;  %L1=sum(Y.*X.^2) şeklinde de yazılabilir
L2=Y*X';
L3=sum(Y);

M=[m11 m12 m13;m21 m22 m23;m31 m32 m33];
L=[L1;L2;L3];

U=M\L;

A=U(1);
B=U(2);
C=U(3);

f=@(x) A*x.^2+B*x+C;
a=-4:0.01:5;
plot(a,f(a),X,Y,'r*')
legend('f(x)=Ax^{2}+Bx+C','Data Points') %^{...} yazdığında üste yazar
xlabel('x values')
ylabel('y values')
title('Least-Squares Quadratic Fit')
grid on