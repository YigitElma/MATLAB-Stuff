%me303 fall 2020
%Data Linearization

%for given data points finding exponential fit of the form f(x)=Cexp(Ax)
%since we know how to make a linear fit, we use the same loqic. Main idea
%is to make some convenient variable changes so that exponential fit will
%be linear fit.
%Hint: leave the x alone

clc
clear all
close all

x=[1 2 3 4 5];
y=[0.6 1.9 4.3 7.6 12.6];
N=length(x);

X=x;
Y=log(y);

m11=sum(X.^2);
m12=sum(X);
m21=m12;
m22=N;
L1=Y*X';
L2=sum(Y);

M=[m11 m12;m21 m22];
L=[L1;L2];

U=M\L;  %M^(-1)*L nin farklı bir gösterim şekli, L/M olarak çalışmaz

A=U(1);
B=U(2);
C=exp(B);
f=@(x) C*exp(A*x);

k=0:0.01:6;
plot(k,f(k),'r',x,y,'b*')
legend('f(x)=Ce^{Ax}','Data Points') %^{...} yazdığında üste yazar
xlabel('x values')
ylabel('y values')
title('Least-Squares Exponential Fit')
grid on
