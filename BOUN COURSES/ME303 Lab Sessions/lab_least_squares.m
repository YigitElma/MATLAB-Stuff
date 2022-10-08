clc
close all
clear all

N = 50;

X = zeros(1,N);
Y = zeros(1,N);

for k=1:N
    X(k) = 0.1*k;
    Y(k) = X(k) + cos(k^(1/2));
end

% scatter (X,Y,'r*')

m11 = X*X';
m12 = N*mean(X);
m21 = m12;
m22 = N;

l1 = Y*X';
l2 = N*mean(Y);

M = [m11 m12 ; m21 m22];
L = [l1 ; l2];

U = M\L;

A = U(1);
B = U(2);

f = @(x) A*x + B;

x = 0:0.01:X(end);
plot(x,f(x),X,Y,'*')
grid on
