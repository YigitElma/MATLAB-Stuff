%ME303 FINAL EXAM  04.01.2021
%QUESTION 3
%CURVE FITTING

clear all
close all
clc

%These are given in the question
x = 1:6;  
y = [ 28 34 36 38 39 39 ];
N = length(x);

%after a bit hand calculation I made these variable changes
Y = x./y;
X = x;

%These are classic curve fitting calculations
m11 = X*X';
m12 = sum(X) ;
m21 = m12 ;
m22 = N ;

r1 = X*Y' ;
r2 = sum(Y) ;

M = [ m11 m12; m21 m22 ];
R = [ r1; r2 ];

%O is the coefficient matrix for A and B 
O = M\R;
B = O(1);
A = O(2);

%to plot my result I defined this function handle
f = @(x) x./(B*x+A);
t = 1:0.001:6;

%this part is for plotting
plot(x,y,'r*',t,f(t),'b')
grid on
legend('Data Points','y=x/(A+Bx)')
xlabel('x')
ylabel('y')

