%ME303 LAB
%Numerical integration using trapezoidal rule

clc
clear all
close all

f=@(x) 2+sin(2*sqrt(x));
a=1;
b=6;
M=40;
h=(b-a)/M;
s=0;

for k=1:(M-1) %don't use end points because they will be multiplied by 1/2 of h
    x=a+h*k;
    s=s+f(x); %taking all middle points in a parantheses(these construct rectangles of heiht h)
end

T=h*(f(a)+f(b))/2+h*s  %end points are used to find area under triangles at the ends







