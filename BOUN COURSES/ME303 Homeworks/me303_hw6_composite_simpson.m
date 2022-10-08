clc
clear all
close all

f=@(x) 2*pi*exp(-x)*sqrt(1+exp(-2*x));
a=0;
b=1;

% in these commented lines, I found value M in another code
% syms x                             %I use smybolic function for taking derivatives
% f=2*pi*exp(-x)*sqrt(1+exp(-2*x));  %orginal function
% g=diff(diff(diff(diff(f))));       %finding minimum M value for a given accuracy requires maximum value of fourth derivative of the function
% g=matlabFunction(g);   %returning back to function handle
% x=0:0.0001:1;  %since extremum points are at the boundaries I can
% plot(x,g(x))   %use the graph of the function to determine absisca of the max y
% g(0)     %I see that maximum occurs at x=0, I took the value g(0) and
           % put into −(b−a)*g(0)*h^4/180=ERROR, after some algebra I found
           % M is greater than 51,... so minimum M can be 52.

M=52;
h=(b-a)/(2*M);

x=a:h:b;
s=zeros(1,M);

for i=1:M        %these are Composite Simpson's rule parts, I didn't make any change
    x1=x(2*i-1);
    x2=x(2*i);
    x3=x(2*i+1);
    f1=f(x1);
    f2=f(x2);
    f3=f(x3);
    s(i)=(h/3)*(f1+4*f2+f3);
end

A_s=sum(s)      %A_s is our result of the integral