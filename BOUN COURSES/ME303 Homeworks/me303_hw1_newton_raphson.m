clc;
clear all,
close all;
format long

epsilon=10^(-8);
delta=10^(-8);
max=100;
p0=0;
d_distance=@(x) 4*x^3-2*x-6; %distance=(3-x)^2+(1-x^2)^2, to find min distance we need d_distance=0
d_d_distance=@(x) 12*x^2-2;  %to find root o d_distance we need d_d_distance for newtons-raphson method

 for i=1:max
     p1=p0-d_distance(p0)/d_d_distance(p0);
     if (abs(p1-p0)<delta || (2*abs(p1-p0)/(abs(p1)+abs(p0)))<delta) && d_distance(p1)<epsilon
         break
     end
     p0=p1;
 end
 x_p=p1;
 y_p=x_p^2;
 p_closest=[x_p y_p]
 dist=sqrt((3-x_p)^2+(1-x_p^2)^2)
    