clear all
close all
clc
format long
V=@(x) 4*x*(8-x)*(12-x)-500; %4*x*(8-x)*(12-x) is the volume,so to find volume=500 we want root of volume-500
eps=10^(-8);
p0=4;  %half of final width=(8-x)is smaller than x for x>4
p1=8;  %2x cannot be greater than 16 for a real box
for i=1:100
    p2 = p1 - V(p1)/((V(p1) - V(p0)) / (p1-p0)); %formula of secant method
    p0 = p1; %subtitude p0 with p1 for next iteration
    p1 = p2; %subtitude p1 with p2 for next iteration
    y = V(p2);  
    if abs(y)<eps
        break
    end   
end
root = p2
vol=V(root)+500 %function was volume-500 so we have to add 500