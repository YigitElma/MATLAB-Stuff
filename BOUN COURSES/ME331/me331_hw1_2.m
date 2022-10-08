clc 
clear all
close all

OS3=@(x) exp(-(pi*x)/sqrt(1-x^2))*100-10; %for overshoot value 10 this func is 0
ksi3=fzero(OS3,0) %finding zero of function
Ts3=@(x) 4/ksi3/x-4;  %for settling time 4 this func=0
w=fzero(Ts3,1)


G=tf([w^2], [1 2*ksi3*w w^2]);  %to check values
step(G,20)
stepinfo(G)
