clc;
clear all;
close all;
format long

density_water = 1;
density_ball = 0.71;
r=15;
epsilon= 10^(-6);
delta= 10^(-6);
max=100;
a = 0;
b = 30;     %max height can be 2*15=30

V_sub1 = @(x) pi*x*x*(3*r-x)/3;  %submerged volume of ball can be calculated by this formula
V_ball=4/3*pi*r^3;    %classic formula for volume of a sphere
f=@(x) V_ball*density_ball-V_sub1(x)*density_water;  %mass=volume*density, we want f(x) to be 0 for the root


a_k=[];
b_k=[];
c_k=[];
yc_k=[];

for i=1:max
    ya=f(a);
    yb=f(b);
    c = b-yb*(b-a)/(yb-ya);
    yc=f(c);
    dx=min(abs(c-a),abs(c-b));
    
    a_k=[a_k,a];
    b_k=[b_k,b];
    c_k=[c_k,c];
    yc_k=[yc_k,yc];
    
    if dx<delta && abs(yc)<epsilon
        break
    end
    if yb*yc<0
        a=c;
    else
        b=c;
    end
end
height=c
V_sub=V_sub1(height)
ratio=V_sub1(c)/V_ball

k=0:1:i-1;
k=k';
a_k=a_k';
b_k=b_k';
c_k=c_k';
yc_k=yc_k';
Data=[k,a_k,b_k,c_k,yc_k];
VarNames={'k','a_k','b_k','c_k','yc_k'};
T=table(Data(:,1),Data(:,2),Data(:,3),...
Data(:,4),Data(:,5),'variablenames',VarNames)    

    
    