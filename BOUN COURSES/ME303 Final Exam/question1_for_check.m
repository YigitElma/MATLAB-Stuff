%ME303 FINAL EXAM  04.01.2021
%QUESTION 1
%SOLVING DIFFERENTIAL EQUATIONS


%NOTE: I used this code to check my calculations in Question 1.

clear all
close all
clc

k=2;
m=4;
c=0.4;
t=0:2;
h=1;

d_y=@(z) z;
d_z=@(y,z) 9.81-k/m*y-c/m*z;

y=zeros(1,3);
z=zeros(1,3);

y(1)=0;
z(1)=0;

for i=1:2
    display('for i')  %to make distinction between i=1 and i=2
    
    y1=d_y( z(i) )  %I didn't put any ; to compare my results
    z1=d_z( y(i),z(i) )
    
    y2=d_y( z(i)+h/2*z1 )
    z2=d_z( y(i)+h/2*y1, z(i)+h/2*z1 )
    
    y3=d_y( z(i)+h/2*z2 )
    z3=d_z( y(i)+h/2*y2, z(i)+h/2*z2 )
    
    y4=d_y( z(i)+h*z3 )
    z4=d_z( y(i)+h*y3, z(i)+h*z3 )
    
    y(i+1)=y(i)+h/6*(y1+2*y2+2*y3+y4);
    z(i+1)=z(i)+h/6*(z1+2*z2+2*z3+z4);
end
y  %to check y and z values at the end
z
plot(t,y)  %this is not really necesssary
    
    
    
    

