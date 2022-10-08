%ME303 LAB SESSION

% y'=f(t,y)=(t-y)/2
% The interval at which y(t) is desired to be approximated is [t0,tf] where
% t0=0, tf=3
% The initial condition: y(0)=y0=1
% The exact solution=y(t)=3*e^(-t/2)+t-2;;
% Find y(tf) at h=1 using the Euler's and Heun's methods
% Find the Final Global Error (FGE) at t=tf for each method.
% It is expected that FGE=Ch^2 for the Heun's method where C is a constant.
% If the step size in Heun’s method is reduced by a factor of 1/2, we can 
% expect that the overall F.G.E. will be reduced by a factor of 1/4.
% Plot y(t) with the approximated solutions from the two methods on [t0,tf]

clc
clear all
close all

t0=0;
tf=3;
y0=1;
h=1;
M=(tf-t0)/h;

t=t0:h:tf;

f=@(t,y) (t-y)/2;
y=@(t) 3*exp(-t/2)+t-2;

y_e=zeros(1,M+1);
y_h=zeros(1,M+1);

y_e(1)=y0;
y_h(1)=y0;

for j=1:M
    %Euler part
    f1_e=f(t(j),y_e(j));
    y_e(j+1)=y_e(j)+h*f1_e;
    %Heun's part
    f1_h=f(t(j),y_h(j));
    f2_h=f(t(j+1),y_h(j)+h*f1_h);
    y_h(j+1)=y_h(j)+h/2*(f1_h+f2_h);
end

y_e_f=y_e(end)
y_h_f=y_h(end)
y_f=y(tf)
FGE_e=abs(y_f-y_e_f)
FGE_h=abs(y_f-y_h_f)

y_average=(y_e+y_h)/2;

tt=t0:0.001:tf;

figure
plot(t,y_e,'r--')
hold on
plot(t,y_h,'g--')
hold on
plot(tt,y(tt),'k')
hold on
plot(t,y_average)
grid on
xlabel('t')
ylabel('y')
legend('y_e','y_h','y','y_a_v_e')

    
    
    
    
    
    