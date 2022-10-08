%ME303 LAB SESSION
%Euler's Method

clc
clear all
close all

% The interval at which y(t) is desired to be obtained is [t0,tf] 
% where t0=0, tf=3
% y'=f(t,y)=(t-y)/2;
% The initial condition: y(0)=y0=1
% The exact solution=y(t)=3*e^(-t/2)+t-2;
% Find y(tf) at h1=1 and h2=h1/2
% Find the Final Global Error (FGE) at t=tf. It is expected 
% that FGE=Ch 
% where C is a constant. 
% If the step size in Euler’s method is reduced by a factor 
% of 1/2, we can expect that the overall F.G.E. will be reduced 
% by a factor of 1/2.
% Plot y(t) with the approximated solutions on [t0,tf]

t0=0;
tf=3;
y0=1;
h1=1;
h2=h1/4;
M1=(tf-t0)/h1;
M2=(tf-t0)/h2;

f=@(t,y) (t-y)/2;
y=@(t) 3*exp(-t/2)+t-2;

t1=t0:h1:tf;
t2=t0:h2:tf;


y1=zeros(1,M1+1);
y2=zeros(1,M2+1);

y0=1;
y1(1)=y0;
y2(1)=y0;

for j=1:M2
    if j<=M1
        y1(j+1)=y1(j)+h1*f(t1(j),y1(j));
        y2(j+1)=y2(j)+h2*f(t2(j),y2(j));
    else
        y2(j+1)=y2(j)+h2*f(t2(j),y2(j));
    end
end

y1_f=y1(end);
y2_f=y2(end);
y_f=y(tf);

FGE_1=abs(y_f-y1_f)
FGE_2=abs(y_f-y2_f)
%ratio
t=t0:0.001:tf;
figure
plot(t1,y1,'r--')
hold on
plot(t2,y2,'g--')
hold on
plot(t,y(t),'k','Linewidth',3)
xlabel('t')
ylabel('y')
legend('y1','y2','y')
grid on
