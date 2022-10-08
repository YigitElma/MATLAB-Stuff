% ME 303 Lab Session 9: Runge-Kutta Method

% y'=f(t,y)=(t-y)/2
% The interval at which y(t) is desired to be approximated is [t0,tf] where
% t0=0, tf=3
% The initial condition: y(0)=y0=1
% The exact solution=y(t)=3*e^(-t/2)+t-2;;
% Find y(tf) at h=1 using the Taylor's and Runge-Kutta methods
% Find the Final Global Error (FGE) at t=tf for each method.
% It is expected that FGE=Ch^4 for the Runge-Kutta method where C is a
% constant.
% If the step size in the Runge-Kutta method is reduced by a factor of 1/2, we can 
% expect that the overall F.G.E. will be reduced by a factor of 1/16.
% Plot y(t) with the approximated solutions from the two methods on [t0,tf]


clear all
close all
clc

N=4;
syms t y
f=(t-y)/2;
d(1)=f;

for k=1:N-1
    d(k+1)=diff(d(k),t)+f*diff(d(k),y);
end
d=matlabFunction(d);
f=matlabFunction(f);
% f = @(t,y) (t-y)/2;

t0=0;
tf=3;
y0=1;
h=1;
M=(tf-t0)/h;
t=t0:h:tf;
y_t=zeros(1,M+1);
y_rk=zeros(1,M+1);

y_t(1)=y0;
y_rk(1)=y0;

y = @(t) 3*exp(-t/2)+t-2;

C=zeros(1,N);
for i=1:N
    C(i)=h^(i)/factorial(i);
end

for j=1:M
    % Taylor Part
    y_t(j+1)=y_t(j)+C*d(t(j),y_t(j))';
    
    % Runge-Kutta Part
    f1_rk=f(t(j),y_rk(j));
    f2_rk=f(t(j)+h/2,y_rk(j)+h/2*f1_rk);
    f3_rk=f(t(j)+h/2,y_rk(j)+h/2*f2_rk);
    f4_rk=f(t(j)+h,y_rk(j)+h*f3_rk);
    
    y_rk(j+1)=y_rk(j)+h*(f1_rk+2*f2_rk+2*f3_rk+f4_rk)/6;
end

y_t_f=y_t(end);
y_rk_f=y_rk(end);
y_f=y(tf);
FGE_t=abs(y_f-y_t_f)
FGE_rk=abs(y_f-y_rk_f)

tt=t0:0.001:tf;

figure
plot(t,y_t,'r--','LineWidth',2)
hold on
plot(t,y_rk,'g--')
hold on
plot(tt,y(tt),'k','LineWidth',3)
xlabel('t')
ylabel('y(t)')
legend('y_t','y_r_k','y')
% legend([' h1 = ' num2str(h1)],[' h2 = ' num2str(h2)],'exact')


