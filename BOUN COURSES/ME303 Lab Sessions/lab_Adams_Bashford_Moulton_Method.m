% ME 303 5 January 2021 Lab Session
% Predictor Corrector Method
% Adams-Bashford-Moulton Method

% y' = -5y + 2t + 30
% y(0) = 1, y(5)=? with dt=0.001 increments
% y1 = 1 y2= 1.02500 y3=1.04988 y4=1.07464 

clc
clear all
close all

f = @(y,t) -5*y + 2*t + 30;
dt = 0.001;
ti = 0;
tf = 5;
N = (tf-ti)/dt + 1;
y = zeros(1,N);
t = ti:dt:tf;
y(1:4) = [1 1.025 1.04988 1.07464 ];

for k = 4:N-1
    fk3 = f(y(k-3),t(k-3));
    fk2 = f(y(k-2),t(k-2));
    fk1 = f(y(k-1),t(k-1));
    fk = f(y(k),t(k));
    
    p = y(k) + dt/24*(-9*fk3 + 37*fk2 - 59*fk1 + 55*fk);
    
    fkp = f(p,t(k+1));
    y(k+1) = y(k) + dt/24*(fk2 - 5*fk1 + 19*fk + 9*fkp );
end

y5 = y(end)
