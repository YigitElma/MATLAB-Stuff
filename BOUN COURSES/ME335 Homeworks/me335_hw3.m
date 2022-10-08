clc
clear all
close all

%given transfer function
G = tf( [0.05116] , [ 1 0.3307 1.744 ] ); 
%I'm defining s as a laplace variable to use in transfer functions
s = tf('s');

% f = 150*G/( 1+150*G );
% step(f)
% %according to part c ii, important constants can be found as
% wn = 2.04; 
% ksi = 0.636;
% p0 = 1;   %this term is chosen arbitrarly
% a=0.05116;
% b=0.3307;
% c=1.744;
% %PID constants are found by pole assignment technique
% ki = wn^2 * p0 / a
% kp = (wn^2 + 2*ksi*wn*p0 - c)/a
% kd = ( 2*wn*ksi + p0 - b )/a
% 
% %according to part c iii, important constants can be found as
% wn2 = 1.45; 
% ksi2 = 0.636;
% p02 = 2;   %this term is chosen arbitrarly
% %PID constants are found by pole assignment technique
% ki2 = wn2^2 * p02 / a
% kp2 = (wn2^2 + 2*ksi2*wn2*p02 - c)/a
% kd2 = ( 2*wn2*ksi2 + p02 - b )/a
% 
% T = (kd*s + kp + ki/s)*G/ ( 1 + (kd*s + kp + ki/s)*G );
% step(T,50)
% stepinfo(T)
% 
% figure
% U = (kd2*s + kp2 + ki2/s)*G/ ( 1 + (kd2*s + kp2 + ki2/s)*G );
% step(U,50)
% stepinfo(U)
% 
% Y = ( 0.0765 + 0.05116*s )/( s^2 + 0.3307*s + 1.744 ) ;
% rlocus(Y)
kp=120;
ki=0;
kd=30.75;
controller=pid(kp,ki,kd);
step(feedback(controller*G,[1]))
figure
pzmap(feedback(controller*G,[1]))







