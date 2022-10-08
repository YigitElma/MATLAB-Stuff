% -------------------------------------------------------------------------
% 
% Yiğit Günsür ELMACIOĞLU  -  2017405120
% 
% -------------------------------------------------------------------------

clc
clear
close all

Q = [ 0 6 12 18 24 30 36 ] ;
H = [ 47.5 46.2 42.5 36.2 26.2 15 0 ] ;
Power = [ 133 142 153 164 172 174 174 ] ;

g = 9.81 ;
rho = 1000 ;
Lpm_to_kgps = 1/60000 ;     % conversion factor from lpm to kg/s

efficiency = Q.*H*rho*g*Lpm_to_kgps./Power*100 ;

plot(Q,efficiency,'ro',Q,Power,'b^',Q,H,'kp')

eff = fit(Q',efficiency','poly2');
Head = fit(Q',H','poly2');
Power_f = fit(Q',Power','poly2');

q_range = 0:0.1:36 ;
hold on 
plot(eff,'r')
hold on 
plot(Head,'k')
hold on 
plot(Power_f,'b')
grid on
xlabel('Q, Flowrate (lpm)')
legend('\eta Data','Power Data','Head Data','\eta Performance Curve','Power Performance Curve','Head Performance Curve')


