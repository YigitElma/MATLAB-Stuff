% -------------------------------------------------------------------------
%
% 09.11.2021
% Yiğit Günsür Elmacıoğlu
% 
% Inverted Double Pendulum Simulation
% 
% -------------------------------------------------------------------------
clc
clear 
close all

global m1;
global m2;
global l1;
global l2;
global g;

m1 = 3 ;
m2 = 1 ;
l1 = 3 ;
l2 = 10 ;
g = 9.81 ;
dt = 0.001 ;
t_final = 15 ;

% initial conditions
teta1 = 10*pi/180 ;
teta2 = 60*pi/180 ;
td1 = 0 ;
td2 = 0 ;

teta_dot = [ td1; td2 ] ;
teta = [ teta1; teta2 ] ;
Y = [teta_dot; teta] ;

x1(1) = l1*sin(Y(3)) ;
y1(1) = -l1*cos(Y(3)) ;

x2(1) = x1(1) + l2*sin(Y(4)) ;
y2(1) = y1(1) - l2*cos(Y(4)) ;

for i = 1:t_final/dt
    
    Y_dot = @(Y) [ Minv(Y)*F(Y) ; Y(1); Y(2) ] ;
    
%   Runge-Kutta Part
    f1 = Y_dot(Y) ;
    f2 = Y_dot(Y + dt*f1*0.5) ;
    f3 = Y_dot(Y + dt*f2*0.5) ;
    f4 = Y_dot(Y + dt*f3) ;
    
    Y = Y + ( f1+2*f2+2*f3+f4)*dt/6 ;
    
    x1(i+1) = l1*sin(Y(3)) ;
    y1(i+1) = -l1*cos(Y(3)) ;
    
    x2(i+1) = x1(i+1) + l2*sin(Y(4)) ;
    y2(i+1) = y1(i+1) - l2*cos(Y(4)) ;
end
    
for i = 1:60:t_final/dt+1
    plot([0 x1(i)],[0 y1(i)],'k','Linewidth',2)
    hold on
    plot([x1(i) x2(i)], [y1(i) y2(i)],'k','Linewidth',2)
    hold on
    scatter(x1(i),y1(i),'ro','filled')
    hold on
    scatter(x2(i),y2(i),'bo','filled')
    hold on
    plot(x1(1:i),y1(1:i),'r--')
    hold on
    plot(x2(1:i),y2(1:i),'r--')
    axis equal
    xlim([min([0,x1,x2]) max([0,x1,x2])])
    ylim([min([0,y1,y2]) max([0,y1,y2])])
    title('Time = ', i*dt)
    
    hold off 
    drawnow
end

function f = F(Y)
global m1;
global m2;
global l1;
global l2;
global g;
f = [ -m2*l2*Y(2)*sin(Y(3)-Y(4)) - (m1+m2)*g*sin(Y(3));
    l1*(Y(1)^2)*sin(Y(3)-Y(4)) - g*sin(Y(4)) ] ;
end

function M = Minv(Y) 
global m1;
global m2;
global l1;
global l2;
global g;
M = inv([ (m1+m2)*l1 m2*l2*cos(Y(3)-Y(4));
    l1*cos(Y(3)-Y(4)) l2 ]) ;
end
    
    