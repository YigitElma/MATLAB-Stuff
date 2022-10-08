% -------------------------------------------------------------------------
%
% 09.11.2021
% Yiğit Günsür Elmacıoğlu
% 
% Sliding Double Pendulum Simulation
% 
% -------------------------------------------------------------------------
clc
clear 
close all

global m1;
global m2;
global m3;
global l1;
global l2;
global g;

m1 = 50 ;     %mass of the car
m2 = 1 ;   %mass of the first mass
m3 = 1 ;   %mass of the second mass
l1 = 1 ;   %first cord length
l2 = 1 ;   %second cord length
g = 9.81 ;   %gravity
l_cord = l1+l2 ;
dt = 0.0001 ;   %time step between Runge-Kutta 
t_final = 10 ;

% initial conditions
teta1 = 170*pi/180 ;
teta2 = 190*pi/180 ;
xd = 0 ;     %x_dot
td1 = 0 ;    %teta1_dot
td2 = 0 ;    %teta2_dot
x1(1) = 0 ;
y1 = 0 ;

gen_vel = [ xd; td1; td2 ] ;           %generalized velocities
gen_coord = [ x1(1); teta1; teta2 ] ;  %generalized coordinates
S = [gen_vel; gen_coord] ;             %state vector

x2(1) = x1(1) + l1*sin(S(5)) ;
y2(1) = l1*cos(S(5)) ;

x3(1) = x2(1) + l2*sin(S(6)) ;
y3(1) = y2(1) + l2*cos(S(6)) ;

S_dot = @(S) [ Minv(S)*F(S) ; S(1); S(2); S(3) ] ;

for i = 1:t_final/dt   
%   Runge-Kutta Part
    f1 = S_dot(S) ;
    f2 = S_dot(S + dt*f1*0.5) ;
    f3 = S_dot(S + dt*f2*0.5) ;
    f4 = S_dot(S + dt*f3) ;
    
    S = S + ( f1+2*f2+2*f3+f4)*dt/6 ;
    
    x1(i+1) = x1(i) + S(1)*dt ;
    
    x2(i+1) = x1(i+1) + l1*sin(S(5)) ;
    y2(i+1) = y1 + l1*cos(S(5)) ;
    
    x3(i+1) = x2(i+1) + l2*sin(S(6)) ;
    y3(i+1) = y2(i+1) + l2*cos(S(6)) ;
end
    
for i = 1:600:t_final/dt+1
    plot([min([0,x1,x2,x3])-l_cord/2 max([0,x1,x2,x3])+l_cord/2],[0 0],'k','Linewidth',4)
    hold on
    plot([x1(i) x2(i)],[y1 y2(i)],'k','Linewidth',2)
    hold on
    plot([x2(i) x3(i)], [y2(i) y3(i)],'k','Linewidth',2)
    hold on
    scatter(x1(i),y1,300*m1,'gs','filled','MarkerEdgeColor','k')
    hold on
    scatter(x2(i),y2(i),m2*100,'ro','filled','MarkerEdgeColor','k')
    hold on
    scatter(x3(i),y3(i),m3*100,'bo','filled','MarkerEdgeColor','k')
    hold on
    plot(x2(1:i),y2(1:i),'r--')
    hold on
    plot(x3(1:i),y3(1:i),'r--')
    axis equal
    xlim([min([0,x1,x2,x3])-l_cord/2 max([0,x1,x2,x3])+l_cord/2])
    ylim([min([0,y1,y2,y3]) max([0,y1,y2,y3])])
    title('Time = ', i*dt)
    
    hold off 
    drawnow
end

function f = F(S)
global m1;
global m2;
global m3;
global l1;
global l2;
global g;
f = [ (m2+m3)*l1*(S(2)^2)*sin(S(5)) + m3*l2*(S(6)^2)*sin(S(6));
      (m2+m3)*l1*g*sin(S(5)) - m3*l1*l2*(S(6)^2)*sin(S(5)-S(6));
      (m3*g*l2*sin(S(6)) + m3*l1*l2*(S(5)^2)*sin(S(5)-S(6)) ) ] ;
end

function M = Minv(S) %inverse mass matrix
global m1;
global m2;
global m3;
global l1;
global l2;
M = inv([ (m1+m2+m3)  (m2*l1*cos(S(5))+m3*l1*cos(S(5)))  (m3*l2*cos(S(6)));
          ((m2+m3)*l1*cos(S(5)))  ((m2+m3)*(l1^2))  (m3*l1*l2*cos(S(5)-S(6))); 
          (m3*l2*cos(S(6)))  (m3*l1*l2*cos(S(5)-S(6)))  (m3*(l2^2))       ]) ;
end
    
    