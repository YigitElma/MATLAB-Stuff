% -------------------------------------------------------------------------
%
% 09.11.2021
% Yiğit Günsür Elmacıoğlu
% 
% Inverted Double Pendulum Control code initialization and visualisation
% 
% -------------------------------------------------------------------------
clc
clear 
close all

m1 = 3 ;
m2 = 1 ;
l1 = 3 ;
l2 = 1 ;
g = 9.81 ;

teta1_i = 175*pi/180 ;
teta2_i = 180*pi/180 ;

teta1d_i = 0 ;
teta2d_i = 0 ;

tf = 10 ;

sim = sim('double_pendulum_simulink.slx') ;
teta1 = sim.teta1.Data(:,1) ;
teta2 = sim.teta2.Data(:,1) ;
time = sim.teta1.Time ;

for i = 1:1:length(teta1)
    x1(i) = l1*sin(teta1(i)) ;
    y1(i) = -l1*cos(teta1(i)) ;
    x2(i) = x1(i) + l2*sin(teta2(i)) ;
    y2(i) = y1(i) - l2*cos(teta2(i)) ;
end
for i = 1:1:length(teta1)
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
    title('Time = ', time(i))
    
    hold off 
    drawnow
end
