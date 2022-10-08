clc
clear
close all

m1 = 0.5 ;     %mass of the car
m2 = 0.536 ;     %mass of the first mass
m3 = 0.176 ;     %mass of the second mass
l1 = 0.379 ;     %first cord length
l2 = 0.315 ;     %second cord length
g = 9.81 ;   %gravity
l_cord = l1+l2 ;
t_final = 15 ;
y1 = 0 ;

teta1_init = 0*pi/180 ;
teta2_init = 0.5*pi/180 ;
x_init = 0 ;
%%
sim   = sim('SIM_sliding_double_pendulum.slx') ;
x1    = sim.x1.Data(:,1) ;
teta1 = sim.teta1.Data(:,1) ;
teta2 = sim.teta2.Data(:,1) ;
time  = sim.teta1.Time ;

for i = 1:length(time)
    x2(i) = x1(i) + l1*sin(teta1(i)) ;
    y2(i) = l1*cos(teta1(i)) ;
    x3(i) = x2(i) + l2*sin(teta2(i)) ;
    y3(i) = y2(i) + l2*cos(teta2(i)) ;
end

x2 = x2' ;
y2 = y2' ;
x3 = x3' ;
y3 = y3' ;

for i = 1:60:length(time)
    subplot(1,2,1)
    plot([min([0;x1;x2;x3])-l_cord/2 max([0;x1;x2;x3])+l_cord/2],[0 0],'k','Linewidth',4)
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
    xlim([min([0;x1;x2;x3])-l_cord/2 max([0;x1;x2;x3])+l_cord/2])
    ylim([min([0;y1;y2;y3]) max([0;y1;y2;y3])])
    title('Time = ', time(i))
    hold off 
    
    subplot(1,2,2)
    plot(time(1:i),teta2(1:i)*180/pi)
    drawnow
end



