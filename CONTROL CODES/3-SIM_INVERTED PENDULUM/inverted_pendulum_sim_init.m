% -------------------------------------------------------------------------
% 
% 
% 09.11.2021
% Yiğit Günsür Elmacıoğlu
% 
% Inverted Pendulum Control and Simulation
% PID control may need tuning depending on the mass and length inputs
% 
% 
% -------------------------------------------------------------------------
clc
clear
close all

m = 5 ;  %mass (kg)
l = 0.5 ;  %length (m)
teta_init = 85*pi/180 ;  %initial angle (degree)
t_f = 5 ;                %total simulation time

% Simulink code runs
simulink = sim('inverted_pendulum_1.slx');

% Outputs of Simulink are taken as matrices
teta = simulink.teta.Data(:,1) ;
time = simulink.teta.Time ;
torque = simulink.torque.Data(:,1) ;
t_torque = simulink.torque.Time ;

% For limit positions of the pendulum
x_min = l*cosd(min(teta)) ;
y_min = l*sind(min(teta)) ;
x_max = l*cosd(max(teta)) ;
y_max = l*sind(max(teta)) ;

for i = 1:5:length(teta)
    x(i) = l*cosd(teta(i)) ;
    y(i) = l*sind(teta(i)) ;
        
    subplot(2,2,[1,3])
    scatter(x(i),y(i),100,'ro','filled')
    hold on
    subplot(2,2,[1,3])
    plot([-1 1],[0 0],'k')
    hold on
    subplot(2,2,[1,3])
    plot([0 x(i)],[0 y(i)],'b')
    hold on
    subplot(2,2,[1,3])
    plot([0 x_min],[0 y_min],'--')
    hold on
    if i >= find(teta==max(teta))
        subplot(2,2,[1,3])
        plot([0 x_max],[0 y_max],'--')
    end
    grid on
    title('Time = ', time(i))
    axis equal
    xlim([-(x_max+0.3) x_max+0.3])
    ylim([-0.1 l+0.1])
    hold off
    subplot(2,2,2)
    plot(t_torque(1:i),torque(1:i))
    xlim([0 max(t_torque)])
    ylim([min(torque) max(torque)])
    title('Torque')
    grid on
    subplot(2,2,4)
    plot(time(1:i),teta(1:i))
    xlim([0 max(time)])
    ylim([min(teta) max(teta)])
    grid on
    title('Teta')
    
    drawnow
end

