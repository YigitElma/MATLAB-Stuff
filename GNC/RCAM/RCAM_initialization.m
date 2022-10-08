% -------------------------------------------------------------------------
% 
% 04.11.2021
% Yiğit Günsür ELMACIOĞLU
% 
% https://www.youtube.com/watch?v=m5sEln5bWuM&list=PLxdnSsBqCrrEx3A6W94s...
% QGClk6Q4YCg-h&index=16&ab_channel=ChristopherLum
% 
% adresindeki uçak modeli kullanılarak oluşturuldu.
% -------------------------------------------------------------------------
clear 
clc
close all

xini = 0 ;
yini = 0 ;

enlem  = 48 + 51/60 + 24/3600 ;  %Paris koordinatları
boylam = 2 + 21/60 + 3/3600 ;
lon0 = deg2rad(boylam) ;
lat0 = deg2rad(enlem) ;
h0 = 580 ;

x0 = [85; 0; 0; 0; 0; 0; 0; 0; 0; xini; yini; h0] ;  % initial state vector
% u = [0; -0.1; 0; 0.08; 0.08] ;                         % initial control vector
u = [0; 0;0;0;0] ;         % initial control vector
  
u1min = -25*pi/180 ;   % aileron
u1max = 25*pi/180 ;

u2min = -25*pi/180 ;   % elevator
u2max = 10*pi/180 ;

u3min = -30*pi/180 ;   % rudder
u3max = 30*pi/180 ;

u4min = 0.5*pi/180 ;   % throttle 1
u4max = 10*pi/180 ;

u5min = 0.5*pi/180 ;   % throttle 2
u5max = 10*pi/180 ;

Xgeodetic0 = [ lat0; lon0; h0 ] ;

TF = 10 ;
%%
simOUT = sim('RCAM_simulink.slx') ;

timeX = simOUT.simX.Time ;
timeU = simOUT.simU.Time ;

X_1 = simOUT.simX.Data(:,1) ;
X_2 = simOUT.simX.Data(:,2) ;
X_3 = simOUT.simX.Data(:,3) ;
X_4 = simOUT.simX.Data(:,4) ;
X_5 = simOUT.simX.Data(:,5) ;
X_6 = simOUT.simX.Data(:,6) ;
X_7 = simOUT.simX.Data(:,7) ;
X_8 = simOUT.simX.Data(:,8) ;
X_9 = simOUT.simX.Data(:,9) ;
X_10 = simOUT.simX.Data(:,10) ;
X_11 = simOUT.simX.Data(:,11) ;
X_12 = simOUT.simX.Data(:,12) ;

U_1 = simOUT.simU.Data(:,1) ;
U_2 = simOUT.simU.Data(:,2) ;
U_3 = simOUT.simU.Data(:,3) ;
U_4 = simOUT.simU.Data(:,4) ;
U_5 = simOUT.simU.Data(:,5) ;

% State values
subplot(3,4,1)
plot(timeX,X_1)
grid on
title('u - x vel body frame')

subplot(3,4,5)
plot(timeX,X_2)
grid on
title('v - y vel body frame')

subplot(3,4,9)
plot(timeX,X_3)
grid on
title('w - z vel body frame')

subplot(3,4,2)
plot(timeX,X_4)
grid on
title('p - body frame')

subplot(3,4,6)
plot(timeX,X_5)
grid on
title('q - body frame')

subplot(3,4,10)
plot(timeX,X_6)
grid on
title('r - body frame')

subplot(3,4,3)
plot(timeX,X_7)
grid on
title('phi')

subplot(3,4,7)
plot(timeX,X_8)
grid on
title('teta')

subplot(3,4,11)
plot(timeX,X_9)
grid on
title('ksi')

subplot(3,4,4)
plot(timeX,X_10)
grid on
title('North Position')

subplot(3,4,8)
plot(timeX,X_11)
grid on
title('East Position')

subplot(3,4,12)
plot(timeX,X_12)
grid on
title('Altitude')

% Control inputs
figure
subplot(5,1,1)
plot(timeU,U_1)
grid on
title('u1 - aileron')

subplot(5,1,2)
plot(timeU,U_2)
grid on
title('u2 - elevator')

subplot(5,1,3)
plot(timeU,U_3)
grid on
title('u3 - rudder')

subplot(5,1,4)
plot(timeU,U_4)
grid on
title('u4 - throttle 1')

subplot(5,1,5)
plot(timeU,U_5)
grid on
title('u5 - throttle 2')








