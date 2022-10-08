close all
clear 
clc

L   = 1;
A_0 = 1e-3;
F   = 100;
E   = 209e9;
s0 = F/A_0 ;

x = 0:0.001:L ;
for i = 1:length(x)
    y(i) = 2*L*F/A_0/E*log( 2*L/(2*L-x(i)) );
    sigma(i) = F/A_0/(1-x(i)/2*L) /s0;
end

% for a
x_a = [0 1] ;
delta_a = [0 4/3*F*L/A_0/E] ;
sigma_a = [4*F/3/A_0/s0 4*F/3/A_0/s0] ;

% for b
x_b = [0 0.5 1] ;
delta_b = [0 4/7*F*L/A_0/E 48/35*F*L/A_0/E] ;
x_bs = [0 0.5 0.5 1] ;
sigma_b = [8*F/7/A_0/s0 8*F/7/A_0/s0 8*F/5/A_0/s0 8*F/5/A_0/s0] ;

% plot(x, y, 'r','LineWidth',2)
% hold on
% plot(x_a, delta_a,'b','LineWidth',2)
% hold on
% plot(x_b, delta_b,'k','LineWidth',2)
% legend('Exact','for a', 'for b')
% xlabel('position x/L')
% ylabel('displacement')
% title('for L = 1m - A_0 = 10cm^2 - P = 100N - E = 209GPa')
% grid on


plot(x, sigma, 'r','LineWidth',2)
hold on
plot(x_a, sigma_a,'b','LineWidth',2)
hold on
plot(x_bs, sigma_b,'k','LineWidth',2)
legend('Exact','for a', 'for b')
xlabel('position x/L')
ylabel('sigma(x)/sigma_0')
title('for L = 1m - A_0 = 10cm^2 - P = 100N - E = 209GPa')
grid on