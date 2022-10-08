clc
clear 
close all

rho = 995.7;
mu = 0.801*10^-3;
L = 30;
D = 0.4;
x0 = 0.015;
rel_roughness = 0.15/D/1000;

m_dot = 15:0.001:30;
V = 4.*m_dot/rho/pi/D^2;
N=length(m_dot);

R = rho.*V*D/mu;

for i=1:N
    f_colebrook = @(x) 2*sqrt(x)*log10(rel_roughness/3.7 + 2.51/R(i)/sqrt(x)) + 1;
    options = optimset('TolFun',0.1);
    f_a(i) = fzero(f_colebrook, x0,options);
end
for i=1:N
    f_haaland = @(x) 1.8*sqrt(x)*log10( (rel_roughness/3.7)^1.11 + 6.9/R(i)) + 1;
    f_b(i) = fzero(f_haaland, x0);
end
    
h_a = f_a.*L/D.*V.^2/2/9.81;
h_b = f_b.*L/D.*V.^2/2/9.81;

plot(m_dot, h_a, 'r', m_dot, h_b, 'b')
grid on
xlabel('Mass Flowrate (kg/s)')
ylabel('Major Head Loss (m)')
title('Major Head Loss vs Mass Flowrate')
legend('from Colebrook','from Haaland')

error = abs(h_a-h_b)./h_a*100;
figure
plot(m_dot, error)
grid on
title('Percent Error')
xlabel('Mass Flowrate (kg/s)')
ylabel('% Error')





