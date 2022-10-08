% -------------------------------------------------------------------------
% 
% Yiğit Günsür ELMACIOĞLU  -  2017405120
% 
% -------------------------------------------------------------------------

clc
clear 
close all

Q_lpm = [ 20 30 40 50] ;
H = [21 18.4 14 7.6];
Lpm_to_m3ps = 1/60000 ;    % conversion factor from lpm to kg/s
psi2Pa = 6894.75729;       % conversion factor from psi to Pascal
g = 9.81 ;
rho = 1000 ;
Q = Lpm_to_m3ps*Q_lpm ;
P_required = 5.8 * psi2Pa 
H_required = P_required / (rho*g) 
Q_operation = 57 * Lpm_to_m3ps ;
Q_operation_lpm = 57 ;

H_available = fit(Q_lpm',H','-a*x^2+b');
H_0 = H_available.b;
a = H_available.a;
func = sprintf('H_{available} = %.2f - %.2f Q^2', H_0, a) ;

H_available_lpm = fit(Q_lpm',H','-a*x^2+b');
H_0_lpm = H_available_lpm.b;
a_lpm = H_available_lpm.a;
func_lpm = sprintf('H_{available/Lpm} = %.2f - %.2f Q^2', H_0_lpm, a_lpm) ;

line = @(x) H_0 - a*x.^2 ;
line_lpm = @(x) H_0_lpm - a_lpm*x.^2 ;

if line(Q_operation) >= H_required
    disp('The chosen pump is appropriate for the job using m^3/s curve fit !')
else 
    disp('The chosen pump is insufficient using m^3/s curve fit!')
end
if line_lpm(Q_operation_lpm) >= H_required
    disp('The chosen pump is appropriate for the job using LPM curve fit !')
else 
    disp('The chosen pump is insufficient using LPM curve fit!')
end

plot(H_available)
hold on
plot(Q_lpm,H,'ko')
grid on
xlabel('\bf Q, Flowrate (lpm)')
ylabel('\bf H_{available} (m)')
legend(func,'Table Data')
