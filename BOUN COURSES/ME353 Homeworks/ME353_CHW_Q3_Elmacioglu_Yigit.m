clc
clear
close all

eta = 0:0.5:8;
d_f = [ 0 0.1659 0.3298 0.4868 0.6298 0.7513 0.846 0.913 0.9555 0.9795 0.9915...
    0.9969 0.999 0.9997 0.9999 1 1 ];
N = length(eta);
h = 0.5;

%disturbance thickness
delta = 4.5 + (0.99-0.9795)/(0.9915-0.9795) * 0.5

int_disp = 1 - d_f;
int_phi = d_f.* (1 - d_f);

%displacement thickness
delta_disp = h/2 * ( int_disp(1)+int_disp(N) ) + h*sum(int_disp(2:N-1)) 
%momentum thickness
delta_phi = h/2 * ( int_phi(1)+int_phi(N) ) + h*sum(int_phi(2:N-1))

%this part is to calculate third derivative of f
dd_f = [ 0.3321 0.3309 0.323 0.3026 0.2668 0.2174 0.1614 0.1078 0.0642 0.034...
    0.0159 0.0066 0.0024 0.0008 0.0002 0.0001 0];
for i=2:N-1
    ddd_f(i) = ( dd_f(i+1)-dd_f(i-1) )/0.5 ;
end
ddd_max=max(abs(ddd_f))
%error caused by Trapezoidal Rule integration (only for displacement thickness)
error_trap_disp = 8*h^2*ddd_max/12
%percent errors relative to values on the book
percent_error_delta = (abs(5-delta)/5)*100
percent_error_disp = (abs(1.721-delta_disp)/1.721)*100
percent_error_phi = (abs(0.664-delta_phi)/0.664)*100

tau_0 = dd_f(1)%shear stress at the surface in terms of (mu.U/x*sqrt(Re_x))
F_drag = 2*tau_0 %drag force in terms of b*L*rho*U^2/sqrt(Re_L)
c_drag = 2*F_drag %friction drag coefficient in terms of 1/sqrt(Re_L) 
percent_error_tau0 = (abs(0.332-tau_0)/0.332)*100
percent_error_F_drag = (abs(0.664-F_drag)/0.664)*100
percent_error_c_drag = (abs(1.328-c_drag)/1.328)*100

plot(d_f, dd_f)
xlabel('Velocity (u/U)')
ylabel('Shear Stress ( mu.U/x*sqrt(Re_x) )')

figure
plot(eta,d_f)
xlabel('Heigth ( x/sqrt(Re_x) )')
ylabel('Velocity (u/U)')

figure
plot(eta,dd_f)
xlabel('Heigth ( x/sqrt(Re_x) )')
ylabel('Shear Stress ( mu.U/x*sqrt(Re_x) )')
