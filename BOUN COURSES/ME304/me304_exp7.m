clear all
close all
clc

load exp_data.mat

temp = 20:10:70;
rho_data = [ 998.2 997.1 995.1	994.1	992.3	988.1 ];

step = exp_data(:,1);
p_o = exp_data(:,2);
p_p = exp_data(:,3);
w = exp_data(:,4);
t = exp_data(:,5);
p_e = exp_data(:,6);

f_p_o = 6.892;  %kPa
f_p_p = 41.356; %kPa
f_w = 12;       %Hz
f_T = 22;       %degree
f_p_e =100;     %Watt

Po = f_p_o*p_o;
Pp = f_p_p*p_p;
W = f_w*w;
T = f_T*t;
Power = f_p_e*p_e;

%water density
rho = interp1(temp,rho_data,T);
%volumetric flow rate
Q = (0.610)*(pi)*(0.017^2)*(1./sqrt(8*rho)).*sqrt(Po*1000);
% pump head (m):
H = 0.048 + ((Pp*1000)./(rho*9.81));
%turbine pressure
Pt = rho.*H*9.81;
%mechanical power
Pmech = Pt.*Q;
%efficiency
E = Pmech./Power*100;

p1 =  -4.108e+07;
p2 =   5.229e+04;
p3 =      0.04098;
f_e = @(x) p1*x.^2 + p2*x + p3;

p4 =  -1.883e+07;
p5 =   1.712e+05;
p6 =       79.84;
f_power = @(x) p4*x.^2 + p5*x + p6;

p7 =  -1.073e+10;
p8 =  -2.234e+07;
p9 =   6.406e+04;
f_pt = @(x) p7*x.^2 + p8*x + p9;

x_axis = 0:0.000000001:9*10^(-4);

c_pt = f_pt(x_axis);
c_power = f_power(x_axis);
c_e = f_e(x_axis);


ylabels{1}='Turbine Pressure (Pa)';
ylabels{2}='Input Power (Watt)';
ylabels{3}='Efficiency (%)';
plotyyy(x_axis,c_pt,x_axis,c_power,x_axis,c_e,ylabels);
title('Pump Performance Characteristics for constant shaft speed')
grid on
% title('Performance')
% ylabels{1}='Turbine Pressure (kPa)';
% ylabels{2}='Input Power (Watt)';
% ylabels{3}='Efficiency (%)';
% plotyyy(Q,Pt,Q,Power,Q,E,ylabels);








