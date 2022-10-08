% ME304 - EXPERIMENTAL ENGINEERING
% Experiment 7 - Centrifugal Pump Characterization

clear all
close all
clc

% orifice differential pressure (kPa):
dpo = [28.9464 23.4328 18.53948 14.19752 10.40692 7.2366 4.61764 2.61896 1.17164 0.27568 0];

% pump differential pressure (kPa):
dpp = [48.38652 57.07128 64.1018 69.47808 73.61368 77.74928 81.05776 82.712 83.95268 85.60692 98.01372];

% rotational speed (rad/s):
omega = [46.92 47.16 47.52 47.76 48 48.36 48.72 48.96 49.2 49.44 49.8];

% inlet temperature (Celcius):
T_in = [25.52 25.74 25.96 26.4 26.62 26.84 26.84 27.06 27.06 27.28 27.94];

% input power (Watt):
P = [266 254 243 230 220 207 193 182 172 166 153];

% density of water (kg/m3):
rho = [997.5928 997.5686 997.5444 997.496 997.4718 997.4476 997.4476 997.4234 997.4234 997.3992 997.3266];

% volumetric flow rate (m3/s):
Q = (0.610)*(pi)*(0.017^2)*(1./sqrt(8*rho)).*sqrt(dpo);

% pump head (m):
H = 0.048 + ((dpp*1000)./(rho*9.81));

% turbine pressure (kPa):
p_t = ((9.81).*(rho).*(H))/1000;

% mechanical power supplied to the fluid (Watt):
P_mech = (p_t*1000).*Q;

% efficiency of the pump:
E = P_mech./P;






