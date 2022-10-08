
% -------------------------------------------------------------------------
%
%
%
% Ege ORAY                   2017405057
% Yiğit Günsür ELMACIOĞLU    2017405120
%
%                       ME424 DESING PROJECT - 1
%                          27 NOVERMBER 2021
% 
% 
%
% -------------------------------------------------------------------------


%   THIS CODE RUNS IN 25 SECONDS ON MY COMPUTER, IF YOU HAVE TROUBLE WITH
%           FAST RESULTS PLEASE CHANGE LOOP LIMITS OF b AS;
%                  for b = 14*module:1:14*module 
clc
clear
close all
format longG

tic

a = 7 ;
b = 0 ;

% Given quantities
w_s = (3000 + 100*a)*2*pi/60 ;  % rad/sec
w_r = 0 ;                       % rad/sec
d_max = 450 + 10*a - 2*5 ;      % mm
w_max = (180 + 10*b)*2*pi/60 ;  % rad/sec arm
w_min = (90 + 5*b)*2*pi/60 ;    % rad/sec
phi = 20*pi/180 ;               % pressure angle  (rad)
SF_b = 1.5 ;    % required safety factor for bending
SF_H = 1.2 ;    % required safety factor for surface fatigue
K_o = 1.75 ;
k_r = 0.814 ;
C_L = 1 ;
k_t = 1 ;
k_ms = 1.4 ;
C_R = 1 ;
C_p = 191;

% For Geometry factor J, we first used a 2 dimensional interpolation but it
% turned to take a lot of time. So, we replaced it with a function of N1
% and N2.
% 
% N1 = [10 15 17 20 24 30 35 40 45 50 60 80 125 275 600] ;
% N2 = [10 17 25 35 50 85 170 1000] ;
% 
% J_values = [0.2 0.25 0.29 0.30 0.32 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.44 0.45;
%             0.2 0.25 0.29 0.31 0.33 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.45 0.46;
%             0.2 0.25 0.29 0.32 0.34 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.45 0.46 0.47;
%             0.2 0.25 0.29 0.33 0.35 0.37 0.39 0.40 0.41 0.42 0.43 0.44 0.46 0.47 0.48;
%             0.2 0.25 0.29 0.33 0.35 0.38 0.40 0.41 0.42 0.43 0.44 0.46 0.47 0.49 0.50;
%             0.2 0.25 0.29 0.34 0.36 0.39 0.41 0.42 0.43 0.44 0.45 0.47 0.48 0.50 0.52;
%             0.2 0.25 0.29 0.34 0.37 0.40 0.42 0.43 0.44 0.45 0.46 0.48 0.50 0.52 0.53;
%             0.2 0.25 0.29 0.35 0.37 0.41 0.42 0.44 0.45 0.46 0.47 0.49 0.51 0.53 0.55 ] ;
%         
% J = interp2(N1, N2, J_values, N_gear(Gear_number), N_gear(Gear_mate) ) ;

module_available = [ 0.2:0.1:1 1.25:0.25:4 4.5:0.5:6 7:1:10 12:2:20 25 32 40 50] ;   % mm 
m = length(module_available) ;
Power_max = 0 ;                 % Watt
loop_count = 1 ;
loop_actual = 1 ;

% For optimum module, we try every possible value.
for i = 1:m
    module = module_available(i) ;   % mm
%   Higher HB factor increases strength with no drawback. To take less computation time, we only use HB=200,400 
    for HB = 200:200:400
        S_u = 3.45*HB ;       % Ultimate strength in MPa
        S_n_prime = S_u/2 ;   % MPa
        S_fe = 2.76*HB - 69 ; % MPa
        C_s = -1.004e-6*HB^2 - 5.446e-5*HB + 0.8114;  % Surface factor for fatigue calculations
%       Gear width must lie between 9module and 14 module
        for b = 9*module:1:14*module 
            % Mounting factor depends on face width
            if  b <= 50
                K_m = 1.6 ;
            elseif b < 405
                K_m = 2.931e-6*b^2 + 0.0003506*b + 1.577 ;
            else
                K_m = 2.2 ;
            end
%           To find integer number of teeth values, we use N as design variable
%           Number of gear teeth of sun (upper limit is for 25m/s velocity
%           requirement
            for N_s = 12:1:50000/(w_s*module)
                d_s = module*N_s ;              % diameter of sun gear in mm
%               Addendum of Planet 1 must not pass housing limits
                for N_p1 = 12:1:( d_max-d_s )/(2*module) - 1
                    d_p1 = module*N_p1 ;        % diameter of planet 1 in mm
%                   Dedendum of ring must obey to housing limit, so N of planet2 has upper limit as follows   
                    for N_p2 = 12:1:( d_max/module -2.5 - N_s - N_p1 )
                        d_p2 = module*N_p2 ;    % diameter of planet 2 in mm
                        w_arm = N_p2*N_s*w_s / ( (N_s+N_p1)*(N_p1+N_p2) ) ;   % angular velocity of arm inrad/sec
                        d_r = module*( N_p1+N_p2+N_s ) ;                      % diameter of the ring in mm
                        N_r = d_r/module ;                                    % Number of gear teeth of ring
                        w_p = w_s*N_s / ( N_p1+N_p2 ) ;                       % angular velocity of planets in rad/sec
                        
                        loop_actual = loop_actual + 1 ;
                        % angular velocity limit check for arm
                        if w_arm < w_min || w_arm > w_max
                            break
                        end
                        % interference check for internal gears (ring and
                        % planet 2)
                        if N_r - N_p2 < 10
                            break
                        end
                        
                        c = (d_p1 + d_s)/2 ;    % central distance between sun and planet 1 in mm
                        c_1 = (d_r - d_p2)/2 ;  % central distance between ring and planet 2 in mm
                        % interference check for external gears (sun and
                        % planet1)
                        if sqrt((d_p1/2*cos(phi))^2 + (c*sin(phi))^2) < d_p1/2 + module ...
                                || sqrt((d_s/2*cos(phi))^2 + (c*sin(phi))^2) < d_s/2 + module
                            break
                        end
                        % 25m/s velocity limit for planet 2 at inner edge
                        if (d_s+d_p1)/2*w_arm + w_p*d_p2/2 > 25000
                            break
                        end
                        % contact ratio for external gears (sun and
                        % planet1)
                        CR_E = (-c*sin(phi) + sqrt((d_s/2+module)^2 - ((d_s/2)*cos(phi))^2 )...
                            + sqrt((d_p1/2+module)^2 - ((d_p1/2)*cos(phi))^2 ) ) / (pi*module*cos(phi)) ;
                        % contact ratio for internal gears (ring and
                        % planet2)
                        CR_I = (c_1*sin(phi) + sqrt((d_p2/2+module)^2 - ((d_p2/2)*cos(phi))^2 )...
                            - sqrt((d_r/2-module)^2 - ((d_r/2)*cos(phi))^2 ) ) / (pi*module*cos(phi)) ;
                        % contact ratio must be bigger than 1.4
                        if CR_E < 1.4 || CR_I < 1.4
                            break
                        end
                        N_gear = [N_s N_p1 N_p2 N_r];                           
                        % for each gear, we make stress and force
                        % calculations which has minimum allowable safety
                        % factor
                        for Gear_number = 1:1:4
                            if Gear_number == 1      % sun gear
                                v = w_s*d_s/2/1000 ; % velocity for K_v 
                                life(Gear_number) = 8 + log10(3*(w_s/w_arm - 1)) ;
                                Gear_mate = 2 ;      % sun is mated with planet 1
                                % pinion and gear diameters
                                if d_s < d_p1
                                    d_p = d_s ;
                                    d_g = d_p1 ;
                                else
                                    d_p = d_p1 ;
                                    d_g = d_s ;
                                end
                            elseif Gear_number == 2  % planet 1
                                v = w_s*d_s/2/1000 ;
                                life(Gear_number) = 8 + log10((w_p+w_arm)/w_arm) ;
                                Gear_mate = 1 ;
                                if d_s < d_p1
                                    d_p = d_s ;
                                    d_g = d_p1 ;
                                else
                                    d_p = d_p1 ;
                                    d_g = d_s ;
                                end
                            elseif Gear_number == 3  % planet 2
                                v = (w_arm*(d_s+d_p1)/2 + w_p*d_p2/2)/2000 ;
                                life(Gear_number) = 8 + log10((w_p+w_arm)/w_arm) ;
                                Gear_mate = 4 ;
                                if d_r < d_p2
                                    d_p = d_r ;
                                    d_g = d_p2 ;
                                else
                                    d_p = d_p2 ;
                                    d_g = d_r ;
                                end
                            else                     % ring gear
                                v = (w_arm*(d_s+d_p1)/2 + w_p*d_p2/2)/2000 ;
                                life(Gear_number) = 8 + log10(3) ;
                                Gear_mate = 3 ;
                                if d_r < d_p2
                                    d_p = d_r ;
                                    d_g = d_p2 ;
                                else
                                    d_p = d_p2 ;
                                    d_g = d_r ;
                                end
                            end
                               
                            if module <= 5
                                C_g = 1;
                            else
                                C_g = 0.85 ;
                            end
                            
                            K_v(Gear_number) = 1 + v/6 ;
                            % we used curve fitting app and custom function
                            % option to have a good fit for J values
                            J(Gear_number) = 3.026*N_gear(Gear_number)^0.03603*N_gear(Gear_mate)^0.008841...
                                - 0.0007267*N_gear(Gear_number) - 0.0001756*N_gear(Gear_mate) - 3.142 ; 
                            
                            % we used curve fitting app and power function
                            % for life coefficient fitting
                            C_Li(Gear_number) = 8.66*(life(Gear_number)^-1.307)+0.3111 ;
                            R = d_g/d_p ;
                            if Gear_number < 3   % for external gear pair (sun and planet1)
                                I(Gear_number) = sin(phi)*cos(phi)/2*R/(R+1) ;
                            else                 % for internal gear pair (ring and planet2)
                                I(Gear_number) = sin(phi)*cos(phi)/2*R/(R-1) ;
                            end
                            
                            S_n(Gear_number) = S_n_prime*C_L*C_g*C_s*k_r*k_t*k_ms ;  % MPa
                            sigma_b_max(Gear_number) = S_n(Gear_number)/SF_b ;       % MPa
                            
                            SH(Gear_number) = S_fe*C_Li(Gear_number)*C_R ;          % MPa
                            sigma_H_max(Gear_number) = SH(Gear_number)/sqrt(SF_H) ; % MPa
                            
                            % Maximum allowable tangential force for bending case in Newton
                            Ft_b = sigma_b_max(Gear_number)/( K_v(Gear_number)*K_o*K_m/(module*b*J(Gear_number)) ) ;    
                            % Maximum allowable tangential force for surface fatigue case in Newton
                            Ft_H = (sigma_H_max(Gear_number)^2)/( (C_p^2)*K_v(Gear_number)*K_o*K_m/( b*d_p*I(Gear_number) ) ); 
                            
                            % Limiting force is the smaller one 
                            lim_F_gear(Gear_number) = min( Ft_H, Ft_b ) ;
                            loop_count = loop_count + 1 ;  
                            d_p_cycle(Gear_number) = d_p ;
                        end
                        
                        Power_s = lim_F_gear(1)*( 3*w_s*d_s )/2000 ;
                        Power_P1 = lim_F_gear(2)*( 3*w_s*d_s )/2000 ;
                        Power_r = lim_F_gear(3)*( 3*d_p2*d_s*w_s )/(2000*d_p1) ;
                        Power_P2 = lim_F_gear(4)*( 3*d_p2*d_s*w_s )/(2000*d_p1) ;
                        Power = min( [Power_s Power_P1 Power_r Power_P2] ) ;
                        limiting_gear = find( Power == [Power_s Power_P1 Power_r Power_P2] ) ;
                        
                        if Power > Power_max
                            Power_max = Power ;
                            Power_all = [Power_s Power_P1 Power_r Power_P2] ;
                            J_all = J ;
                            K_v_all = K_v ;
                            I_all = I ;
                            SH_all = SH ;
                            S_n_all = S_n ;
                            d_p_critical = d_p_cycle ; 
                            module_opt = module ;
                            b_opt = b ;
                            N_planet1_opt = N_p1 ;
                            N_planet2_opt = N_p2 ;
                            N_sun_opt = N_s ;
                            N_ring_opt = N_r ;
                            HB_opt = HB ;
                            w_arm_opt = w_arm ;
                            limiting_gear_final = limiting_gear ;
                            d_sun_opt = d_s ;
                            d_planet1_opt = d_p1 ;
                            d_planet2_opt = d_p2 ;
                            d_ring_opt = d_r ;
                            CR_sun_planet1_opt = CR_E ;
                            CR_ring_planet2_opt = CR_I ;
                            J_critical = J(limiting_gear) ;
                            K_v_critical = K_v(limiting_gear) ;
                            life_critical = life(limiting_gear) ;
                            C_Li_critical = C_Li(limiting_gear) ;
                            I_critical = I(limiting_gear) ;
                            S_n_critical = S_n(limiting_gear) ;
                            w_p_opt = w_p ;
                            c_opt = c ;
                            d_p_opt = d_p ;
                            Ft_critical = lim_F_gear(limiting_gear) ;
                            SH_critical = SH(limiting_gear) ;
                            S_prime_critical = S_n_prime ;
                        end
                    end
                end
            end
        end
    end
end

if module_opt <= 5
    C_g_critical = 1;
else
    C_g_critical = 0.85 ;
end

if  b_opt <= 50
    K_m_critical = 1.6 ;
elseif b_opt < 405
    K_m_critical = 2.931e-6*b_opt^2 + 0.0003506*b_opt + 1.577 ;
else
    K_m_critical = 2.2 ;
end
% ----------------------------- OUTPUTS -----------------------------------

disp('Max Power is in Watts')
Power_max
HB_opt
disp('n is in RPM')
n_s = w_s*30/pi
n_a = w_arm_opt*30/pi
disp('Notice that module and diameter values are in milimeter units')
module_opt
d_sun_opt
d_planet1_opt
d_planet2_opt
d_ring_opt
N_sun_opt
N_planet1_opt
N_planet2_opt
N_ring_opt
d_housing = max( [( d_sun_opt + 2*d_planet1_opt + module_opt + 2*5  ) (d_ring_opt + 2.5*module_opt + 2*5)] )
n_p = w_p_opt*30/pi
v_s = d_sun_opt/2 * w_s /1000
Ft_sun_planet1 = 2*Power_max / (3*w_s*d_sun_opt)*1000
Ft_ring_planet2 = d_planet1_opt/d_planet2_opt * 2*Power_max /(3*w_s*d_sun_opt)*1000
F_arm = Ft_sun_planet1 + Ft_ring_planet2
b_opt
CR_sun_planet1_opt
CR_ring_planet2_opt
delta_interference_planet1 = sqrt((d_planet1_opt/2*cos(phi))^2 + (c_opt*sin(phi))^2) - (d_planet1_opt/2 + module_opt)
delta_interference_sun = sqrt((d_sun_opt/2*cos(phi))^2 + (c_opt*sin(phi))^2) - (d_sun_opt/2 + module_opt)
delta_N_ring_planet2 = N_ring_opt - N_planet2_opt
limiting_gear_final
J_critical
K_v_critical
K_m_critical
K_o_critical = 1.75

%------------------Critical Gear Safety Factors----------------------------
sigma_b_critical = Ft_critical*K_v_critical*K_o_critical*K_m_critical/module_opt/b_opt/J_critical ;
sigma_H_critical = C_p*sqrt( Ft_critical*K_v_critical*K_o_critical*K_m_critical/b_opt/d_p_opt/I_critical ) ;
SF_b_critical = S_n_critical/sigma_b_critical ;
SF_H_critical = ( SH_critical/sigma_H_critical )^2;
%--------------------------------------------------------------------------

%------------------SUN Gear Safety Factors---------------------------------
sigma_b_sun = Ft_sun_planet1*K_v_all(1)*K_o_critical*K_m_critical/module_opt/b_opt/J_all(1) ;
sigma_H_sun = C_p*sqrt( Ft_sun_planet1*K_v_all(1)*K_o_critical*K_m_critical/b_opt/d_p_critical(1)/I_all(1) ) ;
SF_b_sun = S_n_all(1)/sigma_b_sun ;
SF_H_sun = ( SH_all(1)/sigma_H_sun )^2;
%--------------------------------------------------------------------------

%------------------Planet 1 Gear Safety Factors----------------------------
sigma_b_p1 = Ft_sun_planet1*K_v_all(2)*K_o_critical*K_m_critical/module_opt/b_opt/J_all(2) ;
sigma_H_p1 = C_p*sqrt( Ft_sun_planet1*K_v_all(2)*K_o_critical*K_m_critical/b_opt/d_p_critical(2)/I_all(2) ) ;
SF_b_p1 = S_n_all(2)/sigma_b_p1 ;
SF_H_p1 = ( SH_all(2)/sigma_H_p1 )^2;
%--------------------------------------------------------------------------

%------------------RING Gear Safety Factors--------------------------------
sigma_b_ring = Ft_ring_planet2*K_v_all(4)*K_o_critical*K_m_critical/module_opt/b_opt/J_all(4) ;
sigma_H_ring = C_p*sqrt( Ft_ring_planet2*K_v_all(4)*K_o_critical*K_m_critical/b_opt/d_p_critical(4)/I_all(4) ) ;
SF_b_ring = S_n_all(4)/sigma_b_ring ;
SF_H_ring = ( SH_all(4)/sigma_H_ring )^2;
%--------------------------------------------------------------------------

sigma_b_critical
S_prime_critical
C_L_critical = 1
C_g_critical
C_s_critical = -1.004e-6*HB_opt^2 - 5.446e-5*HB_opt + 0.8114
k_r_critical = 0.814
k_ms_critical = 1.4
S_n_critical 
SF_b_critical 
C_p
I_critical
sigma_H_critical
Sf_e_critical = 2.76*HB_opt - 69 
number_contact = 10^life_critical
C_Li_critical
C_R_critical = 1 
SH_critical
SF_H_critical 

SF_array = [SF_b_sun  SF_H_sun ;
            SF_b_p1 SF_H_p1 ;
            SF_b_ring  SF_H_ring ] ;

toc

