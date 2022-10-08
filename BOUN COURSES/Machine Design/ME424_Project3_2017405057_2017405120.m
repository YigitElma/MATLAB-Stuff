%--------------------------------------------------------------------------
%
%         ME424 - DESIGN PROJECT 3
%
%  EGE ORAY                   -  2017405057
%  YİĞİT GÜNSÜR ELMACIOĞLU    -  2017405120
%
%                24/01/2022
%
% -------------------------------------------------------------------------

clc
clear
close all

a = 7 ;
b = 0 ;

%  * * * * * * * * * * * GIVEN VALUES * * * * * * * * * * * * * * * * * * * 
%
n      = 200 + 3*b ;     % RPM
w      = 2*pi*n/60 ;     % Hz
F_min  = 100 + 12*a ;    % N
F_max  = 800 + 12*a ;    % N
d_L    = (20 + b/2)/1000 ;         % m
k      = ( F_max - F_min )/d_L ;   % N/m
d_max  = F_max/k ;       % m
M      = 0.3 + 0.05*a ;  % kg
G      = 79e9 ;          % Pa
E      = 200e9 ;         % Pa
rho    = 7700;           % kg / m^3
d_standard = [ 0.1 0.12 0.16 0.2:0.05:0.7 0.8:0.1:1.2 ...
               1.4:0.2:2.2 2.5 2.8 3:0.5:7 8:14 ] ;  % mm
N_d        = length(d_standard) ;
%
%  * * * * * * * * * * * GIVEN VALUES END * * * * * * * * * * * * * * * * * 

%  Following 2 vector is for storing the allowable wire diameter values for
%  the chosen material. The order of the elements is A227 - A228 - A229 - A232 - A401 
i_d_min    = [ 10 6 10 10 15 ] ;          % index of min available d for materials
i_d_max    = [ N_d 33 N_d N_d-2 N_d-3 ] ; % index of max available d for materials
d_C        = 0.05 ;                       % increment for C iterations
d_clearence   = 0.1*d_L  ;      % m
cost_min = inf ;                % For comparison initial cost is chosen as infinity

SF_fatigue_min = 1.2 ;   % Required minimum safety factor for fatigue failure
SF_surge_min = 5 ;       % Required minimum safety factor for spring surge
S_s_factor = 0.37 ;      % This is taken from the table 12.15 shot-peened and 10^7 cycle


for material = 1:5     % The order of the elements is A227 - A228 - A229 - A232 - A401
    for C = 4:d_C:12
        K_s = 1 + 0.5/C ;
        K_w = (4*C-1)/(4*C-4) + 0.615/C ;
        for i = i_d_min(material):i_d_max(material)
            
            d = d_standard(i) / 1000 ;           % wire diameter (m)
            D = C*d ;                            % coil diameter (m)
            N_a = G*d / (8*k*C^3);               % number of the active coils
            N_t = N_a + 2 ;                      % number of the total coils
            L_s = ( N_t + 1 )*d ;                % solid length (m)
            L_f = L_s + d_clearence + d_max ;    % free length (m)

            m = pi^2 * d^2 * D * N_a * rho / 4;  % mass of the active coils (kg)
         
            S_u = find_Su(material,d*1000) ;     % Ultimate tensile strength (MPa)
            S_us = 0.8*S_u ;                     % Ultimate shear strength (MPa)
            S_sy = 0.45*S_u ;                    % Shear yield strength (MPa)
            S_y = 0.75*S_u ;                     % Yield strength (MPa)

            % --------------- SPRING SURGE CHECK ------------------------ %
            
            % Right boundary condition is solved analytically in the first
            % part of the project. However, the exact solution for the
            % natural frequency requires numerical solution since this is
            % an implicit equation.
            right_boundary_condition = @(f_n) f_n * tan(f_n * sqrt(m/k)) - sqrt(m*k) / M;
            % We could use Newton-Raphson or other algorithm to find the
            % natural frequency but instead we used fzero() function of
            % MATLAB. It gives the root in a reasonable time.
            % Since the implicit equation is an even function on w_n,
            % negative and positive results can be used taking the absolute
            % value. fzero() function always finds the root which is
            % closest to the origin.
            w_n = abs(fzero(right_boundary_condition, 0));

            if w_n < SF_surge_min * w
                continue
            end
            SF_surge = w_n / w ;
            % --------------- SPRING SURGE CHECK END -------------------- %
            %
            %
            % -------------- SPRING BUCKLING CHECK ---------------------- %
            c1        = E/(2*E-2*G) ;
            c2        = 2*pi^2*(E-G)/(2*G+E) ;
            alfa      = 0.5 ;  
            lamda_eff = alfa*L_f/D ;
            d_cr      = L_f*c1*( 1-sqrt(1-c2/lamda_eff^2) ) ;  % m
            delta_S   = L_f - L_s ;     % m

            if 1-c2/lamda_eff^2 < 0
                SF_buckle = inf ;
            else
                SF_buckle = d_cr/d_max ;
            end

            if SF_buckle < 1.5
                continue
            end
            % -------------- SPRING BUCKLING CHECK END ------------------ %
            %
            %
            % ------------------ STATIC FAILURE CHECK ------------------- %
            F_s = delta_S*k;                             % Force required to compress the spring to solid state (L_s) 
            tau_static_max = K_s * 8*F_s*C/pi/d^2 / 1e6; % Corresponding shear stress due to force F_s 
            
            if tau_static_max > S_sy   % This means spring yields
                continue 
            end
            SF_yield = S_sy/tau_static_max ;  % Safety factor agains yielding 
            % -------------- STATIC FAILURECHECK END -------------------- %
            %
            %
            % ----------------- FATIGUE ANALYSIS ------------------------ %
            tau_fatigue_max = K_w * 8*F_max*C/pi/d^2 /1e6;    % Maximum shear stress (MPa)
            tau_fatigue_min = K_w * 8*F_min*C/pi/d^2 /1e6;    % Minimum shear stress (MPa)

            tau_m = (tau_fatigue_min+tau_fatigue_max)/2 ;     % Mean shear stress (MPa)
            tau_a = (tau_fatigue_max-tau_fatigue_min)/2 ;     % Alternating shear stress (MPa)

            if d < 8
                C_size = 1 ;
            else
                C_size = 1.189*(d*1000)^(-0.097) ;
            end

            S_s = S_u*S_s_factor ;
            S_n = S_us*S_s/( 2*S_us - S_s ) ;             % Fatigue strenth for desired life (MPa)
            SF_f = 1/ ( (tau_a/S_n) + (tau_m/S_us) ) ;    % Safety factor against fatigue failure 
       
            if SF_f < SF_fatigue_min    % This means spring does not satifies the requirement 
                continue
            end
            % ----------------- FATIGUE ANALYSIS END -------------------- %


            cost = find_cost(material,d*100,D*100,N_t) ;  % Cost for this iteration which satisfies the requirements 

            if cost < cost_min              % To find the minimum cost 
                cost_min = cost ;
                d_opt = d ;
                C_opt = C ;
                D_opt = D ;
                N_a_opt = N_a ;
                N_t_opt = N_t ;
                L_s_opt = L_s ;
                L_f_opt = L_f ;
                F_s_opt = F_s ;
                K_s_opt = K_s ;
                K_w_opt = K_w ;
                slender_ratio_opt =L_f/D ;
                deflect_ratio_opt = d_max/L_f ;
                SF_buckling_opt = SF_buckle ;
                mass_Na_opt = m ;
                f_n_opt = w_n/2/pi ;
                SF_surge_opt = SF_surge ;
                S_u_opt = S_u ;
                S_us_opt = S_us ;
                S_y_opt = S_y ;
                S_ys_opt = S_sy ;
                tau_s_opt = tau_static_max ; 
                SF_yield_opt = SF_yield ;
                tau_a_opt = tau_a ;
                tau_m_opt = tau_m ;
                S_s_opt = S_s ;
                SF_fatigue_opt = SF_f ;
            end
        end
    end
end

% -----------------------------RESULTS----------------------------------- %
fprintf('The results are in the units of : \n \t pressure - MPa \n \t force  - N \n \t length - m \n \t frequency - Hz\n')
F_min              % minimum force (N)
F_max              % maximum force (N)
d_opt              % optimum wire diameter (m)
C_opt              % optimum C 
D_opt              % optimum coil dimater (m)
d_clearence        % clearence (m)
k                  % stifness of the optimum spring (N/m)
N_a_opt            % number of active coils of the optimum spring
N_t_opt            % number of the total coils
L_s_opt            % solid length og the optimum spring (m)
L_f_opt            % free length of the optimum spring  (m)
F_s_opt            % force reqquired for contraction to solid length (N)
K_s_opt            % static stress factor
K_w_opt            % fatigue stress factor
slender_ratio_opt  % slenderness ratio for the optimum spring
deflect_ratio_opt  % deflection ratio for the optimum spring
SF_buckling_opt    % safety factor against buckling
mass_Na_opt        % mass of the active coils of the optimum spring
f_n_opt            % natural frequency of the optimum gear (Hz)
SF_surge_opt       % safety factor against spring surge
S_u_opt            % ultimate tensile stength of the chosen spring (MPa)
S_us_opt           % ultimate shear strength of the chosen spring (MPa)
S_y_opt            % yield strength of the chosen spring (MPa)
S_ys_opt           % shear yield strength (MPa)
tau_s_opt          % static shear stress (MPa)
SF_yield_opt       % safety factor against yielding
tau_a_opt          % alternating shear stress of the optimum design (MPa)
tau_m_opt          % mean shear stress of the optimum design (MPa)
S_s_opt            % emprical torsional stress (MPa)
SF_fatigue_opt     % safect factor against fatigue failure
% ----------------------------------------------------------------------- %
% 
% 
% ----------------------------------------------------------------------- %
function S_u = find_Su(n,d)           % This function gives ultimate tensile stress for iterated material and diameter
    if n == 1                         % Formulas are taken from the formula sheet
        S_u = 1753.3*d^(-0.1822) ;
    elseif n== 2
        S_u = 2153.5*d^(-0.1625) ;
    elseif n == 3
        S_u = 1831.2*d^(-0.1833) ;
    elseif n == 4
        S_u = 1909.9*d^(-0.1453) ;
    else
        S_u = 2059.2*d^(-0.0934) ;
    end
end
% ----------------------------------------------------------------------- %
function cost = find_cost(n,d,D,N_t)  % This function finds the cost for iterated parameters
    if n == 1                         % Formulas are taken from the given table
        cost = (100*d^2 * N_t * D) ; 
    elseif n == 2
        cost = (200*d^2 * N_t * D) ;
    elseif n == 3
        cost = (130*d^2 * N_t * D) ;
    elseif n == 4
        cost = (250*d^2 * N_t * D) ;
    else
        cost = (400*d^2 * N_t * D) ;
    end
end
% ----------------------------------------------------------------------- %