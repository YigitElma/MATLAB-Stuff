% -------------------------------------------------------------------------
%
%         ME424 - DESIGN PROJECT 2
%
%  EGE ORAY                   -  2017405057
%  YİĞİT GÜNSÜR ELMACIOĞLU    -  2017405120
%
%                01/12/2022
%
% -------------------------------------------------------------------------

clc
clear
close all

a = 2 ;
b = 8 ;

D     = ( 120-4*b )/1000 ;    % m
L     = ( 240+12*a )/1000 ;
t     = ( 4+0.3*a )/1000 ;
h     = 20/1000 ;
L_b   = L + 2*h ;
D_c   = D - 2*t ;
A_ci  = pi*( D_c/2 )^2 ;
A_co  = pi*(D/2)^2 ;
A_c   = A_co - A_ci ;

E = 207e9 ;
v = 0.3 ;
C_load = 1 ;
C_rel = 0.753 ;
C_size = 1 ;
C_surf = 1 ;

P_max   = ( 12 + 0.5*b )*1e6 ;
K_i_min = 0.5 ;
K_i_max = 0.9 ;
SF_min  = 1.6 ;

d_K_i = 0.0001 ;
c_Sf_life = (1-log10(5))/3 ;

Material_class = [ 12.9 10.9 9.8 8.8 5.8 4.8 4.6 ] ;
material_number = length(Material_class) ;
S_pyu = [ 970 1100 1220 ;
          830 940 1040 ;
          650 720 900 ;
          600 660 830 ;
          380 420 520 ;
          310 340 420 ;
          225 240 400 ] ;
d_standard = [ 2:0.5:4 5:8 10:2:16 20:2:24 27:3:42 48 ] ;
p_standard = [ 0.4 0.45 0.5 0.6 0.7 0.8 1 1 1.25 1.5 1.75 2 2 2.5 2.5 3 3 3.5 3.5 4 4 4.5 5 ] ;

d_i_min = [ 1 6 1 13 6 1 6 ] ;
d_i_max = [ 20 20 13 20 16 13 20 ] ;

Price_min = inf ;

for material = 1:material_number
    S_p = S_pyu(material,1)*1e6 ;  %Pa
    S_y = S_pyu(material,2)*1e6 ;
    S_u = S_pyu(material,3)*1e6 ;

    Sn_prime = 0.5*S_u ;
    Sf_1000  = 0.9*S_u ;
    S_n = C_rel*C_size*C_surf*C_load*Sn_prime ;
        
    for i = d_i_min(material):d_i_max(material)
        d   = d_standard(i)/1000 ;     %m
        p   = p_standard(i)/1000 ;
        d_p = d - p*(3/8*sqrt(3)) ;
        d_r = d - p*(17/24*sqrt(3)) ;
        A_b = pi*(d/2)^2 ;
        A_t = pi*( (d_p+d_r)/4 )^2 ;

        Price_array = [ 11 + 2.2*d*1000   ;
                        9 + 1.8*d*1000    ;
                        6 + 1.2*d*1000    ;
                        3 + 0.6*d*1000    ;
                        1.3 + 0.26*d*1000 ;
                        1.1 + 0.22*d*1000 ;
                        1 + 0.2*d*1000  ] ;
        Price = Price_array(material) ;

        k_b = E*A_b/L_b ;
        k_c = E*A_c/L ;

        F_e = P_max*( A_ci - A_b ) ;

        for K_i = K_i_min:d_K_i:K_i_max          
            F_i = K_i*S_p*A_t ;

            delta_Fb = P_max*k_b/(k_b+k_c) * ( pi/4*(D_c^2-d^2) - k_c*v*L/E*(D_c/2/t+1) ) ;
            delta_Fc = P_max*k_b*k_c/(k_b+k_c) * ( L*v/E*(1+D_c/2/t) + pi/4/k_b*(D_c^2-d^2) ) ;

            SF_yield = (S_y*A_t - F_i)/delta_Fb ;
            SF_sep   =  F_i/delta_Fc ;

            if SF_yield < SF_min
                continue
            end

            if SF_sep < SF_min
                continue
            end

            if material < 5
                K_f = 3 ;
            else
                K_f = 2.2 ;
            end

            sigma_a_axial = delta_Fb/2/A_t ;
            sigma_a_p     = P_max/2 ;
            sigma_a       = sqrt(sigma_a_axial^2 + sigma_a_p^2 - sigma_a_axial*sigma_a_p) ;

            sigma_m_axial = (F_i + delta_Fb/2)/A_t ;
            sigma_m_p     = P_max/2 ;
            sigma_m       = max(sigma_m_axial, sigma_m_p) ;

            sigma_init = F_i/A_t ;

            Sf_life  = S_n * 10^( c_Sf_life*(log10(0.9*S_u/S_n)) ) ;

            S_a = (S_u - K_f*sigma_init) / (1 + S_u/Sf_life) ;
            S_m =  S_a + K_f*sigma_init ;

            if  S_m + S_a > S_y
                S_a = Sf_life * (S_u - S_y) / (S_u - Sf_life) ;
                S_m = S_y - sigma_a*K_f ;
            end
  
            SF_fatigue = S_a/sigma_a/K_f ;

            if SF_fatigue < SF_min
                continue
            end

            if Price < Price_min
                Price_min = Price ;
                K_i_opt = K_i ;
                d_opt = d ;
                p_opt = p ;
                material_opt = material ;
                SF_fatigue_opt = SF_fatigue ;
                SF_yield_opt = SF_yield ;
                SF_sep_opt = SF_sep ;
                SF_leak_opt = SF_sep_opt ;
                d_p_opt = d_p ;
                A_t_opt = A_t ;
                F_i_opt = F_i ;
                F_e_max = F_e ;
                Fb_opt  = F_i + delta_Fb ;
                Fc_opt  = F_i - delta_Fc ; 
                delta_Fb_opt = delta_Fb ;
                delta_Fc_opt = delta_Fc ;
                Sn_prime_opt = Sn_prime ;
                S_n_opt = S_n ;
                Sf_life_opt = Sf_life ;
                K_f_opt = K_f ;
                sigma_ea_b = sigma_a*K_f ;
                sigma_a_max = S_a ;
                k_b_opt = k_b ;
                k_c_opt = k_c ;
            end
        end
    end
end

d_opt = d_opt*1000
p_opt = p_opt*1000
d_p_opt = d_p_opt*1000
A_t_opt = A_t_opt*1e6
material_opt_name = Material_class(material_opt)
n = F_i_opt/p_opt*(1/k_b_opt + 1/k_c_opt)
F_i_opt
K_i_opt
sigma_init_b = F_i_opt/A_t_opt/1e6  % K_f olcak mı
sigma_init_c = F_i_opt/A_c/1e6   %
F_e_max
Fb_opt
Fc_opt
delta_sigma_a_b = delta_Fb_opt/A_t_opt/1e6  %
delta_sigma_a_c = delta_Fc_opt/A_c/1e6   %
SF_yield_opt 
F_a_bolt = delta_Fb_opt/2
F_m_bolt = F_i_opt + delta_Fb_opt/2
Sn_prime_opt = Sn_prime_opt/1e6
C_load 
C_size
C_surf
C_rel
S_n_opt = S_n_opt/1e6
Sf_life_opt = Sf_life_opt/1e6
K_f_opt
sigma_a_bolt = delta_Fb_opt/2/A_t_opt/1e6
sigma_ea_b = sigma_ea_b/1e6
sigma_a_max = sigma_a_max/1e6
SF_fatigue_opt
Price_min

