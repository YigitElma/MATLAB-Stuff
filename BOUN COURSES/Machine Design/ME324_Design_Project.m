% 2020-2021-2
% ME324 Design Project
% Yiğit Günsür Elmacıoğlu  &   Ege Oray
%       2017405120         &  2017405057

clc
clear all
close all
%format longG

a = 0;
b = 7;
g = 9.81;            %gravitational acceleration (m/s^2)
Kf = 1.3;            %fatigue stress concentration factor
SF = 2;              %safety factor
p_in = 10^5;         %inlet pressure (Pascal)
p_out = 2*10^6;      %outlet pressure (Pascal)
V_out = (2200+50*a)*10^-6;  %delivered air volume (m^3)
V_in = p_out * V_out / p_in;
eff = 0.5;   %efficiency of the compressor
C = 0.005;   %speed fluctuation coefficient

c_ratio = [ 10 11 14 ];  %respectively AISI-1015 / 1040 / 1095
S_u = 10^6*[ 420.6 620.5 965.3 ];  %respective UTS (Pascal)
S_y = 10^6*[ 313.7 413.7 572.3 ];  %respective Yield stresses (Pascal)
c_surf =[ 0.8 0.77 0.72 ];         %respective surface factors 

W = p_out * V_out * log( V_in / V_out );  %energy required
W_n = W / eff;   %energy needed for motor to supply net energy of W

r_sg = (11+a) / 200;  %radius of small gear
rpm_c = 300 + b*8;    %average RPM of the compressor
T_ave = W_n/(2*pi);   %constant torque exerted to the larger gear
m_f_min = W_n/2/(100^2)/C;   %minimum mass of the flywheel
% m_f_min =270;
% When we plugged a bit higher mass (minimum is 263kg), the minimum diameter
% decreased but the overall cost increased. After multiple trials, we see
% that optimum flywheel mass which minimizes the cost is minimum possible
% flywheel mass as explained in the handwritten part of the report.
n_max = 10;   %maximum value of gear ration
d_min = [ 0.3 0.3 0.3 ];    %this only for comparison purposes
L_max = (50 + 2*b)/100 ;    %max allowable distance between point A and pulley 
L_d = (65 + 2*b)/100;       %distance between points A and D
L_e = L_d +0.14 ;           %distance between points A and E

for n = 1:n_max   
    l = 1;
    for L = 0.3:0.02:L_max
        %according to derived equilibrium equations
        %unknown forces can be found by these equations
        F_z = -2*T_ave/r_sg/n;   %tangential force of the small gear
        F_y = F_z*tand(20);      %radial force of the small gear
        D_y = ( L_e*g*m_f_min - F_y*0.15 - 200*L )/L_d ; %point D reaction
        D_z = -F_z*0.15/L_d ;    %point D reaction
        A_z = -F_z -D_z ;        %point A reaction
        A_y = -D_y -F_y -200 +g*m_f_min ;  %point A reaction
        T_f = -T_ave/n - r_sg*F_z ;   %flywheel torque at compression stroke
        r_f_max = sqrt( pi*T_f / ( C*10*m_f_min*(rpm_c*pi/30)^2 ) );   %maximum radius of the flywheel
        
        h = 0.0001;     %delta x for moment-position graphs
        
        %bending moment between bearing A and small gear
        x1 = 0:h:0.15;  
        M1 = sqrt( A_z^2 + A_y^2 )*x1 ;  

        %bending moment between small gear and pulley
        x2 = 0.15+h:h:L;
        M2 = sqrt( (A_z*x2 + F_z*(x2-0.15) ).^2 +  (A_y*x2 + F_y*(x2-0.15)  ).^2 );

        %bending moment between pulley and bearing D
        x3 = L+h:h:L_d ;
        M3 = sqrt( (A_z*x3 + F_z*(x3-0.15) ).^2 + ...
            (A_y*x3 + F_y*(x3-0.15) + 200*(x3-L)  ).^2 );

        %bending moment between bearing D and flywheel
        x4 = L_d+h:h:L_e;
        M4 = sqrt( (A_z*x4 + F_z*(x4-0.15) + D_z*(x4-L_d) ).^2 + ...
            (A_y*x4 + F_y*(x4-0.15) + 200*(x4-L) + D_y*(x4-L_d)  ).^2 );

        %combined bending moment graph
        M_tot = [ M1 M2 M3 M4 ];
%         x5 = 0:h:L_e;
%         plot(x5,M_tot)  %effects of increasing n and L can be seen
%         hold all
%         title('Moment Diagram')
%         xlabel('distance (m)')
%         ylabel('Bending Moment (Nm)')
%         grid on
        
        %for each n value, we try several L values, for each L value we
        %find peak bending moment
        max_m(l) = max(M_tot);
        %among these peak values (with same n value), we pick the distance 
        %L which minimizes the peak bending moment
        min_m = min(max_m);
        %matrix index of the corresponding moment, which will be needed to
        %find L value
        i_min = find(max_m==min(max_m));
        
        for i = 1:3      %for each material
            k(i) = 1;    %to count iteration number
            d = 0.25;    %initial diameter guess
            while 1  %we use an infinite loop
                %mean and alternating bending stresses
                sigma_m = 0;
                sigma_a = Kf * min_m *32/pi/d^3;

                %mean and alternating shear stresses
                shear_m = Kf * r_sg*F_z *8/pi/d^3;
                shear_a = shear_m;

                c_size = 1.189*(d*1000)^(-0.097);
                S_n = c_size * 1 * 0.753 * c_surf(i) * 0.5 * S_u(i) ;
                
                %equivalent alternating and mean stresses
                eq_a = sqrt( sigma_a^2 + 3*shear_a^2 );
                eq_m = shear_m;
                   
                %since we started from a high diameter value, initial Safety 
                %Factor is higher than 2, so as a break condition, we use Safety
                %Factor smaller than 2 
                if eq_a/S_n + eq_m/S_u(i) >= 1/SF 
                    material = i;  %chosen material
                    L_opt(i) = 0.3 + (i_min-1)*0.02;  %optimum pulley position for each material
                    n_opt(i) = n;  %optimum gear ratio for each material
                    c_eq_a(i) = eq_a;  %corresponding equivalent alternating stresses
                    c_eq_m(i) = eq_m;  %corresponding equivalent mean stresses
                    d_min(i) = d;
                    sig_m(i) = sigma_m;
                    sig_a(i) = sigma_a;
                    sh_a(i) = shear_a;
                    sh_m(i) = shear_m;
                    break
                else
                    d = d - 0.0001; %by decreasing diameter we find minimum
                                    %diameter which has SF=2
                end
               
                if d < 0 %to prevent infinite loop
                    break
                end
                
                if d < d_min(i)
                    d_min(i) = d;  %for each material, minimum allowable diameters
                end
                k(i) = k(i)+1;  %iteration number
            end
        end
        l = l+1;  %to find optimum pulley position
    end
end
n_g = rpm_c;
n_p = n_opt(material)*rpm_c;
power = W_n * rpm_c / 60;
E_motor = W_n/2;
E_fly = W_n/2;
V_ave = r_f_max*n_p*pi/30;
T_m = T_ave/n;
C_load =1;
C_rel = 0.753;
Kfs = Kf;
y_f = S_y./( c_eq_a + c_eq_m );
cost = pi*0.93*7800*0.25*c_ratio.*d_min.^2 + m_f_min;

%%%%%%%%%%
    %Finding compression and idle bending moments as well as torsional moments
    %in both z and y directions for chosen configuration
    L = L_opt(material);
    n = n_opt(material);
    j = 1;
    for F2_z = [0 -2*T_ave/r_sg/n]
        F2_y = F2_z*tand(20);
        D2_y = ( L_e*g*m_f_min - F2_y*0.15 - 200*L )/L_d ;
        D2_z = -F2_z*0.15/L_d ;
        A2_z = -F2_z -D2_z ;
        A2_y = -D2_y -F2_y -200 +g*m_f_min ;
        T2_f = -T_ave/n - r_sg*F2_z ;

        %between bearing A and small gear
        xx1 = [0 0.15];  
        MM1y(j,:) = A2_z*xx1 ;
        MM2z(j,:) = A2_y*xx1 ;

        %between pulley and bearing D
        xx3 = [L L_d] ;
        MM3y(j,:) = A2_z*xx3 + F2_z*(xx3-0.15) ;
        MM4z(j,:) = A2_y*xx3 + F2_y*(xx3-0.15) + 200*(xx3-L) ;

        j = 2;
    end

    F3_z = 0;
    F3_y = F3_z*tand(20);
    D3_y = ( L_e*g*m_f_min - F3_y*0.15 - 200*L )/L_d ;
    D3_z = -F3_z*0.15/L_d ;
    A3_z = -F3_z -D3_z ;
    A3_y = -D3_y -F3_y -200 +g*m_f_min ;
    T3_f = -T_ave/n - r_sg*F3_z ;
    
    A_z_idle = A3_z;
    D_z_idle = D3_z;
    A_y_idle = A3_y;
    D_y_idle = D3_y;
    %bending moments during compression stroke in y direction
    M_Ay_c = MM1y(2,1);
    M_By_c = MM1y(2,2);
    M_pulley_y_c = MM3y(2,1);
    M_Dy_c = MM3y(2,2);
    M_fy_c = 0;
    %bending moments during idle stroke in y direction
    M_Ay_i = MM1y(1,1);
    M_By_i = MM1y(1,2);
    M_pulley_y_i = MM3y(1,1);
    M_Dy_i = MM3y(1,2);
    M_fy_i = 0;
    %bending moments during compression stroke in z direction
    M_Az_c = MM2z(2,1);
    M_Bz_c = MM2z(2,2);
    M_pulley_z_c = MM4z(2,1);
    M_Dz_c = MM4z(2,2);
    M_fz_c = 0;
    %bending moments during idle stroke in z direction
    M_Az_i = MM2z(1,1);
    M_Bz_i = MM2z(1,2);
    M_pulley_z_i = MM4z(1,1);
    M_Dz_i = MM4z(1,2);
    M_fz_i = 0;
    %torsion moments during compression stroke
    T_A_c = 0;
    T_B_c = 2*T_f;
    T_pulley_c = 2*T_f;
    T_D_c = T_f;
    T_f_c = T_f;
    %torsion moments during idle stroke
    T_A_i = 0;
    T_B_i = 0;
    T_pulley_i = T_f;
    T_D_i = T_f;
    T_f_i = T_f;

%%%%%%%%

sprintf('a = %d', a)
sprintf('b = %d', b)
sprintf('Optimum gear ratio is %d', n_opt(material))
sprintf('optimum pulley position is %f m', L_opt(material))
sprintf('Rotational speed of the larger gear is %d RPM', n_g )
sprintf('Rotational speed of the pulley is %d RPM', n_p )
sprintf('Required power of the motor is %f kW', power/1000 )
sprintf('Energy supplied by the motor during compression is %f kJ', E_motor/1000 )
sprintf('Energy supplied by the flywheel during compression is %f kJ', E_fly/1000 )
sprintf('Mass of the flywheel is %f kg', m_f_min)
sprintf('Radius of the flywheel rim is %f m', r_f_max)
sprintf('Average speed of the rim is %f m/s', V_ave )
sprintf('Torque of the flywheel during compression is %f Nm', T_f)
sprintf('Torque of the motor is %f Nm', T_m )
sprintf('Tangential force on the gears is %f N', F_z)
sprintf('Radial force on the gears is %f N', F_y)

sprintf('Reaction of bearing A in y direction during compression is %f N', A_y)
sprintf('Reaction of bearing A in z direction during compression is %f N', A_z)
sprintf('Reaction of bearing D in y direction during compression is %f N', D_y)
sprintf('Reaction of bearing D in z direction during compression is %f N', D_z)

sprintf('Reaction of bearing A in y direction during idle is %f N', A_y_idle)
sprintf('Reaction of bearing A in z direction during idle is %f N', A_z_idle)
sprintf('Reaction of bearing D in y direction during idle is %f N', D_y_idle)
sprintf('Reaction of bearing D in z direction during idle is %f N', D_z_idle)

sprintf('Bending moment on A in y direction during compression is %f Nm', M_Ay_c)
sprintf('Bending moment on gear B in y direction during compression is %f Nm', M_By_c)
sprintf('Bending moment on pulley in y direction during compression is %f Nm', M_pulley_y_c)
sprintf('Bending moment on D in y direction during compression is %f Nm', M_Dy_c)
sprintf('Bending moment on flywheel in y direction during compression is %f Nm', M_fy_c)

sprintf('Bending moment on A in y direction during idle is %f Nm', M_Ay_i)
sprintf('Bending moment on gear B in y direction during idle is %f Nm', M_By_i)
sprintf('Bending moment on pulley in y direction during idle is %f Nm', M_pulley_y_i)
sprintf('Bending moment on D in y direction during idle is %f Nm', M_Dy_i)
sprintf('Bending moment on flywheel in y direction during idle is %f Nm', M_fy_i)

sprintf('Bending moment on A in z direction during compression is %f Nm', M_Az_c)
sprintf('Bending moment on gear B in z direction during compression is %f Nm', M_Bz_c)
sprintf('Bending moment on pulley in z direction during compression is %f Nm', M_pulley_z_c)
sprintf('Bending moment on D in z direction during compression is %f Nm', M_Dz_c)
sprintf('Bending moment on flywheel in z direction during compression is %f Nm', M_fz_c)

sprintf('Bending moment on A in z direction during idle is %f Nm', M_Az_i)
sprintf('Bending moment on gear B in z direction during idle is %f Nm', M_Bz_i)
sprintf('Bending moment on pulley in z direction during idle is %f Nm', M_pulley_z_i)
sprintf('Bending moment on D in z direction during idle is %f Nm', M_Dz_i)
sprintf('Bending moment on flywheel in z direction during idle is %f Nm', M_fz_i)

sprintf('Torsional moment on A during compression %f Nm', T_A_c)
sprintf('Torsional moment on gear B during compression %f Nm', T_B_c)
sprintf('Torsional moment on pulley during compression %f Nm', T_pulley_c)
sprintf('Torsional moment on D during compression %f Nm', T_D_c)
sprintf('Torsional moment on flywheel during compression %f Nm', T_f_c)

sprintf('Torsional moment on A during idle %f Nm', T_A_i)
sprintf('Torsional moment on gear B during idle %f Nm', T_B_i)
sprintf('Torsional moment on pulley during idle %f Nm', T_pulley_i)
sprintf('Torsional moment on D during idle %f Nm', T_D_i)
sprintf('Torsional moment on flywheel during idle %f Nm', T_f_i)

sprintf('Chosen material is number %d which is AISI-1095', material )
sprintf('Corresponding UTS is %f MPa', S_u(material)*10^-6 )
sprintf('Corresponding yield strenght is %f MPa', S_y(material)*10^-6 )
sprintf('Corresponding hardness is %f', S_u(material)/(10^6)/3.45 )
sprintf('Required diameter of the shaft with chosen material is %f cm', d_min(material)*100 )
sprintf('Corresponding surface factor is %f', c_surf(material) )
sprintf('Corresponding size factor is %f', c_size )
sprintf('Corresponding load factor is %d', C_load )
sprintf('Corresponding reliability factor is %f', C_rel )
sprintf('Corresponding S_n is %f MPa', S_n*10^-6 )
sprintf('The most critical point is Point B')
sprintf('Fatigue stress concentration of the keyway is %f for bending', Kf)
sprintf('Fatigue stress concentration of the keyway is %f for torsion', Kfs )
sprintf('Alternating bending stress is %f MPa', sig_a(material)*10^-6)
sprintf('Alternating shear stress is %f MPa', sh_a(material)*10^-6)
sprintf('Mean bending stress is %f MPa', sig_m(material)*10^-6)
sprintf('Mean shear stress is %f MPa', sh_m(material)*10^-6)
sprintf('Equivalent alternating stress is %f MPa', c_eq_a(material)*10^-6)
sprintf('Equivalent mean stress is %f MPa', c_eq_m(material)*10^-6)
sprintf('Factor for local yielding of the chosen design is %f', y_f(material) )
sprintf('Safety factor against fatigue failure is %d', SF)




