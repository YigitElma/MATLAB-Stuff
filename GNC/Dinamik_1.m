% -------------------------------------------------------------------------
% 
% 
% 
% 
% Yiğit Günsür Elmacıoğlu
% 02.09.2021
% 
% Body frame deki hızları bulmak için body frame e göre olan hız değişim değerleri
% bulundu. Ardından numerik integrasyon kullanılarak hızlar, koordinat
% sistemi transformasyonu sonrasında konuma gidildi.
% 
% 
% 
% -------------------------------------------------------------------------
clc
clear all
close all
format longG
 
m = 100;           %mühimmat kütlesi
dt = 0.001;        %iterasyon zaman aralığı
g = 9.807;         %yerçekimi ivmesi
max_time = 100;    %sonsuz döngüyü engellemek için maksimum zaman
time = linspace(0, max_time, max_time/dt);
N = length(time);
C_d = linspace( 0.5, 1, 1000) ;  %rastgele verilmiş CD değerleri
teta_pos = linspace(0, 90, 1000) ;
C_l = linspace( 0, 1, 1000) ;    %rastgele verilmiş CL değerleri
phi_pos = linspace(0, 90, 1000) ;
A = pi * 0.05^2 ;      %sürüklenme ve kaldırma kuyvvetleri için gereken referans yüzey alanı
rho = 1.225;           %hava yoğunluğu
I_yy = 6.225;          %Eylemsizlik momenti
L = 1;                 %CP ve CG arası mesafe
 
V_x_init = 100;    %mühimmat başlangıç hızı (yatay)
V_y_init = 0;      %mühimmat başlangıç hızı (dikey)
 
%başlangıç konumları
x_m_init = -1000;
y_m_init = 1000;
 
x(1) = x_m_init;
y(1) = y_m_init;
 
%yere göre hızlar
V_x(1) = V_x_init;
V_y(1) = V_y_init;
 
%cisim gövde eksenine göre hızlar
v_X(1) = V_x_init;
v_Y(1) = V_y_init;
 
d_teta(1) = 0;
teta(1) = 0;
phi(1) = 0;
i = 1;
 
while y(i) > 0 && i*dt < max_time
    V(i) = sqrt( V_x(i)^2 + V_y(i)^2 );
 
    if phi(i) < 90
        C_d_t(i) = interp1(teta_pos, C_d, abs(phi(i)));
        C_l_t(i) = interp1(phi_pos, C_l, abs(phi(i)));
    else
        C_d_t(i) = interp1(teta_pos, C_d, 180-abs(phi(i)));
        C_l_t(i) = interp1(phi_pos, C_l, 180-abs(phi(i)));
    end
    
    %dönen bir eksen takımı üzerinde hızların dönen eksen üzerindeki
    %değişimini veren ivme değerleri
    a_X(i) = -0.5*rho*V(i)^2 *A*C_d_t(i)*cosd(phi(i))/m - g*sind(teta(i))  - 0.5*rho*V(i)^2*A*C_l_t(i)*sind(phi(i))/m + d_teta(i)*v_Y(i);
    a_Y(i) = 0.5*rho*V(i)^2 *A*C_l_t(i)*cosd(phi(i))/m - g*cosd(teta(i)) -0.5*rho*V(i)^2*A*C_d_t(i)*sind(phi(i))/m - d_teta(i)*v_X(i);
   
    %dönen eksendeki hız değerleri
    v_X(i+1) = v_X(i) + dt*a_X(i) ;
    v_Y(i+1) = v_Y(i) + dt*a_Y(i) ;
    
    %atak açısı (derece)
    phi(i+1) = atan2d( v_Y(i+1),v_X(i+1) );
     
    if phi(i) > 0
        M(i) = 0.5*rho*V(i)^2 *A*C_l_t(i)*L*cosd(phi(i)) + 0.5*rho*V(i)^2 *A*C_d_t(i)*L*sind(phi(i)) ;
    else
        M(i) = - 0.5*rho*V(i)^2 *A*C_l_t(i)*L*cosd(phi(i)) + 0.5*rho*V(i)^2 *A*C_d_t(i)*L*sind(phi(i));
    end
    
    alfa(i) = M(i)/I_yy ;
    d_teta(i+1) = d_teta(i) + alfa(i)*dt ;     %(rad/s)
    teta(i+1) = teta(i) + d_teta(i)*dt*180/pi; %derece
    
    %yer eksen takımına göre hız değerleri
    V_x(i+1) = v_X(i+1)*cosd(teta(i+1)) - v_Y(i+1)*sind(teta(i+1)) ;
    V_y(i+1) = v_X(i+1)*sind(teta(i+1)) + v_Y(i+1)*cosd(teta(i+1)) ;
 
    %konum verileri
    x(i+1) = x(i) + V_x(i)*dt ;
    y(i+1) = y(i) + V_y(i)*dt ;
 
    %hız vektörünün görselleştirilmesi için 
    vector_x(i) = V(i)*cosd(teta(i)+phi(i)) ; 
    vector_y(i) = V(i)*sind(teta(i)+phi(i)) ; 
    
    %roketin duruş pozisyonunun görselleştirilmesi için (50 rastgele bir
    %değer)
    rocket_x(i) = 50*cosd(teta(i));
    rocket_y(i) = 50*sind(teta(i));
 
    i = i+1;
 
end
 
x_max = max(x(1:i)) + 1000
t_total = (i-1) * dt
