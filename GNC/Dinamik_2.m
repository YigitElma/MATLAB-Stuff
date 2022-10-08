% -------------------------------------------------------------------------
% 
% 
% 
% 
% Yiğit Günsür Elmacıoğlu
% 26.08.2021
% 
% Body frame'e göre yazılan ivmeler yer eksenine döürüldü. Bulunan ivmeler
% integre edilip hızlar çekildi. Atak açısı için ters dönüşüm yapıldı.
% Moment kuvveti rastgele inertia ve CP değerleri için çekildi.
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
max_time = 50;     %sonsuz döngüyü engellemek için yazılmıştır, gerektiğinde artırılabilir
time = linspace(0, max_time, max_time/dt);
N = length(time);
C_d = linspace( 0.5, 1, 1000) ;     %değerler rastgele olarak deneme amaçlı verilmiştir
teta_pos = linspace(0, 90, 1000) ;  
C_l = linspace( 0, 1, 1000) ;       %değerler rastgele olarak deneme amaçlı verilmiştir
phi_pos = linspace(0, 90, 1000) ;
A = pi * 0.05^2 ;       %sürüklenme ve kaldırma kuyvvetleri için gereken referans yüzey alanı
rho = 1.225;            %hava yoğunluğu
I_yy = 60.225;          %Eylemsizlik momenti
L = 0.5;                %CP ve CG arası mesafe
 
V_x_init = 100;    %mühimmat başlangıç hızı (yatay)
V_y_init =0;       %mühimmat başlangıç hızı (dikey)
 
%başlangıç konumları
x_m_init = -1000;
y_m_init = 1000;
 
x = zeros(1, N);
y = zeros(1, N);
V_x = zeros(1, N);
V_y = zeros(1, N);
v_X = zeros(1, N);
v_Y = zeros(1, N);
a_X = zeros(1, N);
a_Y = zeros(1, N);
V = zeros(1, N);
 
x(1) = x_m_init;
y(1) = y_m_init;
 
%yer eksenine göre olan hızlar
V_x(1) = V_x_init;
V_y(1) = V_y_init;
 
%Gövde eksenine göre olan hızlar
v_X(1) = V_x_init;
v_Y(1) = V_y_init;
 
d_teta(1) = 0;    %açısal hız (rad/s)
teta(1) = atan2d( V_y_init,V_x_init ); %uçuş yolu açısı
phi(1) = 0;       %atak açısı (radyan)
i = 1;
 
while y(i) > 0
    V(i) = sqrt( v_X(i)^2 + v_Y(i)^2 );
   
    if phi(i) < 90
        C_d_t(i) = interp1(teta_pos, C_d, abs(phi(i)));
    else
        C_d_t(i) = interp1(teta_pos, C_d, 180-abs(phi(i)));
    end
    C_l_t(i) = interp1(phi_pos, C_l, abs(phi(i)));
 
    %gövde ekseninde dekompoze edilmiş ivmeler
    a_X(i) = -0.5*rho*V(i)^2 *A*C_d_t(i)*cosd(phi(i))/m - g*sind(teta(i))  - 0.5*rho*V(i)^2*A*C_l_t(i)*sind(phi(i))/m ;
    a_Y(i) = 0.5*rho*V(i)^2 *A*C_l_t(i)*cosd(phi(i))/m - g*cosd(teta(i)) -0.5*rho*V(i)^2*A*C_d_t(i)*sind(phi(i))/m ;
 
    %dönüşüm matrisiyle yer eksenine göre yazılmış ivmeler
    A_x(i) = a_X(i)*cosd(teta(i)) - a_Y(i)*sind(teta(i)) ;
    A_y(i) = a_X(i)*sind(teta(i)) + a_Y(i)*cosd(teta(i)) ;
    
    %yer eksenine göre hızlar
    V_x(i+1) = V_x(i) + dt*A_x(i) ;
    V_y(i+1) = V_y(i) + dt*A_y(i) ;
    
    %atak açısı bulabilmek için yapılan ters dönüşüm ve gövde eksenindeki
    %hızlar
    v_X(i+1) = V_x(i+1)*cosd(teta(i)) + V_y(i+1)*sind(teta(i)) ;
    v_Y(i+1) = - V_x(i+1)*sind(teta(i)) + V_y(i+1)*cosd(teta(i)) ;
 
    %atak açısı
    phi(i+1) = atan2d( v_Y(i+1),v_X(i+1) );
 
    x(i+1) = x(i) + ( V_x(i+1)+V_x(i) )*dt/2 ;
    y(i+1) = y(i) + ( V_y(i+1)+V_y(i) )*dt/2 ;
 
    if phi(i) > 0
        M(i) = 0.5*rho*V(i)^2 *A*C_l_t(i)*L*cosd(phi(i)) + 0.5*rho*V(i)^2 *A*C_d_t(i)*L*sind(phi(i)) ;
    else
        M(i) = - 0.5*rho*V(i)^2 *A*C_l_t(i)*L*cosd(phi(i)) + 0.5*rho*V(i)^2 *A*C_d_t(i)*L*sind(phi(i));
    end
    alfa(i) = M(i)/I_yy ;
    d_teta(i+1) = d_teta(i) + alfa(i)*dt ;
    teta(i+1) = teta(i) + d_teta(i)*dt*180/pi ;
 
    %grafikler için hız vektörü komponentleri
    vector_x(i) = V(i)*cosd(teta(i)+phi(i)) ; 
    vector_y(i) = V(i)*sind(teta(i)+phi(i)) ; 
    
    %grafikler için roket gövde pozisyonu komponentleri (50 değeri rastgele
    %verilmiştir)
    rocket_x(i) = 50*cosd(teta(i));
    rocket_y(i) = 50*sind(teta(i));
 
    i = i+1;
end
 
