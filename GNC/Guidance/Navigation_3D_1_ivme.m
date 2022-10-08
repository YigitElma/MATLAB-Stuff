% -------------------------------------------------------------------------
% 
% 
% 
% 
% Yiğit Günsür Elmacıoğlu
% 10.09.2021
% 
% Proportional Navigation 
% Önceki kodlara ek olarak PN metodunu 3boyuta uyarlar. 2D de PN için
% sadece pitch açısındaki LOS değişimi önemsenirken, burada hem pitch hem
% de yaw daki LOS açı değişimine bakıyoruz. Ayrıca uygulanan ivmeler de
% gösteriliyor.
% 
% 
% 
% 
% -------------------------------------------------------------------------
clc 
clear all
close all
  
dt = 0.001;          %zaman aralığı
time = 1:dt:90;      %toplam maks simülasyon süresi
N = length(time) ;
 
K_1 = 3 ;            %PN katsayısı 3-5 arası olması öneriliyor
attack_range = 15 ;  %hedefe yakın mesafe kalması için gereken mesafe (metre) (rastgele)
total_path = 0 ;     %çarpışmaya kadar mühimmatın havada aldığı yol burada toplanacak
total_maneuver = 0 ; %çarpışmaya kadar mühimmatın yaptığı açı değişiklikleri burada toplanacak (manevra eforunun göstergesi)
limit_a = 20 ;       %mühimmatın dayanabileceği maksimum merkezcil kuvvet (g cinsinden)
 
x_target = 40*time ; % 400*ones(1,N) ; %     %hedef x pozisyonu rastgele seçildi
y_target = 20*cos(2*time) + 30 ; % 100*ones(1,N) ; %     %hedef y pozisyonu rastgele seçildi
z_target =  20*sin(2*time) + 200 ; % 300*ones(1,N) ; % 
% Hedef yörüngesi rastgele seçilir ama belirlenen yörünge mühimmatın
% hızından daha yüksek bir hıza neden olmamalıdır. (Zamana göre türev alıp Pisagordan hız bulunabilir)
 
x_init = -400;      %mühimmat başlangıç x konumu
y_init = 0;         %mühimmat başlangıç y konumu
z_init = 200;
x(1) = x_init ;     %mühimmat konum matrisleri
y(1) = y_init ;     %mühimmat konum matrisleri
z(1) = z_init ;
 
v = 80;         %mühimmat hızı (m/s)
teta(1) = 0;   %mühimmat atılma açısı
yaw(1) = 0;
% teta açısı mühimmatın hız vektörünün yer ile yaptığı açıdır (Flight Path
% Angle, FPA) 
v_x(1) = v*cosd(teta(1))*cosd(yaw(1)) ;     %mühimmat hız vektörleri
v_y(1) = v*cosd(teta(1))*sind(yaw(1)) ;     %mühimmat hız vektörleri
v_z(1) = v*sind(teta(1)) ;
 
for i = 1:N-1
    
    %mühimmat ve hedef arasındaki yatay, dikey ve toplam mesafe
    d_x(i) = x_target(i) - x(i) ;
    d_y(i) = y_target(i) - y(i) ;
    d_z(i) = z_target(i) - z(i) ;
    distance(i) = sqrt( d_x(i)^2 + d_y(i)^2 + d_z(i)^2 ) ;
 
    %mühimmat ve hedefin birbirini görüş açısı (Line of Sight)
    LOS_teta(i) = atan2d( d_z(i) , sqrt( d_x(i)^2 + d_y(i)^2 ) );
    LOS_yaw(i)  = atan2d( d_y(i),d_x(i)  );
 
    %Proportional Navigation için gerekli olan LOS açısı değişim hızı
    if i == 1
        d_LOS_teta(i) = 0 ;
        d_LOS_yaw(i) = 0 ;
    else
        d_LOS_teta(i) = (LOS_teta(i) - LOS_teta(i-1)) / dt ;
        d_LOS_yaw(i) = (LOS_yaw(i) - LOS_yaw(i-1)) / dt ;
    end
 
    %mühimmatın konumları hızın numerik integrali alınarak bulunuyor
    x(i+1) = x(i) + v_x(i)*dt ;
    y(i+1) = y(i) + v_y(i)*dt ;
    z(i+1) = z(i) + v_z(i)*dt ;
 
    %toplam mesafe her zaman aralığında alınan yolların toplamı
    total_path = total_path + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 + (z(i+1)-z(i))^2 ) ;
    
    a_teta(i) = v * K_1 * d_LOS_teta(i) * pi/180 ;
    a_yaw(i)  = v * K_1 * d_LOS_yaw(i) * pi/180 ;
    
    a_x(i) = -a_yaw(i)*sind(yaw(i))*cosd(teta(i)) - a_teta(i)*sind(teta(i))*cosd(yaw(i)) ;
    a_y(i) = a_yaw(i)*cosd(yaw(i))*cosd(teta(i)) - a_teta(i)*sind(teta(i))*sind(yaw(i)) ;
    a_z(i) = a_teta(i)*cosd(teta(i)) ;
 
    v_x(i+1) = v_x(i) + a_x(i)*dt ;     %mühimmat hız vektörleri
    v_y(i+1) = v_y(i) + a_y(i)*dt ;     %mühimmat hız vektörleri
    v_z(i+1) = v_z(i) + a_z(i)*dt ;     %mühimmat hız vektörleri
    
    v_tot(i+1) = sqrt(v_x(i+1)^2 + v_y(i+1)^2 + v_z(i+1)^2) ;
    
    teta(i+1) = atan2d( v_z(i+1) , sqrt( v_x(i+1)^2 + v_y(i+1)^2 ) );
    yaw(i+1)  = atan2d( v_y(i+1),v_x(i+1)  );
 
    %hedefin hız komponentleri
    vt_x(i) = ( x_target(i+1) - x_target(i) ) / dt ;
    vt_y(i) = ( y_target(i+1) - y_target(i) ) / dt ;
    vt_z(i) = ( z_target(i+1) - z_target(i) ) / dt ;
 
    %aradaki mesafe 0.1 metrenin altına düşerse çarpışma olduğu
    %varsayımı yapıldı
    if distance(i) < 0.1
        disp('  SUCCESSFUL COLLISION !!!')
        message = '  SUCCESSFUL COLLISION !!!' ;
        break    
    end
    %PN mühimmatın dayanabileceği maksimum ivmeden fazla bir değer
    %hesaplarsa döngü biter
    if sqrt( a_teta(i)^2 + a_yaw(i)^2 ) > limit_a*9.807
        disp('  ACCELERATION LIMIT HAS BEEN EXCEEDED  !!!')
        message = '  ACCELERATION LIMIT HAS BEEN EXCEEDED  !!!' ;
        break        
    end
 
end

Y = i ;
t_col = Y*dt ;
t = 1 ;

fig = figure ;
fig.WindowState = 'fullscreen' ;

for i = 1:50:Y
    vector_x = [ x(i) x_target(i) ] ;
    vector_y = [ y(i) y_target(i) ] ;
    vector_z = [ z(i) z_target(i) ] ;
    
    plot3(x(1:i),y(1:i),z(1:i),'b')
    hold on 
    plot3(x_target(1:i),y_target(1:i),z_target(1:i),'k')
    hold on 
    scatter3(x(i),y(i),z(i),'b^')
    hold on
    scatter3(x_target(i),y_target(i),z_target(i),'bo')
    quiver3(x(i),y(i),z(i),v_x(i),v_y(i),v_z(i),'r')
    hold on
    quiver3(x_target(i),y_target(i),z_target(i),vt_x(i),vt_y(i),vt_z(i),'r')
    hold on
    quiver3(x(i),y(i),z(i),a_x(i),a_y(i),a_z(i),'k')
    hold on
    grid on
    axis equal
    
    if distance(i) > attack_range
        plot3(vector_x,vector_y,vector_z,'c')
    else
        plot3(vector_x,vector_y,vector_z,'b')
    end
    
    xlim( [min( [x x_target(1:Y)] )-100 max( [x x_target(1:Y)] )+100] )
    ylim( [min( [y y_target(1:Y)] )-100 max( [y y_target(1:Y)] )+100] )
    zlim( [min( [z z_target(1:Y)] )-100 max( [z z_target(1:Y)] )+100] )
    
%     anim(t) = getframe ;
%     t = t+1 ;
%     
    drawnow
    hold off
    
end

% video = VideoWriter('GNC_1');
% video.FrameRate = 30;  
% open(video);
% writeVideo(video, anim);
% close(video);



