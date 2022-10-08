% -------------------------------------------------------------------------
% Yiğit Günsür Elmacıoğlu
% 08.09.2021
% 
% Proportional Navigation aslında mühimmatın hızına dik olan ivmeyi
% belirtir. Bu kodda açı değişimleri arasında bir ilişki kurulması yerine
% direkt ivme kullanılarak animasyon oluşturulmuştur. İvme için
% belirlenecek maksimum değer (limit_a) aşıldığı takdirde bomba kendini
% parçalayacağından, durdurma koşulu olarak eklenmiştir. Mühimmatın
% merkezcil ivmesi grafikte siyah çizgi ile gösterilmiştir. 
% -------------------------------------------------------------------------
clc 
clear all
close all
  
dt = 0.001;          %zaman aralığı
time = 1:dt:60;      %toplam maks simülasyon süresi
N = length(time) ;
 
K_1 = 3 ;            %PN katsayısı 3-5 arası olması öneriliyor
attack_range = 15 ;  %hedefe yakın mesafe kalması için gereken mesafe (metre) (rastgele)
total_path = 0 ;     %çarpışmaya kadar mühimmatın havada aldığı yol burada toplanacak
limit_a = 20 ;       %mühimmatın dayanabileceği maksimum merkezcil kuvvet (g cinsinden)
 
x_target = 10*time ; % 5000*ones(1,N) ;      %hedef x pozisyonu rastgele seçildi
y_target = 10*cos(2*time) + time*15 + 30; % 200*ones(1,N) ;    %hedef y pozisyonu rastgele seçildi
% Hedef yörüngesi rastgele seçilir ama belirlenen yörünge mühimmatın
% hızından daha yüksek bir hıza neden olmamalıdır. (Zamana göre türev alıp Pisagordan hız bulunabilir)
 
x_init = -100;      %mühimmat başlangıç x konumu
y_init = 0;         %mühimmat başlangıç y konumu
x(1) = x_init ;     %mühimmat konum matrisleri
y(1) = y_init ;     %mühimmat konum matrisleri
 
v = 30;         %mühimmat hızı (m/s)
teta(1) = 90;   %mühimmat atılma açısı
% teta açısı mühimmatın hız vektörünün yer ile yaptığı açıdır (Flight Path
% Angle, FPA) 
v_x(1) = v*cosd(teta(1)) ;     %mühimmat hız vektörleri
v_y(1) = v*sind(teta(1)) ;     %mühimmat hız vektörleri
 
for i = 1:N-1
    
    %mühimmat ve hedef arasındaki yatay, dikey ve toplam mesafe
    d_x(i) = x_target(i) - x(i) ;
    d_y(i) = y_target(i) - y(i) ;
    distance(i) = sqrt( d_x(i)^2 + d_y(i)^2 ) ;
 
    %mühimmat ve hedefin birbirini görüş açısı (Line of Sight)
    LOS_angle(i) = atan2d( d_y(i),d_x(i) );
 
    %Proportional Navigation için gerekli olan LOS açısı değişim hızı
    if i == 1
        LOS_rate(i) = 0 ;
    else
        LOS_rate(i) = (LOS_angle(i) - LOS_angle(i-1))*pi/180 / dt ;
    end
 
    %mühimmatın konumları hızı numerik integrali alınarak bulunuyor
    x(i+1) = x(i) + v_x(i)*dt ;
    y(i+1) = y(i) + v_y(i)*dt ;
 
    %toplam mesafe her zaman aralığında alınan yolların toplamı
    total_path = total_path + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 ) ;
    
    %%%Düz Proportional Navigation
    % Mühimmatın sahip olması gereken dikey ivme
    a_r(i+1) = v * K_1 * LOS_rate(i) ;
   
    %mühimmatın ivme komponentleri (yer eksenine göre)
    a_x(i) = - a_r(i+1)*sind( teta(i) ) ;
    a_y(i) = a_r(i+1)*cosd( teta(i) ) ;
    
    v_x(i+1) = v_x(i) + a_x(i)*dt ;    %görsel simülasyondaki vektörler için mühimmat x hızı
    v_y(i+1) = v_y(i) + a_y(i)*dt ;    %görsel simülasyondaki vektörler için mühimmat y hızı
 
    teta(i+1) = atan2d( v_y(i+1),v_x(i) );
        
    %hedefin hız komponentleri
    vt_x(i) = ( x_target(i+1) - x_target(i) ) / dt ;
    vt_y(i) = ( y_target(i+1) - y_target(i) ) / dt ;
    
    %hedefin ivmesi
    if i == 1
        at(i) = 0;
    else 
        at_x(i) = ( vt_x(i)-vt_x(i-1) ) / dt;
        at_y(i) = ( vt_y(i)-vt_y(i-1) ) / dt;
        at(i) = sqrt( at_x(i)^2 + at_y(i)^2 );
    end
 
    %aradaki mesafe 0.1 metrenin altına düşerse çarpışma olduğu
    %varsayımı yapıldı
    if distance(i) < 0.1
        disp(' SUCCESFULL COLLISION !!!')
        message = ' SUCCESFULL COLLISION !!!' ;
        break    
    end
    %PN mühimmatın dayanabileceği maksimum ivmeden fazla bir değer
    %hesaplarsa döngü biter
    if a_r(i+1) > limit_a*9.807
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
    
    plot(x(1:i),y(1:i),'k')
    hold on 
    plot(x_target(1:i),y_target(1:i),'k')
    hold on 
    scatter(x(i),y(i),'b^')
    hold on
    scatter(x_target(i),y_target(i),'bo')
    quiver(x(i),y(i),v_x(i),v_y(i),'r')
    hold on
    quiver(x_target(i),y_target(i),vt_x(i),vt_y(i),'r')
    hold on
    quiver(x(i),y(i),a_x(i),a_y(i),'k')
    hold on
    grid on
    axis equal
    
    if distance(i) > attack_range
        plot(vector_x,vector_y,'c')
    else
        plot(vector_x,vector_y,'b')
    end
    
    xlim( [min( [x x_target(1:Y)] )-100 max( [x x_target(1:Y)] )+100] )
    ylim( [min( [y y_target(1:Y)] )-100 max( [y y_target(1:Y)] )+100] )
    
    anim(t) = getframe ;
    t = t+1 ;
    
    hold off
    
end

% video = VideoWriter('GNC_1');
% video.FrameRate = 30;  
% open(video);
% writeVideo(video, anim);
% close(video);





