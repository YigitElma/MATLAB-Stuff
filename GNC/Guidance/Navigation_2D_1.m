% -------------------------------------------------------------------------
% 
% 
% 
% Yiğit Günsür Elmacıoğlu
% 27.08.2021
% 
% Proportional Navigation 
% Bu kod farklı navigasyon komutlarının çarpışma için gereken toplam
% zamana, çarpışmaya kadar alınan toplam yola ve açı değişimlerine etkisini
% görmek için yazılmıştır. Mühimmata gerekli manevraları yaptırmak için
% gereken kuvvetler göz ardı edilmiştir. Mühimmat her komutu anında yerine
% getirebiliyor gibi düşünülmüştür.
% Takip komutları için birkaç farklı senaryo test edilmiştir. Bunlar;
% -düz Proportional navigation
% -hedefe az mesafe kaldığında K'yı artıran Proportional Navigation
% -hedefe az mesafe kaldığında takip moduna geçen Proportional Navigation
% -düz takip modu
% ayrıca takip moduna geçişteki ani açı değişikliğini önlemek için ayrı bir
% metot denenmiştir.
% Yapılan kabuller her ne kadar gerçekçi olmasa da genel bir Guidance,
% navigation and control kodu için önceden denenmesi gereken bir aşamaydı.
% 
% 
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
total_maneuver = 0 ; %çarpışmaya kadar mühimmatın yaptığı açı değişiklikleri burada toplanacak (manevra eforunun göstergesi)
 
x_target = 10*time;      %hedef x pozisyonu rastgele seçildi
y_target = 10*cos(pi*time) + time*15 + 30;    %hedef y pozisyonu rastgele seçildi
% Hedef yörüngesi rastgele seçilir ama belirlenen yörünge mühimmatın
% hızından daha yüksek bir hıza neden olmamalıdır. (Zamana göre türev alıp Pisagordan hız bulunabilir)
 
x_init = 10;      %mühimmat başlangıç x konumu
y_init = 0;         %mühimmat başlangıç y konumu
x(1) = x_init ;     %mühimmat konum matrisleri
y(1) = y_init ;     %mühimmat konum matrisleri
 
delay = 0/dt ;
 
v = 35;         %mühimmat hızı (m/s)
teta(1:delay+1) = 80;   %mühimmat atılma açısı
d_teta(1:delay+1) = 0;
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
        LOS_rate(i) = (LOS_angle(i) - LOS_angle(i-1)) / dt ;
    end
 
    %mühimmatın konumları hızı numerik integrali alınarak bulunuyor
    x(i+1) = x(i) + v*cosd(teta(i))*dt ;
    y(i+1) = y(i) + v*sind(teta(i))*dt ;
 
    %toplam mesafe her zaman aralığında alınan yolların toplamı
    total_path = total_path + sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 ) ;
    
    %LOS ve mühimmat hız vektörü arasındaki açı
    derece(i) = abs(teta(i)- LOS_angle(i));        
 
    
    %%%Düz Proportional Navigation
    % LOS açısının değişim hızı ve teta açısı değişimini oransal yaparak bu
    % iki açı arasındaki açı değişimini sıfırlamaya çalışır. Açının değişim
    % hızı sıfırlanması arada hep sabit bir açı olması demektir. çok nadir
    % durumlar haricinde bu durum elbet çarpışma doğurur. Istisnalardan bir
    % tanesi hedef müdimmattan hızlı ise açı sabit ve 0ken çarpışma
    % olmaması.
    Mod = 1 ;
    d_teta(i+1+delay) = K_1 * LOS_rate(i) ;
    teta(i+1+delay) = teta(i+delay) + d_teta(i+1+delay)*dt ;
    total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
 
    
%     %%%%Direkt takip modu
%     % Eğer hedefinden hızlıysan, her zamana hedefe doğru dümdüz gidersen
%     % onu kesinlikle yakalarsın. Bu metot teta açısını her zaman LOS
%     % açısına eşitler. Açı değişimleri çok sert olabileceğinden her durumda
%     % optimum olmayabilir ama kesin çarpışma sağlar.
%     Mod = 2 ;
%     teta(i+1) = LOS_angle(i) ;
%     d_teta(i) = abs(teta(i+1) - teta(i)) ;
%     total_maneuver = total_maneuver + abs(d_teta(i)) ;
    
 
%     %%%Hedefe kısa mesafe kaldığında takip moduna geçiş
%     % PN açı değişimi çabuk sıfırlandığında mühimmat ve hedefin uzun süre
%     % paralele yakın uçmasına sebep olabilir. Bunu engellemek için belli
%     % bir mesafeden sonra takip modu açılması durumu
%     Mod = 3 ;
%     if distance(i) > attack_range
%         d_teta(i) = K_1 * LOS_rate(i) ;
%         teta(i+1) = teta(i) + d_teta(i)*dt ;
%         total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
%     else
%         teta(i+1) = LOS_angle(i) ;
%         d_teta(i) = abs(teta(i+1) - teta(i)) ;
%         total_maneuver = total_maneuver + abs(d_teta(i)) ;
%     end 
%     
    
%     %%%%Hedefe kısa mesafe kaldığında takip moduna geçiş
%     % Takip moduna geçiş anında çok sert bir açı değişimi olabilir. Bu
%     % geçişi daha yumuşak yapmak için önce kademeli bir düşüş sonra direkt
%     % takibe geçiş denemesi
%     Mod = 4 ;
%     if distance(i) > attack_range
%         d_teta(i) = K_1 * LOS_rate(i) ;
%         teta(i+1) = teta(i) + d_teta(i)*dt ;
%         total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
%     else
%         if derece(i) > 1
%             if LOS_angle(i) > teta(i)
%                 teta(i+1) = teta(i) + 0.1 ;
%                 total_maneuver = total_maneuver + 0.1 ;
%             else
%                teta(i+1) = teta(i) - 0.1 ;
%                total_maneuver = total_maneuver + 0.1 ;
%             end
%         else
%             teta(i+1) = LOS_angle(i) ;
%             d_teta(i) = abs(teta(i+1) - teta(i)) ;
%             total_maneuver = total_maneuver + abs(d_teta(i)) ;
%         end
%     end
 
 
%     %%%%Hedefe kısa mesafe kaldığında K sabitini artır
%     % K sabiti LOS değişim hızının kaç katının mühimmata yaptırılacağını
%     % belirler. Fiziksel kısıtlardan ötürü 2-5 arasında seçilen bu değer
%     % uçuş boyunca sabit olmasaydı ne olurdu denemesi
%     Mod = 5 ;
%     if distance(i) > attack_range
%         K = K_1 ;
%         d_teta(i) = K * LOS_rate(i) ;
%         teta(i+1) = teta(i) + d_teta(i)*dt ;
%         total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
%     else
%         K = K_1 +3 ;
%         d_teta(i) = K * LOS_rate(i) ;
%         teta(i+1) = teta(i) + d_teta(i)*dt ;
%         total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
%     end
 
%     %%%%LOS'in değişmeyeceği kadar uzak mesafelerde Takip modu ile
%     % başlayıp sonra PN moduna geçme denemesi
%     Mod = 6;
%     if distance(i) > 500
%         teta(i+1) = LOS_angle(i) ;
%         d_teta(i) = abs(teta(i+1) - teta(i)) ;
%         total_maneuver = total_maneuver + abs(d_teta(i)) ;
%     else
%         d_teta(i) = K_1 * LOS_rate(i) ;
%         teta(i+1) = teta(i) + d_teta(i)*dt ;
%         total_maneuver = total_maneuver + abs(d_teta(i))*dt ;
%     end
        
        
 
    v_x(i+1) = v * cosd(teta(i+1)) ;    %görsel simülasyondaki vektörler için mühimmat x hızı
    v_y(i+1) = v * sind(teta(i+1)) ;    %görsel simülasyondaki vektörler için mühimmat y hızı
 
    %hedefin hız komponentleri
    vt_x(i) = ( x_target(i+1) - x_target(i) ) / dt ;
    vt_y(i) = ( y_target(i+1) - y_target(i) ) / dt ;
 
    %aradaki mesafe 0.1 metrenin altına düşerse çarpışma olduğu
    %varsayımı yapıldı
    if distance(i) < 0.1
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

% video = VideoWriter('GNC_3');
% video.FrameRate = 20;  
% open(video);
% writeVideo(video, anim);
% close(video);



















