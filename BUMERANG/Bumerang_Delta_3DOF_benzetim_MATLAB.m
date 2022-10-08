%Bumerang Roket Takımı
%2021 Teknofest Roket Yarışması
%3DOF Uçuş benzetim kodu
%Yiğit Günsür ELMACIOĞLU
%27.05.2021

%Kendi tasarımımız için hazırlanmıştır

clc
clear all
close all
format longG

load openrocket_data.mat 
%Openrocket ve itki değerleri önceden workspace'e kaydedildi. Burada o
%değerler çağırılıyor. Kodun doğru çalışması için KTR ile teslim edilen
%'openrocket_data.mat' dosyasını bu .m dosyası ile aynı konuma kaydedin.
%Değişken isimlerinin anlamları:
%A_tot = Openrocket'tan gelen toplam ivme değerleri
%A_ver = Openrocket'tan gelen dikey ivme değerleri
%altitude = Openrocket'tan gelen irtifa değerleri
%drag_coef = Openrocket'tan gelen drag coefficient değerleri
%range = Openrocket'tan gelen menzil değerleri
%thrust = Openrocket'tan gelen itki değerleri (yarışma bu değerleri
    %sağlamadığı için kendimiz datayı çekip lineer interpolasyon yaptık.)
%time = Openrocket'tan gelen değerlerin kullandığı zaman matrisi(sadece
    %thrust değerleri için geçerli değil, thrust değerleri 5.29 saniye 
    %boyunca eşit aralıklarla verilmiş)
%V_tot = Openrocket'tan gelen toplam skaler hız değerleri
%V_ver = Openrocket'tan gelen dikey hız değerleri
%Simülasyon apogee'ye kadar olacağından Openrocket verileri sadece apogee
%noktasına kadar alınmıştır
 
h = 1/500;         %iterasyonda kullanılan zaman aralığı (saniye)
t_tot = 180;       %rastgele belirlenmiş bir maksimum simülasyon zamanı (saniye)
                      %(olası bir hatada sonsuz döngüyü engelleme için) 
N = t_tot/h;       %maksimum zaman için oluşacak matris boyutu
t_burnout = 5.29;  %M1545 motoru için yanma süresi (saniye)
M = t_burnout/h;   %burnout'a kadar olan matris boyutu 
ISP = 172.65;      %M1545 motorunun özgül itkisi (saniye)
g = 9.804;         %yer çekimi ivmesi
N_data = length(thrust);   %thrust vektörünün boyutu
t_data = linspace(0, t_burnout ,N_data); %thrust datasını zaman matrisi
t_desired = linspace(0, t_burnout, M);   %thrust datası için istenen zaman matrisi
%iteratif bir metot kullanıldığından her 0.002 saniye için bize thrust
%datası gerekiyor fakat Openrocket verileri bu kurala uymuyor. Bu sebeple, 
%t_data matrisine göre verilmiş data t_desired değerleri için lineer interpolasyona 
%tabi tutuluyor. thrust_ext matrisi interpolasyon sonucu elde edilen itki
%değerlerini kaydediyor.
thrust_ext = interp1(t_data, thrust, t_desired);  
%yine iterasyon sırasında dimension sorunu yaşamamak için T (asıl itki
%değerlerini tutan matris) matrisinde burnout sonrası değerler 0 olarak
%katılıyor.
T = [ thrust_ext zeros(1,N-M) ];  %en son elde edilen asıl itki değerleri (Newton)
%KTR raporunda detaylı anlatıldığı üzere kütle değişimi sabit ISP li bir
%motor için aşağıdaki şekilde bulunur. Bu değer kütle değişim eğrisinin
%türevine eşittir.
m_motor = 4.835;   %yakıt kütlesi (kg)
m_dot = thrust_ext / ( ISP*g );   %kütle değişim hızı (kg/s)
m_coef = m_motor/trapz(t_desired,m_dot);  %motor kütle değişim düzeltme katsayısı
C_d_cfd = [ 0.615, 0.615, 0.588, 0.587, 0.59, 0.6, 0.608 0.61 ]; %sürüklenme katsayısı (CFD ile hesaplandı)
v_cfd = [ 0 25 50 75 125 175 225 275 ];  %CFD datasının hız değerleri (m/s)
A = pi*0.125^2/4;  %ön yüzey alanı (m^2) (12.5cm dış çap gövde için)
ramp_angle = 85;   %atış sırasında rampa açısı (derece)
ramp_length = 6;   %rampa uzunluğu (m)

%kullanılacak değişkenler zaman optimizasyonu için önceden maksimum simülasyon 
%süresine göre yaratıldı. Süre fazla gelse bile değerler başta 0 olduğu
%için sorun yaratmıyor.
a_x = zeros(1,N);         %yatay ivme matrisi(m/s^2)
a_y = zeros(1,N);         %dikey ivme matrisi(m/s^2)
F_x = zeros(1,N);         %yatay kuvvet matrisi(Newton)
F_y = zeros(1,N);         %dikey kuvvet matrisi(Newton)
D = zeros(1,N);           %sürüklenme kuvveti matrisi(Newton)
R = zeros(1,N);           %rampanın rokete uyguladığı tepki kuvveti(Newton)
N_y = zeros(1,N);         %yerin rokete uyguladığı dikey tepki kuvveti(Newton)
x = zeros(1,N);           %menzil matrisi(m)
y = zeros(1,N);           %irtifa matrisi(m)
v_x = zeros(1,N);         %yatay hız matrisi(m/s)
v_y = zeros(1,N);         %dikey hız matrisi(m/s)
Temp = zeros(1,N);        %sıcaklık matrisi(Kelvin)
v_sound = zeros(1,N);     %irtiya bağlı ses hızı matrisi(m/s)
P = zeros(1,N);           %irtifaya bağlı basınç matrisi
rho = zeros(1,N);         %irtifaya bağlı hava yoğunluğu matrisi
teta = zeros(1,N);        %uçuş açısı matrisi(derece)
velocity = zeros(1,N);    %toplam skaler hız matrisi(m/s)
mach_number = zeros(1,N); %mach sayısı matrisi

%başlangıç anı koşulları
v_x(1)=0;  %yatay hız(m/s)
v_y(1)=0;  %dikey hz(m/s)
x(1)=0;    %yatay konum (m)
y(1)=0;    %dikey konum (m)
teta(1) = ramp_angle;    %uçuş yolu açısı (derece)  
m(1) = 27.27;  %roket kütlesi(kg)

for i = 1:N-1
    %rampa tepki kuvveti sadece rampa süresince ve rampaya dik olarak
    %uygulanır. Kod için irtifa rampa yüksekliğinden küçük olduğu zamanlar
    %şeklinde koşullandırıldı.
    if ( y(i) < ramp_length*sind(ramp_angle) )               
        if T(i)*sind(teta(i)) < m(i)*g      
            %ilk birkaç milisaniyelik sürede, itki değeri ağırlıktan düşük,
            %bu iterasyonun direk sonlanmasına sebep olmasın diye yer tepki
            %kuvveti eklendi.
            N_y(i) = m(i)*g - T(i)*sind(teta(i));
        end
        R(i) = ( m(i)*g - N_y(i) )*cosd(ramp_angle);  %rampa tepki kuvveti
    end     %else komutu yok çünkü geri kalan değerler önceden 0 olarak yaratılmıştı
    
    %atmoscoesa fonksiyonu US Standart atmosfer 1976 için basınç,sıcaklık,
    %hava yoğunluğu ve ses hızı verilerini irtifaya bağlı olarak verir.
    [Temp(i), v_sound(i), P(i), rho(i)] = atmoscoesa(y(i));
    velocity(i) = sqrt(v_x(i)^2 + v_y(i)^2);
    mach_number(i) = velocity(i)/v_sound(i);
    C_d(i) = interp1(v_cfd, C_d_cfd, velocity(i)); %anlık C_d değeri
    D(i) = 0.5*rho(i)*(velocity(i)^2)*A*C_d(i);  %Sürüklenme Kuvveti(Newton)
    
    %oluşturulan hareket denklemlerine göre sırasıyla yatay ve dikey kuvvetlerin
    %hesaplanması
    F_x(i) = T(i) * cosd(teta(i)) - D(i) * cosd(teta(i)) - R(i)*sind(ramp_angle);           
    F_y(i) = T(i) * sind(teta(i)) - D(i) * sind(teta(i)) - m(i)*g ...
             + R(i)*cosd(ramp_angle) + N_y(i);  
    a_x(i) = F_x(i) / m(i);
    a_y(i) = F_y(i) / m(i);
    
    %Euler iteratif integrasyon metodu ile hız ve konum verilerinin elde
    %edilmesi
    v_x(i+1) = v_x(i) + a_x(i)*h;
    v_y(i+1) = v_y(i) + a_y(i)*h;
    x(i+1) = x(i) + v_x(i)*h;
    y(i+1) = y(i) + v_y(i)*h;
       
    %eğer yükseklik azalmaya başlarsa apogee'ye ulaşılmıştır ve iterasyon
    %biter
    if y(i+1) < y(i)   
        break
    end    
    %roket rampada olduğu sürece uçuş açısı rampa açısına eşit olmak
    %zorunda! roketin rampada olması irtifanın rampa yüksekliğinden düşük
    %olması ile anlaşılabilir. geri kalan uçuşta roket hız vektörü ile aynı
    %yöne baktığı için arctan kullanılabilir.
    if ( y(i) > ramp_length*sind(ramp_angle) ) 
        teta(i+1) = atand( v_y(i+1)/v_x(i+1) );  
    else
        teta(i+1) = ramp_angle;
    end
    
    %değişen hızlarda kütle azalmasının modellenmesi
    %m_dot değişim hızını verdiği için birimi kg/s, küte değişimi için h
    %yani iterasyonlar arası zaman aralığı ile çarpılır ve kg cinsinden
    %değer bulunur. i<M koşulunun sebebi M'in burnout'a kadar olan
    %iterasyon sayısını temsil etmesi. Kütle değişimindeki fazladan hatayı
    %önlemek için m_dot, m_dot'ın integrali ve gerçek yakıt kütlesinden hesaplanan
    %bir düzeltme katsayısı ile çarpılır.
    if i <= M       
        m(i+1) = m(i) - m_dot(i)*h*m_coef;  
    else 
        m(i+1) = m(i);
    end
end
%iterasyonun kaç defa tekrar edildiği bilgisi süre hesabı için kullanılacak
K = i;

acc = sqrt(a_x.*a_x + a_y.*a_y);      %toplam skaler ivme matrisi(m/s^2)
t_flight = K * h            %apogee'ye kadar geçen zaman (saniye)
apogee = max(y)             %apogee noktası irtifası(m)
burnout_alt = y(M)          %burnout irtifası(m)
max_v = max(velocity)       %maksimum hız(m/s)
max_mach_number = max(mach_number)    %maksimum mach sayısı
max_acc = max(acc)          %maksimum skaler ivme (m/s^2)
hor_disp = x(K)             %yataydaki maksimum menzil(m)

t = linspace(0,t_flight,K); %grafikleri çizdirmede kullanılacak zaman matrisi

%Hatalar  (operocket değerleri doğru kabul edilmiştir)
err_range = abs(range(end)-hor_disp)/range(end)*100     %toplam menzildeki yüzdesel hata
err_alt = abs(altitude(end)-apogee)/altitude(end)*100   %apogee yüksekliğindeki yüzdesel hata
err_max_v = abs(max(V_tot)-max_v)/max(V_tot)*100        %maksimum hızdaki yüzdesel hata
err_t_fligth = abs(time(end)-t_flight)/time(end)*100    %apogee süresindeki yüzdesel hata
err_max_a = abs(max(A_tot)-max_acc)/max(A_tot)*100      %maksimum skaler ivmedeki yüzdesel hata
err_mach = abs(0.77-max_mach_number)/0.77*100           %maksimum mach sayısındaki yüzdesel hata

V_hor = sqrt(V_tot.^2-V_ver.^2);  %Openrocket'in vermediği yatay ivme ve hız verileri
A_hor = sqrt(A_tot.^2-A_ver.^2);
ad = length(V_tot);
for i=1:ad-1
    if V_hor(i+1)<V_hor(i)
        A_hor(i) = -A_hor(i);
    end
end

%yatay hız-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,v_x(1:K),'b',time,V_hor,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Yatay Hız-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Yatay Hız (m/s)')

%roketin yörünge grafikleri (openrocket ve bu kod için)
plot(x(1:K),y(1:K),'b',range,altitude,'r')
grid on
axis equal        
title('Roketin Yörüngesi')
xlabel('X ekseni')
ylabel('Y ekseni')
legend('Kodun verisi','Openrocket verisi')

%Yatay konum-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,x(1:K),'b',time,range,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Menzil-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Menzil (m)')

%dikey konum-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,y(1:K),'b',time,altitude,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('İrtifa-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('İrtifa (m)')

%hız-zaman grafikleri (sadece kod sonuçları)
figure
plot(t,v_x(1:K),'b',t(1:K),v_y(1:K),'k')
hold on
plot(t,velocity(1:K),'r','linewidth',2)
grid on
title('Hız-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Hız (m/s)')
legend('yatay hız','dikey hız','toplam skaler hız')

%dikey hız-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,v_y(1:K),'b',time,V_ver,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Dikey Hız-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Dikey Hız (m/s)')

%Toplam skaler hız-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,velocity(1:K),'b',time,V_tot,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Toplam Skaler Hız-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Toplam Skaler Hız (m/s)')

%yatay ve dikey kuvvet-zaman grafikleri (sadece kod için)
figure
plot(t,F_x(1:K),t,F_y(1:K))
grid on
legend('F_x','F_y')
title('Roket Üzerindeki Kuvvetler-Zaman')
xlabel('Zaman (s)')
ylabel('Kuvvet (N)')

%yatay ivme- zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,a_x(1:K),'b',time,A_hor,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Yatay Yönde İvme-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Yatay İvme (m/s^2)')

%dikey ivme-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,a_y(1:K),'b',time,A_ver,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Dikey ivme-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Dikey İvme (m/s^2)')

%toplam skaler ivme-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,acc(1:K),'b',time,A_tot,'r')
grid on
legend('Kodun verisi','Openrocket verisi')
title('Toplam Skaler İvme-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('İvme (m/s^2)')

%uçuş açısı-zaman grafikleri (openrocket ve bu kod için)
figure
plot(t,teta(1:K))
grid on
legend('teta')
title('Uçuş Yolu Açısı-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Uçuş Yolu Açısı (derece)')

%kütle değişim grafiği
figure
plot(t,m(1:K))
ylim([0 30])
grid on
legend('kütle')
title('Roket kütle-Zaman Grafiği')
xlabel('Zaman (s)')
ylabel('Kütle (kg)')

%hava yoğunluğu grafiği
figure
plot(rho(1:K),y(1:K))
grid on
legend('hava yoğunluğu')
title('Hava Yoğunluğu-İrtifa Grafiği')
xlabel('Hava yoğunluğu (kg/m^3)')
ylabel('İrtifa (m)')

%sıcaklık grafiği
figure
plot(Temp(1:K),y(1:K))
grid on
legend('sıcaklık')
title('Sıcaklık-İrtifa Grafiği')
xlabel('Sıcaklık (Kelvin)')
ylabel('İrtifa (m)')

%basınç-irtifa grafiği
figure
plot(P(1:K),y(1:K))
grid on
legend('basınç')
title('Basınç-İrtifa Grafiği')
xlabel('Basınç (Pascal)')
ylabel('İrtifa (m)')

