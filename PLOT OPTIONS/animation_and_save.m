clc 
clear all
close all

t = 1;
x = 0:0.1:pi*4;
y = sin(x);

fig = figure   %getframe fonksiyonunda piksel ayarlayabil diye

for i=0:0.01:pi*4   %artış miktarını azaltırsan video yavaşlar
    dot = sin(i);
    plot(x,y,i,dot,'r*')
    anim(t) = getframe(fig, [0 0 550 400]);  %fig'in sol alt köşesinin [0 0] yanından 550'ye 400 piksellik resim alır
    t = t+1;  %anim matrisinin indexi düzgün olsun diye
end

video = VideoWriter('deneme5.video');  %denem4.video diye bir avi dosyası yaratır
video.FrameRate = 30;  %saniye başına 30 kare

%açılan avi dosyasına yazılacak matrisi belirtiyorsun
open(video);
writeVideo(video, anim);
close(video);