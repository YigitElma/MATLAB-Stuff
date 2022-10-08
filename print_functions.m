clc;
clear all;
close all;
%yazdırma çeşitleri
%sprintf - fprintf - disp
for i=1:1:10
    if rem(i,2)==0  %remainder kalanı verir
        fprintf('%d çift sayı\n', i);  %\n ile bitmek zorunda
    else
        x=sprintf('%d tek sayı', i);   % disp fonksiyonunda birden çok değişken bastırmak için
        disp(x)
    end
end

       