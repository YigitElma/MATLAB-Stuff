clc
clear all
close all

%Secant metodu kullanarak bir denklemin kökünü numerik olarak bulma
%Öncelikle iki tane tahminde bulun, denklemin gerçek köküyle çok da yakın olmak zorunda değil,
%sonra bunları secant metot denklemine sokarak itere et.
%Newton-Raphson dan farkı iki değer kullanarak çözüme gitmen. Sebebi de
%türevi hesaplarken numerik olarak türev değerini hesaplaman. İki nokta
%arasında eğim yazarak.

f = @(x) exp(-x) - x;
tol = 10^(-15);
maxI = 100;
p0 = 2*pi + 0.6;    %ilk tahminler
p1 = 2*pi - 0.1;

[root, k, y] = Secant(f,tol,maxI,p0,p1)

fplot(f)
grid on
hold on
scatter (p0,f(p0),'b','filled');
scatter (p1,f(p1),'k','filled');
scatter (root,y,'r','filled');

function [root, k, y] = Secant(f, tol, maxI, p0, p1)
xs = [p0 p1];
ys = [f(p0) f(p1)];

for i=1:maxI
    p2 = p1 - f(p1)/((f(p1) - f(p0)) / (p1-p0)); %secant metodu
    p0 = p1;   %her iterasyonda yeni değerler kullanmak için 
    p1 = p2;
    y = f(p2);  
    xs = [xs p2];
    ys = [ys y];
    
    if abs(y)<tol   %istenilen tolerans değerine gelince bitir
        break
    end
    
end

root = xs(end);
k = i;

end
