clear all;
clc;
close all;

%Newton-Raphson Method kullanarak denklemin kökünü numerik olarak bulmak
%başlangıç tahmini yapıp onu itere ederek çözüme gider
%fonksiyonun türevini analitik yöntemlerle bulursun.

format long

A=200;
p0=15;
delta=10^(-6);
max1=100;
f=@(x) x.^3-A;
df=@(x) 3*x.^2;

for k=1:max1
    
    y=f(p0(k));
    p1=p0(k)-f(p0(k))/df(p0(k));
    err=abs(p1-p0(k));
    relerr=err/(abs(p1)+delta);
    p0=[p0 p1];
    
    if (err<delta) && (relerr<delta) && (abs(y)<delta)
        break
    end
    
end
p=p0(end)
err
relerr

