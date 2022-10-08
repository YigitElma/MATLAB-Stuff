%ME303
%Fourier Series-triangular wave

clc
clear all
close all

N=1001;
L=pi;
x=linspace(-2*L,2*L,N);
a0=0;

f3=zeros(1,N);  %3 terim kullanılarak
for i=1:N
    for n=1:2:5
        f3(1,i)=f3(1,i)+(8/(n^2*L^2))*cos(n*L*x(i))/L;
    end
    f3(1,i)=f3(1,i)+a0;
end

plot(x,f3)
grid on
hold on

f10=zeros(1,N);  %10 terim kullanılarak
for i=1:N
    for n=1:2:19
        f10(1,i)=f10(1,i)+(8/(n^2*L^2))*cos(n*L*x(i))/L;
    end
    f10(1,i)=f10(1,i)+a0;
end

plot(x,f10)
grid on
hold on

f50=zeros(1,N);  %50 terim kullanılarak
for i=1:N
    for n=1:2:99
        f50(1,i)=f50(1,i)+(8/(n^2*L^2))*cos(n*L*x(i))/L;
    end
    f50(1,i)=f50(1,i)+a0;
end

plot(x,f50)
grid on
legend('n=3','n=10','n=50')
title('Fourier series representation of f')
xlabel('x')
ylabel('f')