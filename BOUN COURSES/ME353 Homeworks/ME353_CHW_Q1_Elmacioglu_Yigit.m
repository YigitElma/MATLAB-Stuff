clc
clear 
close all

V = [7.07 10 12.25 14.14 15.81 17.32 18.71 20 21.21 22.36 23.45];
M = [50 100 150 200 250 300 350 400 450 500 550];
A_test = 25;
mu_air = 1.8*10^(-5) ;
rho_air = 1.225 ;

%Calculation of pi terms are in the PDF file
pi_1 = M./V/A_test^(3/2) ;
pi_2 = mu_air/rho_air./V/A_test^(1/2) ;

%Next three part have some constants, I found them using Curve Fitting App
%of MATLAB.

%This part is for second degree polynomial fit
p1 =   1.798*10^-5; 
p2 =  -6.426*10^-6;
p3 =   7.079*10^-7;
f = @(x) p1*x.^2 + p2*x + p3 
RmsE_poly2 = 1.021*10^(-8)

%This part is for 2 term power fit
a =  2.349*10^-8;
b =  -1;  
g = @(x) a*x.^(b)
RmsE_power = 5.549*10^-11

%This part is for exponential fit
d = 6.259*10^-7;
w = -9.022;
h = @(x) d.*exp(w.*x)
RmsE_exp = 1.777*10^-8

x=0.05:0.000001:0.2;
plot(pi_1,pi_2,'r*',x,f(x),'b',x,g(x),'k',x,h(x),'--r')
grid on
xlabel('pi_1')
ylabel('pi_2')
legend('Data Points','Second Degree Polynomial Fit','Power Fit','Exponential Fit')
    
%In this section of the code, I will find v-rho, v-m, v-A graphs
A = 5:0.01:200;
rho = 1:0.01:10;
m = 10:0.1:600;
%Derivation of the formula is on my notes
final = @(area, density, mass) 1./area * sqrt(mu_air.*mass/a./density );

figure
for i=50:50:550
    plot(A,final(A,rho_air,i))
    grid on
    hold on
end
title('Velocity vs Area')
xlabel('Area')
ylabel('Velocity')
legend('m=50','m=100','m=150','m=200','m=250','m=300','m=350',...
    'm=400','m=450','m=500','m=550')

figure
for i=50:50:550
    plot(rho,final(A_test,rho,i))
    grid on
    hold on
end
title('Velocity vs Density')
xlabel('Density')
ylabel('Velocity')
legend('m=50','m=100','m=150','m=200','m=250','m=300','m=350',...
    'm=400','m=450','m=500','m=550')

figure
plot(m,final(A_test,rho_air,m)) 
grid on
title('Velocity vs Mass')
xlabel('Mass')
ylabel('Velocity')







