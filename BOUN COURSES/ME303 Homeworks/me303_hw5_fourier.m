clc
clear all
close all
tic
T1=2; %for part a, T1=1-(-1)=2
T2=8; %for part b, T2=4-(-4)=8
N=1000;
dx=linspace(-1,1,N);
dt=linspace(-4,4,N);

syms x  %later in the code I'll use int() function which requires syms x command
n=1:5;
m=1:15;

%these are for part a
f=@(x) x.^4; %original function
fa_5=@(x) x.^4*cos(2*pi*n.*x/T1); %later I'll use this for a_n calculation, formula comes from fourier equations
fa_15=@(x) x.^4*cos(2*pi*m.*x/T1); %later I'll use this for a_m calculation, formula comes from fourier equations
%fb_5=@(x) x^4*sin(2*pi*n.*x/T1);  I would use this for calculating b_n's
%but since integral will be 0 don't force computer to calculate
%fb_15=@(x) x^4*sin(2*pi*m.*x/T1); I would use this for calculating b_n's
%but since integral will be 0 don't force computer to calculate
f_fourier_5=zeros(1,N);
f_fourier_15=zeros(1,N);

a_0=int(f,x,-1,1)/T1;
a_n=int(fa_5,x,-1,1)*2/T1; %since a_n has different values for different n's, I make computer calculate them by int() function
a_m=int(fa_15,x,-1,1)*2/T1; %since a_n has different values for different n's, I make computer calculate them by int() function
%b_n=zeros(1,5);  since sin is odd function and x^4 is even function this integral will always be 0
%b_m=zeros(1,15); since sin is odd function and x^4 is even function this integral will always be 0

f_fourier_5(1,:)=a_0;
for j=1:5  %for each term there is sum operation a_n*cos(...)+b_n*sin(...)
    f_fourier_5=f_fourier_5+a_n(j)*cos(2*pi*j.*dx./T1);
end
f_fourier_15(1,:)=a_0;
for k=1:15 %for each term there is sum operation a_m*cos(...)+b_m*sin(...)
    f_fourier_15=f_fourier_15+a_m(k)*cos(2*pi*k.*dx./T1);
end

plot(dx,f_fourier_5,dx,f_fourier_15,dx,f(dx))
grid on
legend('f-fourier-5','f-fourier-15','f=x^{4}')
figure

%these are for part b
g=@(x) x.^3; %original function
%later these will be used by int() to calculate coefficients
%gc_5=@(x) x.^3*cos(2*pi*n.*x/T2); 
%gc_15=@(x) x.^3*cos(2*pi*m.*x/T2);
gd_5=@(x) x.^3*sin(2*pi*n.*x/T2);  
gd_15=@(x) x^3*sin(2*pi*m.*x/T2); 
g_fourier_5=zeros(1,N);
g_fourier_15=zeros(1,N);    
    
c_0=int(g,x,-4,4)/T2;
%c_n=zeros(1,5) int(gc_5,x,-4,4)*2/T2; since x^3 is odd and cos is even this integral il always 0
%c_m=zeros(1,15) int(gc_15,x,-4,4)*2/T2; since x^3 is odd and cos is even this integral il always 0
d_n=int(gd_5,x,-4,4)*2/T2; %coefficients of fourier transform
d_m=int(gd_15,x,-4,4)*2/T2; %coefficients of fourier transform
    
g_fourier_5(1,:)=c_0;
for j=1:5 %for each term there is sum operation c_n*cos(...)+d_n*sin(...)
    g_fourier_5=g_fourier_5+d_n(j)*sin(2*pi*j.*dt./T2);
end
g_fourier_15(1,:)=c_0;
for k=1:15 %for each term there is sum operation c_m*cos(...)+d_m*sin(...)
    g_fourier_15=g_fourier_15+d_m(k)*sin(2*pi*k.*dt./T2);
end

plot(dt,g_fourier_5,dt,g_fourier_15,dt,g(dt))
legend('g-fourier-5','g-fourier-15','g=x^{3}')
toc
