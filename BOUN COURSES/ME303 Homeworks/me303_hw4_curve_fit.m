clc
clear all
close all

x=[1 2 3 4 5 6 7 8];
y=[10 7.78 4.96 2.83 1.49 0.738 0.351 0.164];
N=length(x);
L=12;

%for convenient variable changes I tried to leave x alone, coefficient of x
%will be equal to A, constant term will be B and what is other side of the
%equation will be change of variable y.

%%This section of the code is for f=Cxe^{-Dx} type of fit
X=x;         %since exponential function is not linear, I made some variable 
Y=log(y./x); %changes to make linear fit of log(y/x)=ln(C)-Dx

m11=sum(X.^2); %these are coefficients of matrix needed to calculate  
m12=sum(X);    %linear fit of a data set.
m21=m12;
m22=N;
r1=Y*X';
r2=sum(Y);

M=[m11 m12;m21 m22];
R=[r1;r2];
U=M\R;     % U(1)=-D  ;    U(2)=ln(C)

C=exp(U(2)); %above matrix operation found A and B of linear func=Ax+B
D=-U(1);     %but since we want C and D, I return back to orginal variables
f=@(x) C*x.*exp(-D*x); %instead of f1, i wrote f(x) for estetic reasons

%%This section of the code is for g=L/(1+Ee^{Fx}) type of fit
T=x;           %since g(x) is no way linear, I made some variable changes
Z=log(L./y-1); %to make linear fit of log(L/y-1)=Fx+ln(E)

s11=sum(T.^2);
s12=sum(T);
s21=s12;
s22=N;
w1=Z*T';
w2=sum(Z);

S=[s11 s12; s21 s22];
W=[w1;w2];
O=S\W;   %O(1)=F   ;    O(2)=ln(E)

E=exp(O(2));  %above matrix operation found A and B of linear func=Ax+B
F=O(1);       %but since we want E and F, I return back to orginal variables
g=@(x) L./(1+E*exp(F*x));  %instead of f2, i wrote g(x) for estetic reasons

E_2_f1=sqrt(sum((y-f(x)).^2)/N);  %root mean square error of function f(x)
E_2_f2=sqrt(sum((y-g(x)).^2)/N);  %root mean square error of function g(x)

K=[C D; E F]
Error=[E_2_f1 E_2_f2]

if E_2_f1<E_2_f2   %this if statement will print better option among f and g, which has smaller error
    display('f(x) (a.k.a f1) is better fit with lower error')
else
    display('g(x) (a.k.a f2) is better fit with lower error')
end

a=1:0.001:8;  %this is just for plotting functions continuously
plot(a,f(a),'b',x,y,'r*',a,g(a))
legend('f=Cxe^{-Dx}','Data Points','g=L/(1+Ee^{Fx})')
xlabel('x')
ylabel('y')
title('Data Linearization')
grid on
