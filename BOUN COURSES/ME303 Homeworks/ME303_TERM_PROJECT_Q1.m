%Yiğit Günsür ELMACIOĞLU
%2017405120
%18 JANUARY 2021

%ME303 TERM PROJECT 2020 FALL
%QUESTION 1
%MODELLING THE SPREAD OF A VIRUS

clc
clear all
close all

%Given quantities
ti=0;         %first day
tf=100;       %last day
N_bed=1000;   %total number of beds
h=1/500;
N= (tf-ti)/h;
t=ti:h:tf;
c=4;      %encounters per day
b=0.2;    %transmission probability per encounter
a=0.125;  %infected to infectious
y=0.1;    %infected to symptomatic
w=0.2;    % 1 / recovery time

%I will use these matrices to store Runge-Kutta results
S=zeros(1,N+1);  %Susceptibles
E=zeros(1,N+1);  %Exposed
I=zeros(1,N+1);  %Infected
M=zeros(1,N+1);  %Medically Symptomatic
R=zeros(1,N+1);  %Recovered
%Initial Conditions
S(1)=10000;
E(1)=10;
I(1)=0;
M(1)=0;
R(1)=0;

%I will use these matrices to store Euler results
S_e=zeros(1,N+1);
E_e=zeros(1,N+1);
I_e=zeros(1,N+1);
M_e=zeros(1,N+1);
R_e=zeros(1,N+1);
%Initial Conditions
S_e(1)=10000;
E_e(1)=10;
I_e(1)=0;
M_e(1)=0;
R_e(1)=0;

%Differantial Equations that relate S,E,I,M and R
%I will use them for both Runge-Kutta and Euler Methods
d_S = @(S,E,I,M,R) -c*b*I / (S+E+I+M+R) * S;
d_E = @(S,E,I,M,R)  c*b*I / (S+E+I+M+R) * S - a*E;
d_I = @(S,E,I,M,R)  a*E - y*I;
d_M = @(S,E,I,M,R)  y*I -  w*M;
d_R = @(S,E,I,M,R)  w*M;
%tic;
for i=1:N
    % EULER PART 
    d_S_e= d_S( S_e(i),E_e(i),I_e(i),M_e(i),R_e(i) );
    d_E_e= d_E( S_e(i),E_e(i),I_e(i),M_e(i),R_e(i) );
    d_I_e= d_I( S_e(i),E_e(i),I_e(i),M_e(i),R_e(i) );
    d_M_e= d_M( S_e(i),E_e(i),I_e(i),M_e(i),R_e(i) );
    d_R_e= d_R( S_e(i),E_e(i),I_e(i),M_e(i),R_e(i) );
    
    S_e(i+1)= S_e(i) + d_S_e*h;
    E_e(i+1)= E_e(i) + d_E_e*h;
    I_e(i+1)= I_e(i) + d_I_e*h;
    M_e(i+1)= M_e(i) + d_M_e*h;
    R_e(i+1)= R_e(i) + d_R_e*h; 
    
    % RUNGE-KUTTA PART
    f1_S= d_S( S(i),E(i),I(i),M(i),R(i) );
    f1_E= d_E( S(i),E(i),I(i),M(i),R(i) );
    f1_I= d_I( S(i),E(i),I(i),M(i),R(i) );
    f1_M= d_M( S(i),E(i),I(i),M(i),R(i) );
    f1_R= d_R( S(i),E(i),I(i),M(i),R(i) );
        
    f2_S= d_S( S(i)+h/2*f1_S ,E(i)+h/2*f1_E, I(i)+h/2*f1_I, M(i)+h/2*f1_M, R(i)+h/2*f1_R );
    f2_E= d_E( S(i)+h/2*f1_S ,E(i)+h/2*f1_E, I(i)+h/2*f1_I, M(i)+h/2*f1_M, R(i)+h/2*f1_R );
    f2_I= d_I( S(i)+h/2*f1_S ,E(i)+h/2*f1_E, I(i)+h/2*f1_I, M(i)+h/2*f1_M, R(i)+h/2*f1_R );
    f2_M= d_M( S(i)+h/2*f1_S ,E(i)+h/2*f1_E, I(i)+h/2*f1_I, M(i)+h/2*f1_M, R(i)+h/2*f1_R );
    f2_R= d_R( S(i)+h/2*f1_S ,E(i)+h/2*f1_E, I(i)+h/2*f1_I, M(i)+h/2*f1_M, R(i)+h/2*f1_R );
    
    f3_S= d_S( S(i)+h/2*f2_S ,E(i)+h/2*f2_E, I(i)+h/2*f2_I, M(i)+h/2*f2_M, R(i)+h/2*f2_R );
    f3_E= d_E( S(i)+h/2*f2_S ,E(i)+h/2*f2_E, I(i)+h/2*f2_I, M(i)+h/2*f2_M, R(i)+h/2*f2_R );
    f3_I= d_I( S(i)+h/2*f2_S ,E(i)+h/2*f2_E, I(i)+h/2*f2_I, M(i)+h/2*f2_M, R(i)+h/2*f2_R );
    f3_M= d_M( S(i)+h/2*f2_S ,E(i)+h/2*f2_E, I(i)+h/2*f2_I, M(i)+h/2*f2_M, R(i)+h/2*f2_R );
    f3_R= d_R( S(i)+h/2*f2_S ,E(i)+h/2*f2_E, I(i)+h/2*f2_I, M(i)+h/2*f2_M, R(i)+h/2*f2_R );
    
    f4_S= d_S( S(i)+h*f3_S ,E(i)+h*f3_E, I(i)+h*f3_I, M(i)+h*f3_M, R(i)+h*f3_R );
    f4_E= d_E( S(i)+h*f3_S ,E(i)+h*f3_E, I(i)+h*f3_I, M(i)+h*f3_M, R(i)+h*f3_R );
    f4_I= d_I( S(i)+h*f3_S ,E(i)+h*f3_E, I(i)+h*f3_I, M(i)+h*f3_M, R(i)+h*f3_R );
    f4_M= d_M( S(i)+h*f3_S ,E(i)+h*f3_E, I(i)+h*f3_I, M(i)+h*f3_M, R(i)+h*f3_R );
    f4_R= d_R( S(i)+h*f3_S ,E(i)+h*f3_E, I(i)+h*f3_I, M(i)+h*f3_M, R(i)+h*f3_R );
    
    S(i+1)= S(i) + h/6 * (f1_S + 2*f2_S + 2*f3_S + f4_S );
    E(i+1)= E(i) + h/6 * (f1_E + 2*f2_E + 2*f3_E + f4_E );
    I(i+1)= I(i) + h/6 * (f1_I + 2*f2_I + 2*f3_I + f4_I );
    M(i+1)= M(i) + h/6 * (f1_M + 2*f2_M + 2*f3_M + f4_M );
    R(i+1)= R(i) + h/6 * (f1_R + 2*f2_R + 2*f3_R + f4_R );       
    
end
% T=M-M_e;      %I used this part to see differences between Euler and Runge-Kutta
% a=max(abs(T)) %In the report, I found by this part that there is
                %negligble difference between RK and Euler.

% time_e=toc   %Commenting Euler and RK respectively this toc gives me
               %their time consumption

% [m_m,d_m]=max(M);    %This part gave me the peak values of M and I also                      
% [m_i,d_i]=max(I);    %their time of occurance
% m_m  %max of M
% d_m=d_m/(N+1)*100    %day of the peak
% m_i  %max of I
% d_i=d_i/(N+1)*100    %day of the peak

%Plotting Euler Solutions 
plot(t,S_e)
hold on
plot(t,E_e)
hold on
plot(t,I_e)
hold on
plot(t,M_e)
hold on
plot(t,R_e)
hold on
plot([0 100],[N_bed N_bed],'r--','Linewidth',2)
legend('Susceptibles (S)','Exposed (E)','Infected (I)','Medically Symptomatic (M)','Recovered (R)','Available Beds (AB)')
title('Euler Solution for given conditions and h=1/500')
xlabel('Time(days)')
ylabel('Population')

%Plotting Runge-Kutta Solution
figure
plot(t,S)
hold on
plot(t,E)
hold on
plot(t,I)
hold on
plot(t,M)
hold on
plot(t,R)
hold on
plot([0 100],[N_bed N_bed],'r--','Linewidth',2)
legend('Susceptibles (S)','Exposed (E)','Infected (I)','Medically Symptomatic (M)','Recovered (R)','Available Beds (AB)')
title('Runge-Kutta Solution for given conditions and h=1/500')
xlabel('Time(days)')
ylabel('Population')


%Plotting ODE45 Solution
figure
%tic
[t,Y]=ode45(@odefun,[0 100],[10000; 10; 0; 0; 0]);
%time_ode45=toc
plot(t,Y(:,1))
hold on
plot(t,Y(:,2))
hold on
plot(t,Y(:,3))
hold on
plot(t,Y(:,4))
hold on
plot(t,Y(:,5))
plot([0 100],[N_bed N_bed],'r--','Linewidth',2)
legend('Susceptibles (S)','Exposed (E)','Infected (I)','Medically Symptomatic (M)','Recovered (R)','Available Beds (AB)')
title('ODE45 Solution for given conditions')
xlabel('Time(days)')
ylabel('Population')

function Y=odefun(t,K)
    c=4;       %constants cannot be transfered to ODE45 directly
    b=0.2;     %this way was more practical for me
    a=0.125; 
    y=0.1;    
    w=0.2;    

    Y(1,1) = -c*b*K(3) / (K(1)+K(2)+K(3)+K(4)+K(5)) * K(1);
    Y(2,1) =  c*b*K(3) / (K(1)+K(2)+K(3)+K(4)+K(5)) * K(1) - a*K(2);
    Y(3,1) =  a*K(2) - y*K(3);
    Y(4,1) =  y*K(3) -  w*K(4);
    Y(5,1) =  w*K(4);
end
