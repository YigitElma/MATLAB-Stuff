clc
clear all
tic
%From textbook
density=995.7; %kg/m^3-water @ 30C (Table B.2)
viscosity=7.975*10^-4; %kg/m.s-water @ 30C (Table B.2)
roughness=0.15; %mm-galvanized iron (Table 8.1)
%----
%From question
L=30; %m-pipe length
D=40*10^-2; %m-pipe diameter
area=pi*D^2*0.25; %m^2-pipe cross-section
%----
%Solution (a)
k=1;
for m=15:1:30 %kg/s-mass flowrate
 V=m/(density*area); %m/s
 Re=density*V*D/viscosity; %>4000 --> Turbulent Flow
 rel=roughness*10^-3/D; %relative roughness
 f_old=1;
 for i=1:10^6
 a=-2*log10(rel/3.7+2.51/(Re*sqrt(f_old)));
 f_new=(1/a)^2;
 error=abs(f_new-f_old)/f_new;
 if error<10^-1
 f=f_new;
 break
 else
 f_old=f_new;
 end
 end
 h_L_a(k)=f*(L/D)*(V^2/(2*9.81)); %m
 mass_a(k)=m; %kg/s
 k=k+1;
end
figure
plot(mass_a,h_L_a)
grid on
xlabel('Mass Flowrate (kg/s)')
ylabel('Major Head Loss (m)')
title('Colebrook Formula - 10^-^1 Relative Error')
%----
%Solution (b)
k=1;
for m=15:1:30 %kg/s-mass flowrate
 V=m/(density*area); %m/s
 Re=density*V*D/viscosity; %>4000 --> Turbulent Flow
 rel=roughness*10^-3/D; %relative roughness
 a=-1.8*log10((rel/3.7)^1.11+6.9/Re);
 f=(1/a)^2;
 h_L_b(k)=f*(L/D)*(V^2/(2*9.81)); %m
 mass_b(k)=m; %kg/s
 k=k+1;
end

figure
plot(mass_b,h_L_b)
grid on
xlabel('Mass Flowrate (kg/s)')
ylabel('Major Head Loss (m)')
title('Haaland Formula')
%----
%Solution (c)
diff=100*(h_L_a-h_L_b)./h_L_a;
figure
plot(mass_b,diff)
grid on
xlabel('Mass Flowrate (kg/s)')
ylabel('Percent Error (%)')
title('Percent Error of Haaland Formula')
toc
