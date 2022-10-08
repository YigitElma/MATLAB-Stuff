% -------------------------------------------------------------------------
% 
% Yiğit Günsür ELMACIOĞLU  -  2017405120
% 
% -------------------------------------------------------------------------

clc
clear
close all

k    = 1.667 ;
c_p  = 520.3 ;
R    = 208.1 ;
D    = 0.08 ;
V_1  = 70 ;
T_1  = 520 ;
P_1  = 350e3 ;
f    = 0.005 ;
T_exit_max = 520 ;
T_exit_min = 400 ;
d_T = 10 ;

c_in = sqrt(k*R*T_1) ;
M_in = V_1/c_in 

M = M_in ;
l_star = D/f * ( (1-M^2)/(k*M^2) + (k+1)/(2*k)*log((k+1)/2*M^2 / (1+(k-1)/2*M^2)) ) 
T_0 = T_1 * ( 1 + (k-1)/2*M^2 ) 
p_star = P_1 / ( 1/M_in * sqrt( (k+1)/2/(1+(k-1)/2*M_in^2) ) ) 
 
i = 1 ;
for T = T_exit_min:10:T_exit_max    
    M_exit = sqrt((T_0/T - 1)*2/(k-1)) ;
    P_exit = p_star * ( 1/M_exit * sqrt( (k+1)/2/(1+(k-1)/2*M_exit^2) ) ) ;
    d_S(i) = c_p*log(T/T_1) - R*log(P_exit/P_1) ;
    i = i + 1 ;
end
x = linspace(T_exit_min,T_exit_max,i-1);

% % FANNO LINE with ITERATION -----------------------------------------------
% i = 1 ;
% for T = T_exit_min:1:T_exit_max    
%     M_exit = sqrt((T_0/T - 1)*2/(k-1)) ;
%     P_exit = p_star * ( 1/M_exit * sqrt( (k+1)/2/(1+(k-1)/2*M_exit^2) ) ) ;
%     d_S_cont(i) = c_p*log(T/T_1) - R*log(P_exit/P_1) ;
%     i = i + 1 ;
% end
% x_cont = linspace(T_exit_min,T_exit_max,i-1);
% % -------------------------------------------------------------------------

% FANNO LINE with ANALYTICAL FORMULA --------------------------------------
Fanno = @(x) (c_p-R)*log(x)+R/2*log(T_0-x)+(-(c_p-R)*log(T_1)+R/2*log(2*c_p)-R*log(V_1));
x_Fanno = T_exit_min:1:T_exit_max ;
% -------------------------------------------------------------------------

plot(d_S,x,'b.','MarkerSize',25)
hold on
% plot(d_S_cont, x_cont)
% hold on
plot(Fanno(x_Fanno), x_Fanno,'r')
for i=1:13
    data =  int2str(d_S(i)) ;
    text(d_S(i),x(i),data,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',13)
end
grid on
xlabel('\bf Entropy (J/kg.K)')
ylabel('\bf Temperature (Kelvin)')
legend('Entropy Change')
title('T-s Diagram')







