% ME 303 
% Lagrange Polynomial vs. Newton Polynomial

clear all
close all
clc

X=0:1:4;
f=@(x) cos(x);
Y=f(X);
xc=-5+min(X):0.001:max(X)+5;

%%% ---- LAGRANGE POLYNOMIAL ---- %%%

n=length(X);
L=zeros(n,n);

for k=1:n
    V=1;
    for j=1:n
        if k~=j
            V=conv(V,poly(X(j)))/(X(k)-X(j));
        end
    end
    L(k,:)=V;
end

C_l=Y*L;

yc_lag=zeros(1,length(xc));
for i=1:length(C_l)
    yc_lag=yc_lag+C_l(i)*xc.^(n-i);
end

%%% ---- NEWTON POLYNOMIAL ---- %%%

D(:,1)=Y';
for j=2:n
    for k=j:n
        D(k,j)=(D(k,j-1) - D(k-1,j-1))/...
            (X(k)-X(k-j+1));
    end
end

C_n=D(n,n);
for k=(n-1):-1:1
    C_n=conv(C_n,[1,-X(k)]);
    C_n(end)=C_n(end)+D(k,k);
end

yc_new=zeros(1,length(xc));
for i=1:length(C_n)
    yc_new=yc_new+C_n(i)*xc.^(n-i);
end

plot(xc,f(xc),'r')
hold on
plot(xc,yc_lag,'g--','Linewidth',2)
hold on
plot(xc,yc_new,'b')
hold on
plot(X,Y,'r*')
legend('f(x)','p_5_l_a_g_r_a_n_g_e(x)','p_5_n_e_w_t_o_n(x)','Data Points')
xlabel('x')
ylabel('y')
grid on
