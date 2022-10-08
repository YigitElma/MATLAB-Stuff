clear all
clc
close all

f=@(x) 2*sin(pi*x/6);

X=[0 1 3];
X_find=[2 2.4];

Y=f(X);

w=length(X);
n=w-1;
L=zeros(w,w);
for k=1:n+1
    V=1;
    for j=1:n+1
        if k~=j
            V=conv(V,poly(X(j)))/(X(k)-X(j));
        end
    end
    L(k,:)=V;
end
Y
L
C=Y*L
p_n=@(x) C(1)*x.^2+C(2)*x+C(3);
F=zeros(1,length(X_find));
P_n=zeros(1,length(X_find));
E=zeros(1,length(X_find));

for i=1:length(X_find)
    F(i)=f(X_find(i));
    P_n(i)=p_n(X_find(i));
    E(i)=abs(F(i)-P_n(i));
end
F
P_n
E
% x=0:0.1:6;
% plot(x,f(x),'r')
% hold on
% plot(x,p_n(x),'b')
% hold on
% grid on
% plot(X,Y,'r*')

