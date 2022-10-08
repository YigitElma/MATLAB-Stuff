clc
clear all
close all

x=[7 5 8 13 17 24];
y=[2 -14 -13 27 -16 4];
n=length(x);
D=zeros(n,n);
D(:,1)=y';

for j=2:n
    for k=j:n
        D(k,j)=(D(k,j-1)-D(k-1,j-1))/(x(k)-x(k-j+1));
    end
end

C=D(n,n);
for k=(n-1):-1:1
    C=conv(C,[1,-x(k)]);
    C(end)=C(end)+D(k,k);
end

xc=min(x):0.001:max(x);
yc=zeros(1,length(xc));

for i=1:length(C)
    yc=yc+C(i)*xc.^(n-i);
end
plot(xc,yc)
grid on
hold on
xlabel('x')
ylabel('y')
for i=1:n
    plot(x(i),y(i),'r*')
    text(x(i),y(i),num2str(i),'fontsize',16,'color',...
        [0 1 0],'verticalalignment','top');
end

