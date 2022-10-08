clear all; 
clc;
close all;

format long

a=0.1;
b=0.9;
delta=10^(-10);
epsilon=delta;
max=100;

r=1;
V=4/3*pi*r^3;
V1=@(x) pi*(r-x).^2.*(2*r+x)/3;
V2=@(x) V-V1(x);
f=@(x) V2(x)-3*V1(x);

x=0:0.1:r;
plot(x,f(x))
grid on

a_k=[];
b_k=[];
c_k=[];
yc_k=[];

for i=1:max
    ya=f(a);
    yb=f(b);
    c=b-yb*(b-a)/(yb-ya);
    dx=min(abs(c-a),abs(c-b));
    yc=f(c);
    

    a_k=[a_k,a];
    b_k=[b_k,b];
    c_k=[c_k,c];
    yc_k=[yc_k,yc];
    
    if (dx<delta) && (abs(yc)<epsilon)
        break
    end
    
    if yb*yc>0
        b=c;       
    else
        a=c;     
    end      
end
dx
root=c
Volume1=V1(root)
Volume2=V2(root)
ratio=Volume2/Volume1


k=0:1:i-1;
k=k';
a_k=a_k';
b_k=b_k';
c_k=c_k';
yc_k=yc_k';
Data=[k,a_k,b_k,c_k,yc_k];
VarNames={'k','a_k','b_k','c_k','yc_k'};
T=table(Data(:,1),Data(:,2),Data(:,3),...
Data(:,4),Data(:,5),'variablenames',VarNames)



    



