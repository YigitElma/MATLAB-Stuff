clc
clear all
close all

X_0=1.50; % x_0 value
X_N=2.70; % x_N value

N=1:3; % The degree of the Lagrange polynomials
f=@(x) x.^x; % The unknown function f(x) that is assumed to be known 
x_range=1.000:0.001:3.000; % The range of the x-axis

h=zeros(1,length(N)); % The increment size of the x values
C=zeros(length(N),length(N)+1); % The coefficient matrix
P=zeros(length(N),length(x_range)); % The polynom matrix

for i=1:length(N)
    h(i)=(X_N-X_0)/N(i);
    X=[1.50];
    Y=[f(1.50)];
    for j=1:N(i)   %construction of X and Y for each N value
        X=[X 1.5+h(i)*j];
        Y=f(X);
    end
    [P_N, c]=lagrange(X,Y,N(i),x_range);
    for k=1:i+1 %puts c vectors into rows of C matris
        C(i,k)=c(k);
    end  
    P(i,:)=P_N; %puts P_N vector into rows of P matris
end
P_1=P(1,:); %defines P_1 as first row of P matris
P_2=P(2,:); %I could define them in previous loop but calling the P_1,
P_3=P(3,:); %P_2 etc would require if and elseif's.
C
plot(x_range,P_1)
hold on
plot(x_range,P_2)
hold on
plot(x_range,P_3)
hold on
plot(x_range,f(x_range),'linewidth',2)
grid on
legend('P_1','P_2','P_3','f')

function [P_N, c]=lagrange(X,Y,N,x_range)
w=length(X);
L=zeros(w,w);

for k=1:w   %loop has to iterate number of data points times
    V=1;
    for j=1:w   %finding Lagrange functions
        if k~=j
            V=conv(V,poly(X(j)))/(X(k)-X(j));
        end
    end
    L(k,:)=V;
end
c=Y*L;    %definiton of coefficient vector
P_N=zeros(1,length(x_range)); %each collumn is for different collumn of x_range
for i=1:N+1  %for any value of N find the values of polynomial of degree N
    P_N=P_N+c(i)*x_range.^(N+1-i);  %with coefficients c(i).there is N+1 because in a Nth degree polynomial there is N+1 terms(x^0)
end           %At each itteration add one degree smaller term of polynomial

%P_N=polyval(c,x_range) de diyebilirdik ama ödev bunu istiyordu ek puan
%için

end
