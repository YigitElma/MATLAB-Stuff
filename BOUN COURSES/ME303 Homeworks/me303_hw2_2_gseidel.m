clear all
close all
clc
%A*x=B lineer sisteminin Gauss Seidel yöntemiyle çözülmesi
A=[3 7 16 2 2;2 4 8 24 8;1 15 5 6 2;2 1 -10 7 -21;-13 1 -3 -4 -2];
B=[14;23;4;17;9];
X_star=[0;0;0;0;0];
[P,R,C]=equilibrate(A); %finds matrix P which makes absolute value of diagonal elements maximum
Per=P;  %permutation matrix is P

delta=10^-10;

N=length(A);
A_ordered=Per*A; %to converge to result we must take this product
B_ordered=Per*B; %also this

for i=1:100
    for j=1:N  %formula of Gauss Seidel method
        if j==1
            X(1) = (B_ordered(1)-A_ordered(1,2:N)*X_star(2:N))/A_ordered(1,1);
        elseif j==N
            X(N) = (B_ordered(N)-A_ordered(N,1:N-1)*(X(1:N-1))')/A_ordered(N,N);
        else
            X(j)=(B_ordered(j)-A_ordered(j,1:j-1)*(X(1:j-1))'...
                - A_ordered(j,j+1:N)*X_star(j+1:N)) / A_ordered(j,j);
        end
    end
    
    X_star = X'; %X is row matrix, we take transpose
    B_star=A_ordered*X_star; %B_star is calculated by X_star which we calculate by G.Seidel
   
    err=0; %this is just to give a value to err which is smaller than everything
    for k=1:N   %to find the max(abs(B*k-Bk)) element
        if abs(B_star(k)-B(k))>err
            err=abs(B_star(k)-B(k));
        end
    end

    if err<delta
        break
    end
end
i
X_star
Per