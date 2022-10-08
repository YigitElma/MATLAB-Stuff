%ME303 FINAL EXAM  04.01.2021
%QUESTION 4
%SOLVING LINEER EQUATIONS

clear all
close all
clc

%these are given in the question
A = [ 4 -1 1 1; 4 -8 1 5; 5 2 -6 8; -2 1 5 -4 ];
B = [ 10; -6; 24; 3 ];
delta = 10^(-6);  %this is my choice for accuracy

%I used this built in function to re-order matrix A so that it become
%diagonally strict. diagonals are the biggest number in the row. I could
%write a function for that but this way is more practical
[P R C] = equilibrate(A);
Per = P;

%to converge a solution I have to order A and B
a_ordered = Per*A;
b_ordered = Per*B;

%I made initial guess X_guess to find  solution, it's arbitrary
X_guess = [0;0;0;0];
N = length(B);
X = zeros(1,N);

for i=1:1000  %this part is classic Gauss-Seidel method
    for j=1:N
        if j==1
            X(1) = ( b_ordered(1)-a_ordered(1,2:N)*X_guess(2:N) )/a_ordered(1,1);
        elseif j==N
            X(N) = ( b_ordered(N)-a_ordered(N,1:N-1)*X(1:N-1)' )/a_ordered(N,N);
        else
            X(j) = ( b_ordered(j)-a_ordered(j,1:j-1)*X(1:j-1)'...
                -a_ordered(j,j+1:N)*X_guess(j+1:N) )/a_ordered(j,j);
        end
    end
    X_guess = X';  %since X is row matrix I have to transpose it
    B_new = A*X_guess;  %after each itteration I get new B_new
    d = max( abs(B_new-B) ); %difference between actual and new ones

    if d < delta %if difference is smaller than my delta break the loop
        break
    end
end

%these are the result I found
x = X_guess(1)
y = X_guess(2)
z = X_guess(3)
t = X_guess(4)


