function Per=permutation(A)
%26.11.2020
%NxN lik bir matrisin diagonaldeki girdilerini en büyük yapcak şekilde
%bir permütasyon matrisi verir.Satırlar sıralanmış olur.
%Per*A=B işlemi sonucunda B'nin diagonali her satırdaki en büyük değere 
%sahiptir

N=length(A);
order=[1:N];

for i=1:N
    check=0;  %başlangıç için düşük değer
    for j=1:N    %her satırdaki en büyük sayıyı bulma
        if abs(A(i,j))>check   
            check=abs(A(i,j));  %satırdaki en büyük değer check içinde kayıtlı
            order(i)=j;   %en büyük değerin konumu orderda kayıtlı
        end
    end
end
Per=zeros(N,N);
for i=1:N
    for j=1:N
        if order(j)==i  %en büyük sayının bulunduğu konum 1 olmalı
            Per(i,j)=1; %permütasyon matrisi özelliğinden
        end
    end
end
% Programı geliştirmek adına en büyük değeri aynı kolonlarda olan girdileri
% düzenlemek için bir şey yapılabilir.
% örnek: a=[5 9 7;6 78 11;81 9 47]; 
% ans =
%      0     0     1
%      1     1     0
%      0     0     0
% sonucunu vermemeli. 78, 9 dan daha büyük olduğu için ikinci satıra onu
% almalı gibi.
end
