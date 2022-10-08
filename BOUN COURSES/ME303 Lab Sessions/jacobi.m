function [X_jacobi k_jacobi]=jacobi(A, B, P)

%AX=B şeklindeki lineer sistemin çözümünü bulmak
%P başlangıçta kullandığımız tahmin değeri
%X_jacobi asıl sonuç
%k_jacobi iterasyon sayısı
%Jacobi metodunun mantıklı sonuç vermesi için A matrisinin diagonallerinin
%diğerlerine göre daha büyük sayılar olması lazım. olmadığı durumlarda
%düzeltmek için ödevlerdeki soruya bak. Ya da "permutation" diye yazdığım fonksiyona bak


delta=10^(-9);
eps=10^(-9);
maxI=1000;
N=length(B);

for i=1:maxI
    for j=1:N   
        X(j)=(B(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
    end
    err=abs(norm(X'-P));
    relerr=err/(norm(X')+eps);
    P=X';
    
    if err<delta && relerr<delta
        break
    end
end
k_jacobi=i;
X_jacobi=X';