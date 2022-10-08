function [X_gseidel,k_gseidel] = gseidel(A,B,P)

%AX=B lineer sistemini iteratif çözme
%temel mantık Jacobi ile aynı sadece Jacobi'de her iterasyon bulduğun
%değerleri sonraki iterasyonda kullanıyorsun. Gauss Seidel'da ise bir
%iterasyonda bulduğun değeri o iterasyon içindeki daha büyük j değerleri
%için kullanıyorsun. Böylece sonuca daha hızlı converge ediyor.

delta = 10^-9;
eps = 10^-9;
maxI = 1000;
N = length(B);

for i=1:maxI
    for j=1:N
        if j==1
            X(1) = (B(1)-A(1,2:N)*P(2:N))/A(1,1);  %j=1 değeri için önceden bulunmuş çözüm yok
        elseif j==N
            X(N) = (B(N)-A(N,1:N-1)*(X(1:N-1))')/A(N,N); %j=N için bütün çözümler önceden bulunmuş
        else
            X(j)=(B(j)-A(j,1:j-1)*(X(1:j-1))'... %aradakiler için asıl formül geçerli, karışık olan
                - A(j,j+1:N)*P(j+1:N)) / A(j,j);
        end
    end

    err = abs(norm(X'-P));
    relerr = err / (norm(X')+eps);
    P = X';

    if (err<delta) && (relerr<delta)
        break
    end
end

k_gseidel = i;
X_gseidel = X';

end

