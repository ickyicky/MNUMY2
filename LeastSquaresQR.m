function a = LeastSquaresQR(x, y, n)
% LeastSquaresQR oblicza współczynniki wielomianu
%  stopnia n, najdoładniej przybliżające wartości
%  zadanej próbki x, y. Metoda wykrozystująca rozkład
%  QR.

A = zeros(length(x), n + 1);

for i = 1 : length(x)
    for j = 1 : n + 1
        A(i, j) = x(i) ^ (j -1);
    end
end

[Q, R] = FactorizeQR(A);
a = R \ (Q' * y);
a = fliplr(a');
end