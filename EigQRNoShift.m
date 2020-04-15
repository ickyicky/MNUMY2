function [eigvalues,i,success] = EigQRNoShift(A,tol,imax)
% EigQRNoShift znajduje wartości własne (D) algorytmem wykorzystującym
%   rozkład QR do sprowadzenia przez podobieństwo macierzy A
%   do macierzy diagonalnej, na której przekątnej znajdują
%   się wartości własne macierzy A. Funkcja zwraca te wartości
%   oraz ilość iteracji, która była potrzebana do
%   znalezienia wyniku. Na wejściu funckja przyjmuje argumenty:
%   A - macierz nieosobliwa
%   tol - tolerancja elementów zerowanych (górna granica)
%   imax - maksymalna ilość iteracji
%   Algorytm bez przesunięcia.

i = 1;
success = true;

while max(max(abs(A - diag(diag(A))))) > tol & success
    [Q1, R1] = FactorizeQR(A);
    A = R1 * Q1;
    i = i + 1;
    if i > imax
        disp('przekroczono dopuszczoną ilość maksymalnych iteracji');
        success = false;
    end
end

eigvalues = diag(A);

end