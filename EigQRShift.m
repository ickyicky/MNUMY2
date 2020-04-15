function [eigvalues,i,success] = EigQRShift(A,tol,imax)
% EigQRShift znajduje wartości własne (D) algorytmem wykorzystującym
%   rozkład QR do sprowadzenia przez podobieństwo macierzy A
%   do macierzy diagonalnej, na której przekątnej znajdują
%   się wartości własne macierzy A. Funkcja zwraca te wartości
%   oraz ilość iteracji, która była potrzebana do
%   znalezienia wyniku. Na wejściu funckja przyjmuje argumenty:
%   A - macierz nieosobliwa
%   tol - tolerancja elementów zerowanych (górna granica), domyślnie
%   imax - maksymalna ilość iteracji
%   Algorytm z przesunięciem. W kolejnych iteracjach k
%   opuszczamy ostatni wierz i kolumnę macierzy Ak (deflacja).

n = size(A,1);
eigvalues = zeros(n,1);
i = 1;
InitialSubmatrix = A;
success = true;

for k = n : -1 : 2
    Dk = InitialSubmatrix;
    while max(abs(Dk(k,1:k-1))) > tol & success
        Dd = Dk(k-1:k, k-1:k);
        
        b = -(Dd(1,1) + Dd(2,2));
        c = Dd(2,2) * Dd(1,1) - Dd(2,1) * Dd(1,2);
        [ev1, ev2] = quadpolynroots(1,b,c);
        
        if abs(ev1 - Dd(2,2)) < abs(ev2 - Dd(2,2))
            shift = ev1;
        else
            shift = ev2;
        end
        
        Dk = Dk - eye(k) * shift;
        [Q,R] = FactorizeQR(Dk);
        Dk = R * Q + eye(k) * shift;
        i = i + 1;
        
        if i > imax
            disp('przekroczono dopuszczoną ilość maksymalnych iteracji');
            success = false;
        end
    end
    eigvalues(k) = Dk(k,k);
    if k > 2
        InitialSubmatrix = Dk(1:k-1, 1:k-1);
    else
        eigvalues(1) = Dk(1,1);
    end
end
end