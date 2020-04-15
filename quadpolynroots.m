function [x1,x2] = quadpolynroots(a,b,c)
% quadpolynroots zwraca pierwiastki wielomianu
%   stopnia drugiego o współczynnikach a,b,c

delta = sqrt(b^2 - 4*a*c);
x1 = (-b - delta) / (2*a);
x2 = (-b + delta) / (2*a);
end