x = [-5 : 5];
y = [-5.8643,
    -6.7445 ,
    -5.2378,
    -3.2868,
    -2.2393,
    -0.5084,
    -1.2237,
    -0.7893,
    -4.8761,
    -11.0466,
    -20.0868,
];

n_max = 6;
xval = [-5 : 0.1 : 5];

figure(1);
plot(x, y, 'o', 'DisplayName', 'Próbki');
title('Metoda najmniejszych kwardatów');
hold on;
legend('show', 'Location', 'northwest');
legend('boxoff');

errors_eukl = zeros(1, n_max + 1);
errors_cz = zeros(1, n_max + 1);

for n = 0 : n_max
    a = LeastSquares(x,y,n);
    yval = polyval(a, xval);
    %plot(xval, yval, 'DisplayName', sprintf('Stopień wielomianu: %d', n));
    errors_eukl(n + 1) = norm(polyval(a, x) - y');
    errors_cz(n + 1) = max(abs(polyval(a, x) - y'));
end
plot(xval, yval, 'DisplayName', sprintf('Stopień wielomianu: %d', n));
disp(a);
hold off;
saveas(1, './plots/22/least_squares.png');
saveas(1, './plots/22/least_squares.fig');

disp('Normy Euklidesowe dla motody z użyciem równań normalnych:');
disp(errors_eukl);
disp('Normy Czebyszewa');
disp(errors_cz);
csvwrite('./data/22/errors_eukl.csv', errors_eukl);
csvwrite('./data/22/errors_cz.csv', errors_cz);

figure(2);
plot(x, y, 'o', 'DisplayName', 'Próbki');
title('Metoda najmniejszych kwardatów');
hold on;
legend('show', 'Location', 'northwest');
legend('boxoff');

for n = 0 : n_max
    a = LeastSquaresQR(x,y,n);
    yval = polyval(a, xval);
    %plot(xval, yval, 'DisplayName', sprintf('Stopień wielomianu: %d', n));
    errors_eukl(n + 1) = norm(polyval(a, x) - y');
    errors_cz(n + 1) = max(abs(polyval(a, x) - y'));
end
plot(xval, yval, 'DisplayName', sprintf('Stopień wielomianu: %d', n));
disp(a);
hold off;
saveas(2, './plots/22/least_squares_qr.png');
saveas(2, './plots/22/least_squares_qr.fig');

disp('Normy Euklidesowe dla motody z użyciem rozkładu QR:');
disp(errors_eukl);
disp('Normy Czebyszewa');
disp(errors_cz);
csvwrite('./data/22/errors_qr_eukl.csv', errors_eukl);
csvwrite('./data/22/errors_qr_cz.csv', errors_cz);