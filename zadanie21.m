sizes = [5, 10, 20];
n = 30;

tol = 10^-10;
imax = 10^5;

result = zeros(3);

for i = 1 : size(sizes,2)
    err_b = zeros(n,1);
    err_bs = zeros(n,1);
    err_as = zeros(n,1);
    
    i_b = zeros(n,1);
    i_bs = zeros(n,1);
    i_as = zeros(n,1);
    
    s_b = zeros(n,1);
    s_bs = zeros(n,1);
    s_as = zeros(n,1);
    
    for j = 1 : n
        A = rand(sizes(i));
        B = A + A';
        
        meigb = sort(eig(B));
        meiga = sort(eig(A));
        [eigb, i_b(j), s_b(j)] = EigQRNoShift(B, tol, imax);
        [eigbs, i_bs(j), s_bs(j)] = EigQRShift(B, tol, imax);
        [eigas, i_as(j), s_as(j)] = EigQRShift(A, tol, imax);
        err_b(j) = norm(meigb - sort(eigb));
        err_bs(j) = norm(meigb - sort(eigbs));
        err_as(j) = norm(meiga - sort(eigas));
    end
    
    result(:,1) = [mean(i_b), mean(i_bs), mean(i_as)]';
    result(:,2) = abs([sum(s_b - 1), sum(s_bs - 1), sum(s_as - 1)]');
    result(:,3) = [max(err_b), max(err_bs), max(err_as)]';
    csvwrite(sprintf('./data/21/results_%d.csv', sizes(i)), result);
end