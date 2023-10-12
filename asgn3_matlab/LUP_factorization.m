function [L, U, P] = LUP_factorization(A)
    n = size(A, 1);
    L = eye(n);
    P = eye(n);
    U = A;
    for k = 1:n-1
        [U, P] = pivot(U, k, P);
        [U, L] = elimination(U, k, L);
    end
end
