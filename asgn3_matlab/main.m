function main()
    A = [-3, 0, 24, 6;
        3, 9, -9, 24;
        6, 6, -6, 24;
        2, 5, 1, 29];

    disp("A = ");
    disp(A);

    B = [-15;
                39;
                30;
                31];
    
    disp("B = ");
    disp(B);

    [L, U, P] = LUP_factorization(A);

    disp("L = ");
    disp(L);
    disp("U = ");
    disp(U);
    disp("P = ");
    disp(P);

    X = ggepp(A, B);

    disp("X = ");
    disp(X);



end

    