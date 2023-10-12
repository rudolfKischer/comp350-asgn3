function A_I = inverse(A)
    n = size(A,1);
    I = eye(n);
    A_I = ggepp(A,I);
end

