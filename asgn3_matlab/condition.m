function c = condition(A)
    normA = norm(A, 2);
    A_I = inverse(A);
    normA_I = norm(A_I,2);
    c = eps * normA * normA_I;
end

