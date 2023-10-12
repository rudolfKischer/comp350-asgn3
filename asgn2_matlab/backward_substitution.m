function x = backward_substitution(A, b)
    n = size(A, 1);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:end) * x(i+1:end)) / A(i, i);
    end
end

