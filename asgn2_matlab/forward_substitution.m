function x = forward_substitution(A, b)
    n = size(A, 1);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1)) / A(i, i);
    end
end

