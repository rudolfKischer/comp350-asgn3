function [A, P] = pivot(A, k, P)
    [n, ~] = size(A);
    [~, q] = max(abs(A(k:n, k)));
    q = q + k - 1; % Adjust index
    if q == k
        return;
    end
    A([k, q], k:n) = A([q, k], k:n);
    P([k, q], :) = P([q, k], :);
end
