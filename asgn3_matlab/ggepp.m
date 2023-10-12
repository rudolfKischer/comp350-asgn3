
 function X = ggepp(A, B)
     [n, ~] = size(A);
     [L,U,P] =  lupp(A);
     [~, p] = size(B);
     Y = zeros(n, p);
     X = zeros(n, p);
     for i = 1:p
         y = forward_substitution(L, P * B(:, i));
         x = backward_substitution(U, y);
         Y(:, i) = y;
         X(:, i) = x;
     end
end

% function X = ggepp(A, B)
%     [n, ~] = size(A);
%     L = eye(n);
%     P = eye(n);
%     U = A;
%     [~, p] = size(B);
%     for k = 1:n-1
%         [U, P] = pivot(U, k, P);
%         [U, L] = elimination(U, k, L);
%     end
%     Y = zeros(n, p);
%     X = zeros(n, p);
%     for i = 1:p
%         y = forward_substitution(L, P * B(:, i));
%         x = backward_substitution(U, y);
%         Y(:, i) = y;
%         X(:, i) = x;
%     end
% end
