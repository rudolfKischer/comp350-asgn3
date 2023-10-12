function q3()
    % Q3

    %  use Lupp.m and LU factorization to calculate the derminant of a 10 x
    %  10 hilbert matrix

    % recall that the determinant of a triagular matrix is the product of
    % the diagonal entries
    % Also recall that the determinant of a product of matrices is the
    % product of the determinants of the matrices
    % |AB| = |A||B|

    % We can take advantage of the fact that we can easily factorize A

    % A = LU, so |A| = |LU| = |L||U|

    A = hilb(10);

    [L,U,P] = lupp(A);

    det_u = 1;
    det_l = 1;
    for i = 1:10
        det_u = det_u * U(i,i);
        det_l = det_l * L(i,i);
    end

    det_a = det_u * det_l;


    disp('The determinant of the 10 x 10 hilbert matrix is:');
    fprintf('|A| = |LU| = %d\n', det_a);
end

function [L,U,P] = lupp(A)
    % lupp: LU factorization with partial pivoting
    % 
    % input:  A  
    % output: L, U and P such that PA = LU
    %
    n = size(A,1);
    P = eye(n);
    
    for k = 1:n-1
       [maxval, maxindex] = max(abs(A(k:n,k)));
       q = maxindex + k - 1;
       if maxval == 0, error('A is singular'), end
       if q ~= k
           A([k,q],:) = A([q,k],:); 
           P([k,q],:) = P([q,k],:);
       end
       i = k+1:n;
       A(i,k) = A(i,k)/A(k,k);
       A(i,i) = A(i,i) - A(i,k)*A(k,i); 
    end
    
    L = tril(A,-1) + eye(n);
    U = triu(A);
end

