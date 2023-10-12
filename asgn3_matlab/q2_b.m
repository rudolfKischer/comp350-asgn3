function q2_b()
    % Solve for X_c 
    % AX = B
    % B = A * X_t
    % x_t = rand(10,4)
    % A = hilbert

    iters = 10;
    iterations = (1:iters)';
    epsilon = ones(iters, 1) * eps;
    distances = zeros(iters, 4);
    condition_numbers = zeros(iters, 1);
    relative_residuals = zeros(iters, 4); 

    for i = 1:iters
        A = hilb(10);
        X_t = rand(10,4);
        B = A * X_t;
        
        % I'm using MATLAB's mldivide function here as a placeholder.
        % Replace it with your 'ggepp' function if you have it.
        X_c = A \ B;

        % Part 1
        % distances: computed per column j
        for j = 1:4
            distances(i, j) = norm(X_c(:,j) - X_t(:,j), 2) / norm(X_t(:,j), 2);
        end

        % Part 2
        condition_number = eps * cond(A, 2);
        condition_numbers(i) = condition_number;

        % Part 3
        % compute the relative residual per column j
        for j = 1:4
            residual = B(:,j) - A * X_c(:,j);
            relative_residuals(i, j) = norm(residual, 2) / (norm(A, 2) * norm(X_c(:,j), 2));
        end
    end

    % Separate tables: one for epsilon and distances, another for the rest.
    T1 = table(iterations, condition_numbers, distances(:,1), distances(:,2), distances(:,3), distances(:,4), ...
              'VariableNames', {'i', 'cond', 'error(1)', 'error(2)', 'error(3)', 'error(4)'});
    
    T2 = table(iterations, epsilon, ...
               relative_residuals(:,1), relative_residuals(:,2), relative_residuals(:,3), relative_residuals(:,4), ...
              'VariableNames', {'i', 'epsilon', ...
              'RelRes(1)', 'RelRes(2)', 'RelRes(3)', 'RelRes(4)'});
    
    disp("10 Iterations: Condition Numbers , Euclidiean Error Distances j=1:4");
    disp(T1)
    disp("10 Iterations: Machine Epsilon , Relative Residual j=1:4");
    disp(T2)





    
    % Part 3 i)
    % As seen in class the condition number given by the 2-norm of A X A^-1
    % can be used to estimate the relative error in the solution of a linear
    % system. The condition number multiplied by the machine epsilon 
    % is greater than the distance between the computed solution and the 
    % true solution. This is seen in the table above. Each distance is less than 
    % the condition number multiplied by the machine epsilon.

    % Part 3 ii)
    % we saw in class that 
    % norm(r) < (machine epsilon) * norm(A) * norm(X_c)
    %
    %r = B - AX_c
    %
    % Note this can be re arranged to our equation
    % norm( B - AX_c, 2 ) / (norm(A, 2) * norm(x_c, 2)) < (machine epsilon)
    %
    % We can see from our results that for all values of relative residuals
    % they are indeed smaller the machine epsilon
    %



    % Part 4

    % use the code to compute the inverse of a n x n non singular matrix A

    % note that A inverse is deifned as A x A_inv = I
    % where I is the identity matrix 
    % this means to find A_inv we need to solve for this equation
    %
    
    A = hilb(10);
    format;
    disp("A =");
    disp(A);

    n = size(A,1);
    I = eye(n);
    A_I = ggepp(A,I);
    
    disp("A_I =");
    disp(A_I);

    % Note the computational cost is the same as for GEPP which is 
    % 2/3n^3 + 1/2n^2


end



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

 function x = backward_substitution(A, b)
    n = size(A, 1);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:end) * x(i+1:end)) / A(i, i);
    end
 end

 function x = forward_substitution(A, b)
    n = size(A, 1);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1)) / A(i, i);
    end
end


