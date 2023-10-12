function q2_b()

    % Solve for X_c 
    % AX = B
    % B = A * X_t
    % x_t = rand(10,4)
    % A = hilbert

    iters = 10

    iterations = 1:iters;
    epsilon = ones(1, iters) * eps;
    distances = zeros(1, iters);
    condition_numbers = zeros(1, iters);
    relative_residuals = zeros(1, iters);



    for i = 1:iters
        A = hilb(10);
        X_t = rand(10,4);
        B = A * X_t;

        X_c = ggepp(A, B);

        % Part 1
        % norm(X_c(:j))-X_t(:j))) / norm(X_t(:j))

        % ratio
        distance = norm(X_c(:,1:4) - X_t(:,1:4), 2) / norm(X_t(:,1:4), 2);
        distances(i) = distance;

        % Part 2
        % machine_epsilon norm(A) norm(A^-1)
        conddition_number = eps * cond(A, 2);
        condition_numbers(i) = conddition_number;


        % Part 3
        % compute the relative residual

        residual = B(:,1:4) - (A * X_c(:,1:4));
        relative_residual = norm(residual, 2) / (norm(A, 2) * norm(X_c(:,1:4), 2));
        relative_residuals(i) = relative_residual;
    end

    % table
    T = table(iterations', epsilon', distances', condition_numbers', relative_residuals', 'VariableNames', {'Iteration', 'Epsilon', 'Distance', 'Condition_Number', 'Relative_Residual'});

    disp(T)

    
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


