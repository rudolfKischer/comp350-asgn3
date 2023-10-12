function q2_b()

    % Solve for X_c 
    % AX = B
    % B = A * X_t
    % x_t = rand(10,4)
    % A = hilbert

    iters = 10

    distances = zeros(1, iters);
    condition_numbers = zeros(1, iters);
    relative_residuals = zeros(1, iters);
    iteration_number



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
    T = table(distances', condition_numbers', relative_residuals', 'VariableNames', {'Distance', 'Condition_Number', 'Relative_Residual'});

    disp(T)

    
    % Part 3 i)

    % Part 3 ii)

    % Part 4


    






end