function [optimal_a, optimal_b, optimal_c, min_eigenvalue_real] = minimize_eigenvalue(A)
    % A is the fixed 15x15 matrix
    
    % Initial guesses for a, b, c
    initial_guess = [0, 0, 0]; % Modify based on your starting values

    % Optimization using fminsearch (or fmincon if constraints are needed)
    options = optimset('Display', 'iter');  % Set options for verbose output
    ub=[.75,.75,.75];
    lb=-ub;

   [optimal_params, min_eigenvalue_real] = fmincon(@(params) objective_function(A, params), ...
                                                initial_guess, [], [], [], [], lb, ub, [], options);

    
    % Extract optimized a, b, and c
    optimal_a = optimal_params(1);
    optimal_b = optimal_params(2);
    optimal_c = optimal_params(3);
end

function f = objective_function(A, params)
    a = params(1);
    b = params(2);
    c = params(3);
    
    % Call your jacobian function to get eigenvalues
    eigenvalues = Jacobian(A, a, b, c);
    
    % Get the smallest real part of the eigenvalues
    biggest_real_part = max(real(eigenvalues));
    
    % Objective function to minimize
    f = biggest_real_part;
end

