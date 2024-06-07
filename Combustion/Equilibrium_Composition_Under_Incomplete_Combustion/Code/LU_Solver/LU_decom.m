function [n_new,iter,delta_k] = LU_decom(n_guess,Kp_Tcc,lamda,a,b)
addpath('./LU_Solver/');

% Set convergence tolerance
tol = 1e-7;
max_iter = 100;
iter = 0;

while true
    % Increment iteration count
    iter = iter + 1;
    
    % Evaluate the functions at the current guess
    [f1_s,f2_s,f3_s,f4_s,f5_s] = fk(n_guess(1),n_guess(2),n_guess(3),n_guess(4),n_guess(5),Kp_Tcc,lamda,a,b);
    f_s = [f1_s; f2_s; f3_s; f4_s; f5_s];
    
    % Compute the matrix at the current guess
    [J] = Jfk(n_guess(1),n_guess(2),n_guess(3),n_guess(4));
    
    % Compute the right-hand side vector
    b_s = -f_s;
    
    % Perform LU decomposition
    [L, U, P] = lu(J);
    
    % Solve for delta_k using forward and backward substitution
    y = L \ (P * b_s);
    delta_k = U \ y;
    
    % Update the guess values
    n_new = n_guess + delta_k;
    
    % Check for convergence
    if norm(delta_k) < tol || norm(f_s) < tol || iter >= max_iter
        break;
    end
    
    % Update the guess for the next iteration
    n_guess = n_new;
end
end