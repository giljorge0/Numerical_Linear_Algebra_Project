%% script_part1de.m
% Tests the Newton-Schulz iteration on the four matrix families required
% in Questions 1(d) and 1(e).
%
% For each family the script:
%   - builds the matrix A
%   - sets X0 = A.' / ||A||_F^2
%   - calls newton_schulz with combined stopping criteria
%   - prints a table of n, iterations, final step error, final residual
%   - for small n, compares X with inv(A)
%   - estimates computational convergence order (Question 1e)

format long;
maxIter = 200;
tol     = 1e-12;

%% =========================================================
%  1. Vandermonde matrix  V([1 3 5 7 9])
%% =========================================================
fprintf('\n========== Vandermonde V([1,3,5,7,9]) ==========\n');

pts = [1, 3, 5, 7, 9];
A   = vander(pts);       % MATLAB built-in: columns v.^(n-1), ..., v.^0
n   = size(A, 1);

X0  = ns_initial_guess(A);
[X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol);

fprintf('n = %d | Iterations = %d\n', n, iter);
fprintf('||X_{k+1}-X_k||_F (final) = %.6e\n', err_hist(end));
fprintf('||R_k||_F         (final) = %.6e\n', res_hist(end));

% Compare with MATLAB inv
A_inv_matlab = inv(A);
fprintf('||X - inv(A)||_F          = %.6e\n', norm(X - A_inv_matlab, 'fro'));

fprintf('\nComputed inverse (Newton-Schulz):\n');  disp(X);
fprintf('inv(A) via MATLAB:\n');                    disp(A_inv_matlab);

% Convergence order
p = convergence_order(res_hist);
fprintf('Computational convergence orders:\n');
disp(p.');

%% =========================================================
%  2. Tridiagonal matrix An  (several values of n)
%% =========================================================
fprintf('\n========== Tridiagonal matrix An ==========\n');
ns_tri = [5, 10, 50, 100, 200];

fprintf('%-6s  %-10s  %-18s  %-18s\n', ...
        'n', 'iter', '||X_{k+1}-X_k||_F', '||R_k||_F');
fprintf('%s\n', repmat('-', 1, 60));

for n = ns_tri
    A  = build_An(n);
    X0 = ns_initial_guess(A);
    [X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol);

    fprintf('%-6d  %-10d  %-18.6e  %-18.6e\n', ...
            n, iter, err_hist(end), res_hist(end));

    if n <= 10
        fprintf('  ||X - inv(A)||_F = %.6e\n', norm(X - inv(A), 'fro'));
        fprintf('  Convergence orders: ');
        disp(convergence_order(res_hist).');
    end
end

%% =========================================================
%  3. Hilbert matrix Hn  (several values of n)
%% =========================================================
fprintf('\n========== Hilbert matrix Hn ==========\n');
ns_hilb = [4, 6, 8, 10, 12];

fprintf('%-6s  %-10s  %-18s  %-18s\n', ...
        'n', 'iter', '||X_{k+1}-X_k||_F', '||R_k||_F');
fprintf('%s\n', repmat('-', 1, 60));

for n = ns_hilb
    A  = hilb(n);        % MATLAB built-in Hilbert matrix
    X0 = ns_initial_guess(A);
    [X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol);

    fprintf('%-6d  %-10d  %-18.6e  %-18.6e\n', ...
            n, iter, err_hist(end), res_hist(end));

    if n <= 6
        fprintf('  ||X - inv(A)||_F = %.6e\n', norm(X - inv(A), 'fro'));
        fprintf('  Convergence orders: ');
        disp(convergence_order(res_hist).');
    end
end
% Note: for large n, Hilbert matrices are extremely ill-conditioned;
% convergence may stall or the method may fail. Comment on this.

%% =========================================================
%  4. Lehmer matrix  Ln  (n = 10, 100, 200, 300, 400, 500)
%% =========================================================
fprintf('\n========== Lehmer matrix Ln ==========\n');
ns_lehmer = [10, 100, 200, 300, 400, 500];

fprintf('%-6s  %-10s  %-18s  %-18s\n', ...
        'n', 'iter', '||X_{k+1}-X_k||_F', '||R_k||_F');
fprintf('%s\n', repmat('-', 1, 60));

for n = ns_lehmer
    A  = build_Lehmer(n);
    X0 = ns_initial_guess(A);
    [X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol);

    fprintf('%-6d  %-10d  %-18.6e  %-18.6e\n', ...
            n, iter, err_hist(end), res_hist(end));

    if n == 10
        fprintf('  ||X - inv(A)||_F = %.6e\n', norm(X - inv(A), 'fro'));
        fprintf('  Convergence orders: ');
        disp(convergence_order(res_hist).');
    end
end
