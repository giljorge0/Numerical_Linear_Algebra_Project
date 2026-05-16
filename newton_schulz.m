function [X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol)
% NEWTON_SCHULZ  Approximates A^{-1} via the Newton-Schulz iteration:
%                   X_{k+1} = 2*X_k - X_k * A * X_k
%
% Inputs:
%   A       - nonsingular square matrix (n x n)
%   X0      - initial approximation of A^{-1} (n x n)
%   maxIter - maximum number of iterations
%   tol     - convergence tolerance (used for both step and residual)
%
% Outputs:
%   X        - computed approximation of A^{-1}
%   iter     - number of iterations performed
%   err_hist - history of ||X_{k+1} - X_k||_F  (step errors)
%   res_hist - history of ||I - A*X_k||_F       (residuals)
%
% Stopping criteria:
%   (1)  ||X_{k+1} - X_k||_F < tol
%   (2)  ||I - A * X_{k+1}||_F < tol
%   Both must be satisfied, or maxIter is reached.

    n        = size(A, 1);
    I        = eye(n);
    X        = X0;
    err_hist = zeros(maxIter, 1);
    res_hist = zeros(maxIter, 1);

    % Compute residual of initial iterate
    R = I - A * X;
    res_hist(1) = norm(R, 'fro');
    err_hist(1) = NaN;          % no previous iterate for step error

    for k = 1 : maxIter
        X_new = 2*X - X*A*X;

        % ------ stopping quantities ------
        step = norm(X_new - X, 'fro');
        R    = I - A * X_new;
        res  = norm(R, 'fro');

        % store (shift by 1 because index 1 holds the initial residual)
        if k + 1 <= maxIter
            err_hist(k + 1) = step;
            res_hist(k + 1) = res;
        end

        X = X_new;

        if step < tol && res < tol
            iter = k;
            err_hist = err_hist(1 : k + 1);
            res_hist = res_hist(1 : k + 1);
            return;
        end
    end

    iter     = maxIter;
    err_hist = err_hist(1 : maxIter);
    res_hist = res_hist(1 : maxIter);

    warning('newton_schulz: maximum iterations (%d) reached.', maxIter);
end
