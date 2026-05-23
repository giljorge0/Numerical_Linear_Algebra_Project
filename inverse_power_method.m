function [lambda_min, x, iter, lambda_hist] = inverse_power_method(A, x0, tol, maxIter)
% INVERSE_POWER_METHOD  Estimates the smallest eigenvalue (in modulus) of a
%                       symmetric matrix A by applying the power iteration
%                       to A^{-1}.
%
%   A^{-1} is NEVER formed explicitly.  Instead, the product A^{-1}*x is
%   evaluated by solving the linear system  A*z = x  at each step, using a
%   pre-computed LU factorisation  P*A = L*U.
%
% Inputs:
%   A       - n x n real symmetric nonsingular matrix
%   x0      - initial vector (||x0||_inf = 1 after normalisation)
%   tol     - relative stopping tolerance on eigenvalue change
%   maxIter - maximum number of iterations
%
% Outputs:
%   lambda_min  - approximation to the smallest eigenvalue of A (in modulus)
%   x           - associated eigenvector estimate
%   iter        - number of iterations performed
%   lambda_hist - history of eigenvalue estimates of A (not of A^{-1})

    % Pre-compute LU factorisation once  (P*A = L*U)
    [L, U, P] = lu(A);

    % Normalise initial vector
    x = x0 / norm(x0, inf);

    mu_prev     = Inf;                  % eigenvalue of A^{-1}
    lambda_hist = zeros(maxIter, 1);

    for k = 1 : maxIter
        % Solve A*z = x  via forward/back substitution
        z = U \ (L \ (P * x));

        % Rayleigh quotient approximation: eigenvalue of A^{-1}
        [~, i] = max(abs(x));
        mu      = z(i) / x(i);         % eigenvalue of A^{-1}

        lambda      = 1 / mu;           % corresponding eigenvalue of A
        lambda_hist(k) = lambda;

        % Stopping criterion (relative change on mu, the iterated quantity)
        if abs(mu - mu_prev) / (abs(mu) + eps) < tol
            iter        = k;
            lambda_hist = lambda_hist(1:k);
            x           = z / mu;
            return;
        end

        x       = z / mu;              % normalise for next iteration
        mu_prev = mu;
    end

    iter        = maxIter;
    lambda_hist = lambda_hist(1:maxIter);
    warning('inverse_power_method: maximum iterations (%d) reached.', maxIter);
end
