function [lambda, x, iter, lambda_hist] = power_method(A, x0, tol, maxIter)
% POWER_METHOD  Estimates the dominant eigenvalue of a symmetric matrix A
%               using the power (von Mises) iteration.
%
%   Algorithm (as stated in Question 4):
%     x^(0) with ||x^(0)||_inf = 1
%     z^(k) = A * x^(k-1)
%     lambda^(k) = z^(k)_i / x^(k-1)_i    (i: any index where x^(k-1)_i ~= 0)
%     x^(k) = z^(k) / lambda^(k)
%     stop when |lambda^(k) - lambda^(k-1)| / |lambda^(k)| < tol
%
% Inputs:
%   A       - n x n real symmetric matrix
%   x0      - initial vector (will be normalised so ||x0||_inf = 1)
%   tol     - relative stopping tolerance on eigenvalue change
%   maxIter - maximum number of iterations
%
% Outputs:
%   lambda      - approximation to the dominant eigenvalue
%   x           - approximation to the associated eigenvector
%   iter        - number of iterations performed
%   lambda_hist - history of eigenvalue estimates

    % Normalise initial vector
    x = x0 / norm(x0, inf);

    lambda_prev  = Inf;
    lambda_hist  = zeros(maxIter, 1);

    for k = 1 : maxIter
        z = A * x;

        % Choose index i where x is largest in magnitude (most stable choice)
        [~, i] = max(abs(x));
        lambda  = z(i) / x(i);

        lambda_hist(k) = lambda;

        % Stopping criterion (relative change)
        if abs(lambda - lambda_prev) / (abs(lambda) + eps) < tol
            iter        = k;
            lambda_hist = lambda_hist(1:k);
            x           = z / lambda;          % final eigenvector estimate
            return;
        end

        x          = z / lambda;               % normalise: ||x||_inf = 1
        lambda_prev = lambda;
    end

    iter        = maxIter;
    lambda_hist = lambda_hist(1:maxIter);
    warning('power_method: maximum iterations (%d) reached.', maxIter);
end
