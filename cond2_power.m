function [c2, lam_max, lam_min, iter_max, iter_min] = cond2_power(A, x0, tol, maxIter)
% COND2_POWER  Estimates the 2-norm condition number of a real symmetric
%              matrix A using the power method and the inverse power method.
%
%   For a symmetric matrix:
%       cond_2(A) = ||A||_2 * ||A^{-1}||_2 = |lambda|_max / |lambda|_min
%
% The dominant eigenvalue is found via POWER_METHOD.
% The smallest eigenvalue (in modulus) is found via INVERSE_POWER_METHOD,
% which solves linear systems with the LU factorisation of A (never forms A^{-1}).
%
% Inputs:
%   A       - n x n real symmetric nonsingular matrix
%   x0      - initial vector for both methods (length n)
%   tol     - convergence tolerance (relative) for both methods
%   maxIter - maximum iterations for both methods
%
% Outputs:
%   c2       - estimated condition number cond_2(A)
%   lam_max  - estimated dominant eigenvalue
%   lam_min  - estimated smallest eigenvalue (in modulus)
%   iter_max - iterations used by power method
%   iter_min - iterations used by inverse power method

    [lam_max, ~, iter_max] = power_method(A, x0, tol, maxIter);
    [lam_min, ~, iter_min] = inverse_power_method(A, x0, tol, maxIter);

    c2 = abs(lam_max) / abs(lam_min);
end
