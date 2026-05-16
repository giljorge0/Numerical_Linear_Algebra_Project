function X0 = ns_initial_guess(A)
% NS_INITIAL_GUESS  Returns the initial approximation X0 = A.' / ||A||_F^2
%
% This choice guarantees rho(I - A*X0) < 1, ensuring convergence of the
% Newton-Schulz iteration.  See proof in the project report (Question 1b).
%
% Input:
%   A  - real n x n nonsingular matrix
%
% Output:
%   X0 - initial approximation to A^{-1}

    X0 = A.' / norm(A, 'fro')^2;
end
