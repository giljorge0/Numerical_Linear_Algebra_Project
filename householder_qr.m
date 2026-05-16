function [Q, R] = householder_qr(A)
% HOUSEHOLDER_QR  Computes the QR decomposition of A via Householder
%                 reflections.
%
%   [Q, R] = HOUSEHOLDER_QR(A)
%
%   Given A ∈ R^{m x n} with m >= n, returns orthogonal Q ∈ R^{m x m}
%   and upper-triangular R ∈ R^{m x n} such that  A = Q * R.
%
% Algorithm:
%   At step k we choose a Householder reflector H_k = I - 2*v*v' that
%   zeroes out entries k+1,...,m in column k.  The overall factorisation
%   is:
%       H_n * ... * H_1 * A = R
%   so
%       Q = H_1 * H_2 * ... * H_n     (each H_k is symmetric & orthogonal)
%
%   We accumulate Q^T = H_n * ... * H_1 applied to the identity, then
%   transpose at the end.
%
% Inputs:
%   A - real m x n matrix, m >= n
%
% Outputs:
%   Q - m x m orthogonal matrix
%   R - m x n upper-triangular matrix

    [m, n] = size(A);
    assert(m >= n, 'HOUSEHOLDER_QR: require m >= n.');

    R = A;           % will be overwritten column by column
    Q = eye(m);      % accumulates Q^T; transposed at the end

    for k = 1 : n
        % Extract the sub-column that must be zeroed (except first entry)
        x = R(k:m, k);
        len_x = m - k + 1;

        % Householder vector: v = x + sign(x(1))*||x||*e_1
        % The sign choice avoids catastrophic cancellation.
        s = sign(x(1));
        if s == 0, s = 1; end        % treat zero as positive

        e1    = zeros(len_x, 1);
        e1(1) = 1;
        v     = x + s * norm(x) * e1;

        % If x is already a multiple of e1, the reflector is trivial
        nv = norm(v);
        if nv < eps * norm(x)
            continue;
        end
        v = v / nv;                  % normalise  (||v|| = 1)

        % Apply H_k = I - 2*v*v' to the active sub-block of R:
        %   R(k:m, k:n) = R(k:m, k:n) - 2*v*(v'*R(k:m, k:n))
        R(k:m, k:n) = R(k:m, k:n) - 2 * v * (v.' * R(k:m, k:n));

        % Accumulate Q^T:  Q = H_k * Q  (applied from the left in subspace)
        Q(k:m, :) = Q(k:m, :) - 2 * v * (v.' * Q(k:m, :));
    end

    % Force strict upper-triangular structure (eliminate rounding noise below diagonal)
    for k = 1 : n
        R(k+1:m, k) = 0;
    end

    Q = Q.';         % Q^T was accumulated; transpose to get Q
end
