function L = build_Lehmer(n)
% BUILD_LEHMER  Constructs the n x n Lehmer matrix
%
%   (L)_{ij} = min(i,j) / max(i,j)
%
% Equivalently: (L)_{ij} = i/j for j >= i,  (L)_{ij} = j/i for j < i.
% The matrix is symmetric, positive-definite, and has entries in (0,1].
%
% Input:
%   n - matrix dimension (positive integer)
%
% Output:
%   L - n x n Lehmer matrix

    [I, J] = meshgrid(1:n, 1:n);   % I(i,j) = i,  J(i,j) = j
    I = I.';  J = J.';              % fix orientation: rows=i, cols=j
    L = min(I, J) ./ max(I, J);
end
