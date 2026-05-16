function A = build_An(n)
% BUILD_AN  Constructs the n x n tridiagonal matrix
%
%        [10  3  0  ...  0 ]          (diagonal = 10, except A(1,1)=9)
%        [ 3 10  3  ...  0 ]
%   An = [          ...    ]
%        [ 0  ... 3  10  3 ]
%        [ 0  ...  0  3 10 ]
%
% As stated in the project, the (1,1) entry is 9 and all other diagonal
% entries are 10, with off-diagonal entries equal to 3.
%
% Input:
%   n - matrix dimension (positive integer)
%
% Output:
%   A - n x n tridiagonal matrix

    d  = 10 * ones(n, 1);
    d(1) = 9;                       % top-left corner is 9, not 10
    od = 3  * ones(n-1, 1);

    A = diag(d) + diag(od, 1) + diag(od, -1);
end
