%% script_part2.m
% Applies the Householder QR decomposition to rectangular Lehmer and
% Hilbert matrices  A ∈ R^{100 x n},  n = 10, 20, ..., 100.
%
% For each matrix and each n the script computes:
%   ||Q*R - A||_F        (reconstruction error)
%   ||Q'*Q - I||_F       (orthogonality error)
% and collects them in vectors for plotting and tabulation.

format long;

ns = 10 : 10 : 100;     % n = 10, 20, ..., 100
m  = 100;

err_QR_Lehmer  = zeros(length(ns), 1);
err_QQ_Lehmer  = zeros(length(ns), 1);
err_QR_Hilbert = zeros(length(ns), 1);
err_QQ_Hilbert = zeros(length(ns), 1);

fprintf('\n=== Householder QR: Lehmer matrices  L_{100 x n} ===\n');
fprintf('%-5s  %-20s  %-20s\n', 'n', '||QR - L||_F', '||Q''Q - I||_F');
fprintf('%s\n', repmat('-', 1, 50));

for idx = 1 : length(ns)
    n = ns(idx);

    % Build 100 x n Lehmer matrix (submatrix of the full n-square one)
    % The (i,j) entry is min(i,j)/max(i,j) for i=1..m, j=1..n
    [I, J] = meshgrid(1:n, 1:m);      % I(r,c)=c, J(r,c)=r  -> fix below
    % meshgrid: rows vary with 2nd arg, cols vary with 1st arg
    % so I has column index values, J has row index values
    row_idx = (1:m).';
    col_idx = 1:n;
    A = min(row_idx, col_idx) ./ max(row_idx, col_idx);

    [Q, R] = householder_qr(A);

    e1 = norm(Q * R - A, 'fro');
    e2 = norm(Q.' * Q - eye(m), 'fro');

    err_QR_Lehmer(idx) = e1;
    err_QQ_Lehmer(idx) = e2;

    fprintf('%-5d  %-20.6e  %-20.6e\n', n, e1, e2);
end

fprintf('\n=== Householder QR: Hilbert matrices  H_{100 x n} ===\n');
fprintf('%-5s  %-20s  %-20s\n', 'n', '||QR - H||_F', '||Q''Q - I||_F');
fprintf('%s\n', repmat('-', 1, 50));

for idx = 1 : length(ns)
    n = ns(idx);

    % Build 100 x n Hilbert submatrix:  H(i,j) = 1/(i+j-1)
    row_idx = (1:m).';
    col_idx = 1:n;
    A = 1 ./ (row_idx + col_idx - 1);

    [Q, R] = householder_qr(A);

    e1 = norm(Q * R - A, 'fro');
    e2 = norm(Q.' * Q - eye(m), 'fro');

    err_QR_Hilbert(idx) = e1;
    err_QQ_Hilbert(idx) = e2;

    fprintf('%-5d  %-20.6e  %-20.6e\n', n, e1, e2);
end

%% ---- Plots ----
figure(1);
semilogy(ns, err_QR_Lehmer,  'b-o', 'DisplayName', '||QR - L||_F');
hold on;
semilogy(ns, err_QQ_Lehmer,  'r-s', 'DisplayName', '||Q''Q - I||_F');
xlabel('n'); ylabel('Frobenius norm error'); grid on;
title('Householder QR – Lehmer matrices L_{100 \times n}');
legend('Location','northwest'); hold off;

figure(2);
semilogy(ns, err_QR_Hilbert, 'b-o', 'DisplayName', '||QR - H||_F');
hold on;
semilogy(ns, err_QQ_Hilbert, 'r-s', 'DisplayName', '||Q''Q - I||_F');
xlabel('n'); ylabel('Frobenius norm error'); grid on;
title('Householder QR – Hilbert matrices H_{100 \times n}');
legend('Location','northwest'); hold off;

%% ---- Comparison with MATLAB's built-in QR for n=100 ----
fprintf('\n=== Sanity check at n=100: comparing with MATLAB qr() ===\n');

% Lehmer
row_idx = (1:m).'; col_idx = 1:100;
A = min(row_idx, col_idx) ./ max(row_idx, col_idx);
[Q_hs, R_hs] = householder_qr(A);
[Q_ml, R_ml] = qr(A);
fprintf('Lehmer 100x100:\n');
fprintf('  Our QR:     ||QR - A||_F = %.6e,  ||Q''Q - I||_F = %.6e\n', ...
        norm(Q_hs*R_hs - A,'fro'), norm(Q_hs.'*Q_hs - eye(m),'fro'));
fprintf('  MATLAB qr:  ||QR - A||_F = %.6e,  ||Q''Q - I||_F = %.6e\n', ...
        norm(Q_ml*R_ml - A,'fro'), norm(Q_ml.'*Q_ml - eye(m),'fro'));

% Hilbert
A = 1 ./ (row_idx + col_idx - 1);
[Q_hs, R_hs] = householder_qr(A);
[Q_ml, R_ml] = qr(A);
fprintf('Hilbert 100x100:\n');
fprintf('  Our QR:     ||QR - A||_F = %.6e,  ||Q''Q - I||_F = %.6e\n', ...
        norm(Q_hs*R_hs - A,'fro'), norm(Q_hs.'*Q_hs - eye(m),'fro'));
fprintf('  MATLAB qr:  ||QR - A||_F = %.6e,  ||Q''Q - I||_F = %.6e\n', ...
        norm(Q_ml*R_ml - A,'fro'), norm(Q_ml.'*Q_ml - eye(m),'fro'));
