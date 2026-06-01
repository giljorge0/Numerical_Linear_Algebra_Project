%% Projeto de Álgebra Linear Numérica – LMAC 2025/26
%  Autores: Gil Jorge 110062 | Gonçalo Girante
%  Instituto Superior Técnico – Departamento de Matemática
%
%  Este ficheiro único contém todo o código do projeto (funções e scripts),
%  estruturado em secções numeradas correspondentes às questões do enunciado.
%  Pode ser aberto como Live Script no MATLAB (File > Open > seleccionar este .m).
%
%  NOTA: Antes de executar a Questão 3, coloque dois ficheiros de imagem na
%  mesma pasta e actualize os nomes nas variáveis img_gray_file e img_color_file.

format long;

%% ================================================================
%  QUESTÃO 1 – Método de Newton-Schulz
%% ================================================================
% A sucessão X_{k+1} = 2X_k - X_k A X_k converge quadraticamente para A^{-1}
% quando rho(I - A X_0) < 1.
%
% PROVA 1a)  R_{k+1} = I - A X_{k+1}
%                     = I - A(2X_k - X_k A X_k)
%                     = I - 2 A X_k + (A X_k)^2
%                     = (I - A X_k)^2  =  R_k^2
% Portanto  ||R_k|| <= ||R_0||^{2^k} -> 0  se  rho(R_0) < 1.
%
% PROVA 1b)  Com X_0 = A'/||A||_F^2 :
%   v.p. de (I - A X_0) = 1 - sigma_i^2 / ||A||_F^2  in [0,1)
%   pois  ||A||_F^2 = sum(sigma_i^2) > sigma_i^2  para A não singular com n>=2.
%   Logo rho(I - A X_0) < 1.  QED.

%% --- 1c) Parâmetros e escolha inicial ----------------------------------
maxIter = 200;
tol     = 1e-12;

%% --- 1d/1e) Vandermonde -----------------------------------------------
fprintf('\n--- Vandermonde V([1,3,5,7,9]) ---\n');
A_vand  = vander([1, 3, 5, 7, 9]);
[X, it, errs, ress] = ns_iterate(A_vand, maxIter, tol);
print_ns_result('Vandermonde', size(A_vand,1), it, errs, ress);
fprintf('||X - inv(A)||_F = %.6e\n', norm(X - inv(A_vand),'fro'));
fprintf('Computed inverse (Newton-Schulz):\n'); disp(X);
fprintf('inv(A) via MATLAB:\n');               disp(inv(A_vand));
fprintf('Convergence orders: '); disp(comp_conv_order(ress).');

%% --- Tridiagonal An ---------------------------------------------------
fprintf('\n--- Tridiagonal An ---\n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [5, 10, 50, 100, 200]
    A = tridiag_An(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n <= 10
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Conv. orders: '); disp(comp_conv_order(ress).');
    end
end

%% --- Hilbert Hn -------------------------------------------------------
fprintf('\n--- Hilbert Hn ---\n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [4, 6, 8, 10, 12]
    A = hilb(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n <= 6
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Conv. orders: '); disp(comp_conv_order(ress).');
    end
end

%% --- Lehmer Ln --------------------------------------------------------
fprintf('\n--- Lehmer Ln ---\n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [10, 100, 200, 300, 400, 500]
    A = lehmer(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n == 10
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Conv. orders: '); disp(comp_conv_order(ress).');
    end
end

%% ================================================================
%  QUESTÃO 2 – Decomposição QR de Householder
%% ================================================================
% Algoritmo:  para k = 1,...,n
%   x = R(k:m, k)
%   v = x + sign(x_1)*||x||*e_1     (sinal evita cancelamento catastrófico)
%   v = v/||v||
%   R(k:m,k:n) -= 2*v*(v'*R(k:m,k:n))
%   Qt(k:m,:)  -= 2*v*(v'*Qt(k:m,:))   % acumula Q^T
%  No fim: Q = Qt'  =>  A = Q*R

fprintf('\n=== QUESTÃO 2 – QR de Householder ===\n');
fprintf('\n--- Lehmer L_{100xn} ---\n');
fprintf('%-5s %-20s %-20s\n','n','||QR-L||_F','||Q''Q-I||_F');
for n = 10:10:100
    ri = (1:100).'; ci = 1:n;
    A  = min(ri,ci) ./ max(ri,ci);
    [Q, R] = house_qr(A);
    fprintf('%-5d %-20.6e %-20.6e\n', n, ...
            norm(Q*R-A,'fro'), norm(Q'*Q-eye(100),'fro'));
end

fprintf('\n--- Hilbert H_{100xn} ---\n');
fprintf('%-5s %-20s %-20s\n','n','||QR-H||_F','||Q''Q-I||_F');
for n = 10:10:100
    ri = (1:100).'; ci = 1:n;
    A  = 1./(ri + ci - 1);
    [Q, R] = house_qr(A);
    fprintf('%-5d %-20.6e %-20.6e\n', n, ...
            norm(Q*R-A,'fro'), norm(Q'*Q-eye(100),'fro'));
end

%% ================================================================
%  QUESTÃO 3 – Compressão de Imagens via SVD
%% ================================================================
% Qualidade:  q = sum(sigma_i^2, i=1..p) / sum(sigma_i^2, i=1..r)  x 100%
%
% Substitua os nomes de ficheiro pelos seus ficheiros de imagem.

img_gray_file  = 'image_gray.jpg';   % <- coloque o seu ficheiro aqui
img_color_file = 'image_color.jpg';  % <- coloque o seu ficheiro aqui
pct_keep       = [1, 5, 10, 25, 50, 75];

%% 3a) Imagem a preto e branco ------------------------------------------
fprintf('\n=== QUESTÃO 3a – Imagem grayscale ===\n');
A_raw = imread(img_gray_file);
if ndims(A_raw)==3, A_raw = rgb2gray(A_raw); end
A = double(A_raw);
[m_g, n_g] = size(A);
[U,S,V] = svd(A,'econ');
sv       = diag(S);
r        = length(sv);
E_tot    = sum(sv.^2);
fprintf('Dimensão: %d x %d | Nº valores singulares: %d\n', m_g, n_g, r);
fprintf('%-10s %-8s %-12s\n','% SV kept','p','Qualidade %');
figure('Name','Q3a – Grayscale');
subplot(2,4,1); imshow(uint8(A)); title('Original');
for idx = 1:length(pct_keep)
    p   = max(1, round(pct_keep(idx)/100*r));
    Ap  = U(:,1:p)*diag(sv(1:p))*V(:,1:p)';
    q   = sum(sv(1:p).^2)/E_tot*100;
    fprintf('%-10.0f %-8d %-12.4f\n', pct_keep(idx), p, q);
    subplot(2,4,idx+1); imshow(uint8(Ap));
    title(sprintf('%d%% SV (q=%.1f%%)', pct_keep(idx), q));
end
sgtitle('Compressão SVD – Grayscale');

%% 3b) Imagem a cores ---------------------------------------------------
fprintf('\n=== QUESTÃO 3b – Imagem a cores ===\n');
A_raw = imread(img_color_file);
if ndims(A_raw)~=3, error('Esperada imagem RGB.'); end
[m_c,n_c,~] = size(A_raw);
ch   = {'R','G','B'};
fprintf('%-10s %-8s %-12s %-12s %-12s %-12s\n', ...
        '% SV kept','p','Q_R %','Q_G %','Q_B %','Q_mean %');
figure('Name','Q3b – Colour');
subplot(2,4,1); imshow(A_raw); title('Original');
% Pre-compute SVDs for each channel
Uch = cell(3,1); Svch = cell(3,1); Vch = cell(3,1); Etot_ch = zeros(3,1);
for c=1:3
    Ac = double(A_raw(:,:,c));
    [Uch{c},Sc,Vch{c}] = svd(Ac,'econ');
    Svch{c}    = diag(Sc);
    Etot_ch(c) = sum(Svch{c}.^2);
end
r_c = length(Svch{1});
for idx = 1:length(pct_keep)
    p = max(1, round(pct_keep(idx)/100*r_c));
    Acomp = zeros(m_c,n_c,3,'uint8');
    q_ch  = zeros(1,3);
    for c=1:3
        Ac_rec = Uch{c}(:,1:p)*diag(Svch{c}(1:p))*Vch{c}(:,1:p)';
        Ac_rec = min(max(Ac_rec,0),255);
        Acomp(:,:,c) = uint8(Ac_rec);
        q_ch(c) = sum(Svch{c}(1:p).^2)/Etot_ch(c)*100;
    end
    fprintf('%-10.0f %-8d %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
            pct_keep(idx), p, q_ch(1), q_ch(2), q_ch(3), mean(q_ch));
    subplot(2,4,idx+1); imshow(Acomp);
    title(sprintf('%d%% SV (Q≈%.1f%%)',pct_keep(idx),mean(q_ch)));
end
sgtitle('Compressão SVD – Cores');

%% ================================================================
%  QUESTÃO 4 – Número de condição via método das potências
%% ================================================================
% cond_2(A) = |lambda|_max / |lambda|_min   (A simétrica)
% |lambda|_max : método das potências
% |lambda|_min : método das potências inversas (LU, sem inverter A)

tol4    = 1e-10;
maxit4  = 5000;

%% 4b-i) Hilbert -------------------------------------------------------
fprintf('\n=== QUESTÃO 4b-i – Hilbert Hn ===\n');
fprintf('%-5s %-20s %-20s %-20s %-8s %-8s\n', ...
        'n','cond2 (potências)','cond2 (MATLAB)','Erro relativo','it_max','it_min');
for n = 5:14
    A  = hilb(n);
    x0 = ones(n,1);
    [c2, lmax, lmin, itM, itm] = cond2_sym(A, x0, tol4, maxit4);
    c2_ml = cond(A);
    fprintf('%-5d %-20.6e %-20.6e %-20.6e %-8d %-8d\n', ...
            n, c2, c2_ml, abs(c2-c2_ml)/c2_ml, itM, itm);
end

%% 4b-ii) Lehmer -------------------------------------------------------
fprintf('\n=== QUESTÃO 4b-ii – Lehmer Ln ===\n');
fprintf('%-6s %-20s %-20s %-20s %-8s %-8s\n', ...
        'n','cond2 (potências)','cond2 (MATLAB)','Erro relativo','it_max','it_min');
for n = [10, 100, 200, 300, 400, 500]
    A  = lehmer(n);
    x0 = ones(n,1);
    [c2, ~, ~, itM, itm] = cond2_sym(A, x0, tol4, maxit4);
    if n <= 200
        c2_ml = cond(A);
        fprintf('%-6d %-20.6e %-20.6e %-20.6e %-8d %-8d\n', ...
                n, c2, c2_ml, abs(c2-c2_ml)/c2_ml, itM, itm);
    else
        fprintf('%-6d %-20.6e %-20s %-20s %-8d %-8d\n', ...
                n, c2, '(omitido)', '(omitido)', itM, itm);
    end
end

%% ================================================================
%  FUNÇÕES LOCAIS
%  (devem ficar após o fim do corpo do script)
%% ================================================================

% ----------------------------------------------------------------
function [X, iter, err_hist, res_hist] = ns_iterate(A, maxIter, tol)
% NS_ITERATE  Newton-Schulz iteration.  X_{k+1} = 2Xk - Xk*A*Xk
% Initial guess: X0 = A' / ||A||_F^2   (garantia de convergência, ver 1b)
    n   = size(A,1);
    I   = eye(n);
    X   = A.' / norm(A,'fro')^2;       % X0 from question 1b
    err_hist = NaN(maxIter+1,1);
    res_hist = zeros(maxIter+1,1);
    res_hist(1) = norm(I - A*X,'fro');

    for k = 1:maxIter
        X_new = 2*X - X*A*X;
        step  = norm(X_new - X,'fro');
        res   = norm(I - A*X_new,'fro');
        err_hist(k+1) = step;
        res_hist(k+1) = res;
        X = X_new;
        if step < tol && res < tol
            iter     = k;
            err_hist = err_hist(1:k+1);
            res_hist = res_hist(1:k+1);
            return;
        end
    end
    iter     = maxIter;
    err_hist = err_hist(1:maxIter+1);
    res_hist = res_hist(1:maxIter+1);
    warning('ns_iterate: max iterations reached.');
end

% ----------------------------------------------------------------
function p = comp_conv_order(res_hist)
% COMP_CONV_ORDER  Ordem de convergência computacional:
%   p_k = log(||Rk+1||/||Rk||) / log(||Rk||/||Rk-1||)
    res_hist = res_hist(res_hist > 0);   % remove zeros/NaN
    m = length(res_hist);
    if m < 3, p = NaN; return; end
    p = zeros(m-2,1);
    for k = 2:m-1
        num = log(res_hist(k+1)/res_hist(k));
        den = log(res_hist(k)  /res_hist(k-1));
        if abs(den) < eps, p(k-1) = NaN; else, p(k-1) = num/den; end
    end
end

% ----------------------------------------------------------------
function A = tridiag_An(n)
% TRIDIAG_AN  Matriz tridiagonal n x n do enunciado (d=10 excepto A(1,1)=9)
    d    = 10*ones(n,1); d(1) = 9;
    od   = 3*ones(n-1,1);
    A    = diag(d) + diag(od,1) + diag(od,-1);
end

% ----------------------------------------------------------------
function L = lehmer(n)
% LEHMER  Matriz de Lehmer n x n:  L(i,j) = min(i,j)/max(i,j)
    r = (1:n).'; c = 1:n;
    L = min(r,c) ./ max(r,c);
end

% ----------------------------------------------------------------
function [Q, R] = house_qr(A)
% HOUSE_QR  Decomposição QR por reflexões de Householder.
%   A = Q*R,  Q ortogonal m x m,  R triangular superior m x n.
    [m, n] = size(A);
    R  = A;
    Qt = eye(m);      % acumula Q^T

    for k = 1:n
        x   = R(k:m, k);
        s   = sign(x(1)); if s==0, s=1; end
        e1  = [1; zeros(m-k,1)];
        v   = x + s*norm(x)*e1;
        nv  = norm(v);
        if nv < eps*norm(x), continue; end
        v   = v/nv;

        R(k:m,  k:n) = R(k:m,  k:n)  - 2*v*(v.'*R(k:m,  k:n));
        Qt(k:m, :)   = Qt(k:m, :)    - 2*v*(v.'*Qt(k:m, :));
    end
    % zero sub-diagonal noise
    for k=1:n, R(k+1:m,k)=0; end
    Q = Qt.';
end

% ----------------------------------------------------------------
function [c2, lam_max, lam_min, it_max, it_min] = cond2_sym(A, x0, tol, maxit)
% COND2_SYM  cond_2(A) = |lambda_max| / |lambda_min|  (A simétrica)
%   Usa método das potências (max) e das potências inversas com LU (min).
    [lam_max, ~, it_max] = pow_max(A, x0, tol, maxit);
    [lam_min, ~, it_min] = pow_min(A, x0, tol, maxit);
    c2 = abs(lam_max)/abs(lam_min);
end

% ----------------------------------------------------------------
function [lam, x, iter] = pow_max(A, x0, tol, maxit)
% POW_MAX  Método das potências: valor próprio dominante de A.
    x        = x0 / norm(x0,inf);
    lam_prev = Inf;
    for k = 1:maxit
        z         = A*x;
        [~, i]    = max(abs(x));
        lam       = z(i)/x(i);
        x         = z/lam;
        if abs(lam-lam_prev)/(abs(lam)+eps) < tol
            iter = k; return;
        end
        lam_prev = lam;
    end
    iter = maxit;
    warning('pow_max: max iterations reached.');
end

% ----------------------------------------------------------------
function [lam_min, x, iter] = pow_min(A, x0, tol, maxit)
% POW_MIN  Método das potências inversas: menor v.p. (em módulo) de A.
%   Resolve A*z = x via factorização LU (sem calcular A^{-1}).
    [L, U, P] = lu(A);
    x         = x0 / norm(x0,inf);
    mu_prev   = Inf;
    for k = 1:maxit
        z         = U \ (L \ (P*x));
        [~, i]    = max(abs(x));
        mu        = z(i)/x(i);          % v.p. de A^{-1}
        x         = z/mu;
        if abs(mu-mu_prev)/(abs(mu)+eps) < tol
            lam_min = 1/mu; iter = k; return;
        end
        mu_prev = mu;
    end
    lam_min = 1/mu;
    iter    = maxit;
    warning('pow_min: max iterations reached.');
end

% ----------------------------------------------------------------
function print_ns_result(name, n, iter, errs, ress)
% PRINT_NS_RESULT  Imprime sumário de uma corrida do Newton-Schulz.
    fprintf('Matriz: %s  |  n=%d  |  iter=%d\n', name, n, iter);
    fprintf('  ||Xk+1-Xk||_F (final) = %.6e\n', errs(end));
    fprintf('  ||Rk||_F      (final) = %.6e\n', ress(end));
end
