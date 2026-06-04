format long;
% 1 – Método de Newton-Schulz
%% 1c Parâmetros e escolha inicial
maxIter = 200;
tol     = 1e-12;
%% 1d/1e) Vandermonde
fprintf('\n--- Vandermonde V([1,3,5,7,9]) ---\n');
A_vand  = vander([1, 3, 5, 7, 9]);
[X, it, errs, ress] = ns_iterate(A_vand, maxIter, tol);
print_ns_result('Vandermonde', size(A_vand,1), it, errs, ress);
fprintf('||X - inv(A)||_F = %.6e\n', norm(X - inv(A_vand),'fro'));
fprintf('Inversa calculada (Newton-Schulz):\n'); disp(X);
fprintf('inv(A) via MATLAB');                 disp(inv(A_vand));
fprintf('Ordens de convergência: '); disp(comp_conv_order(ress).');
%% Tridiagonal An
fprintf('\n Tridiagonal An \n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [5, 10, 50, 100, 200]
    A = tridiag_An(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n <= 10
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Ordens de conv.: '); disp(comp_conv_order(ress).');
    end
end
%% Hilbert
fprintf('\n Hilbert Hn\n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [4, 6, 8, 10, 12]
    A = hilb(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n <= 6
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Ordens de conv.: '); disp(comp_conv_order(ress).');
    end
end
%%Lehmer Ln
fprintf('\nLehmer Ln\n');
fprintf('%-6s %-8s %-18s %-18s\n','n','iter','||Xk+1-Xk||_F','||Rk||_F');
for n = [10, 100, 200, 300, 400, 500]
    A = lehmer(n);
    [X, it, errs, ress] = ns_iterate(A, maxIter, tol);
    fprintf('%-6d %-8d %-18.6e %-18.6e\n', n, it, errs(end), ress(end));
    if n == 10
        fprintf('  ||X-inv(A)||_F = %.6e\n', norm(X-inv(A),'fro'));
        fprintf('  Ordens de conv.: '); disp(comp_conv_order(ress).');
    end
end
%%  Gráficos de Convergência (Exemplo com Matriz Tridiagonal n=10)
% Configurar uma matriz bem condicionada com boa convergência
A_plot = tridiag_An(10);
[~, it_plot, errs_plot, ress_plot] = ns_iterate(A_plot, 200, 1e-12);
figure('Name', 'Convergência Newton-Schulz', 'Position', [100, 100, 900, 400]);
% Subplot 1: Erros ao longo do tempo (Escala Logarítmica)
subplot(1, 2, 1);
semilogy(0:it_plot, errs_plot, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
semilogy(0:it_plot, ress_plot, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on;
xlabel('Iteração (k)');
ylabel('Erro');
title('Evolução do Erro (Escala Logarítmica)');
legend('||X_{k+1} - X_k||_F', '||R_k||_F', 'Location', 'southwest');
% Subplot 2: Verificação da Convergência Quadrática
% Gráfico de R_{k+1} vs R_k numa escala log-log. Um declive de 2 prova a convergência quadrática.
subplot(1, 2, 2);
% Extrair R_k e R_{k+1} (ignorando valores zero/NaN no final)
R_k = ress_plot(1:end-1);
R_k1 = ress_plot(2:end);
% Traçar os dados experimentais
loglog(R_k, R_k1, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
hold on;
% Criar uma linha de referência teórica (Declive = 2 no espaço log-log)
% Estimamos a constante C usando as duas primeiras iterações
C_est = R_k1(1) / (R_k(1)^2);
ref_x = logspace(log10(min(R_k)), log10(max(R_k)), 20);
ref_y = C_est * (ref_x.^2);
% Traçar a linha de referência
loglog(ref_x, ref_y, 'b--', 'LineWidth', 1.5);
grid on;
xlabel('||R_k||_F');
ylabel('||R_{k+1}||_F');
title('Verificação Quadrática (R_{k+1} vs R_k)');
legend('Dados Experimentais', 'Referência Quadrática (Declive 2)', 'Location', 'northwest');
sgtitle('Análise de Convergência do Método de Newton-Schulz');
%2 – Decomposição QR de Householder
fprintf('\n Lehmer L_{100xn} \n');
fprintf('%-5s %-20s %-20s\n','n','||QR-L||_F','||Q''Q-I||_F');
for n = 10:10:100
    ri = (1:100).'; ci = 1:n;
    A  = min(ri,ci) ./ max(ri,ci);
    [Q, R] = house_qr(A);
    fprintf('%-5d %-20.6e %-20.6e\n', n, ...
            norm(Q*R-A,'fro'), norm(Q'*Q-eye(100),'fro'));
end
fprintf('\n Hilbert H_{100xn}\n');
fprintf('%-5s %-20s %-20s\n','n','||QR-H||_F','||Q''Q-I||_F');
for n = 10:10:100
    ri = (1:100).'; ci = 1:n;
    A  = 1./(ri + ci - 1);
    [Q, R] = house_qr(A);
    fprintf('%-5d %-20.6e %-20.6e\n', n, ...
            norm(Q*R-A,'fro'), norm(Q'*Q-eye(100),'fro'));
end
%3 – Compressão de Imagens via SVD
pct_keep = [1, 5, 10, 25, 50, 75];
%% a) Imagem a preto e branco
fprintf('\n 3a – Imagem grayscale \n');
[file_gray, path_gray] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.bmp', 'Imagens'}, 'Selecione a imagem a preto e branco');
if isequal(file_gray, 0), error('Seleção de imagem P&B cancelada.'); end
A_raw = imread(fullfile(path_gray, file_gray));
if ndims(A_raw)==3, A_raw = rgb2gray(A_raw); end
A = double(A_raw);
[m_g, n_g] = size(A);
[U,S,V] = svd(A,'econ');
sv       = diag(S);
r        = length(sv);
E_tot    = sum(sv.^2);
fprintf('Dimensão: %d x %d | Nº valores singulares: %d\n', m_g, n_g, r);
fprintf('%-10s %-8s %-12s %-14s\n','% SV kept','p','Qualidade %','Storage ratio');
figure('Name','Q3a – Grayscale');
subplot(2,4,1); imshow(uint8(A)); title('Original');
for idx = 1:length(pct_keep)
    p   = max(1, round(pct_keep(idx)/100*r));
    Ap  = U(:,1:p)*diag(sv(1:p))*V(:,1:p)';
    q   = sum(sv(1:p).^2)/E_tot*100;
    stor = (m_g + n_g + 1)*p / (m_g*n_g);  
    fprintf('%-10.0f %-8d %-12.4f %-14.6f\n', pct_keep(idx), p, q, stor);
    subplot(2,4,idx+1); imshow(uint8(Ap));
    title(sprintf('%d%% SV (q=%.1f%%)', pct_keep(idx), q));
end
sgtitle('Compressão SVD – Grayscale');
%% 3b) Imagem a cores
fprintf('\n 3b – Imagem a cores \n');
[file_color, path_color] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.bmp', 'Imagens'}, 'Selecione a imagem a cores');
if isequal(file_color, 0), error('Seleção de imagem a cores cancelada.'); end
A_raw = imread(fullfile(path_color, file_color));
if ndims(A_raw)~=3, error('Esperada imagem RGB.'); end
[m_c,n_c,~] = size(A_raw);
fprintf('%-10s %-8s %-12s %-12s %-12s %-12s\n', ...
        '% SV kept','p','Q_R %','Q_G %','Q_B %','Q_mean %');
figure('Name','Q3b – Colour');
subplot(2,4,1); imshow(A_raw); title('Original');
% svd de cada canal
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
%  QUESTÃO 4 – Número de condição via método das potências
% cond_2(A) = |lambda|_max / |lambda|_min   (A simétrica)
% |lambda|_max : método das potências
% |lambda|_min : método das potências inversas (LU, sem inverter A)
tol4    = 1e-10;
maxit4  = 5000;
%% 4b-i) Hilbert -------------------------------------------------------
fprintf('\n=== QUESTÃO 4b-i – Hilbert Hn ===\n');
fprintf('%-5s %-20s %-20s %-20s %-8s %-8s\n', ...
        'n','cond2 (potências)','cond2 (MATLAB)','Erro relativo','it_max','it_min');
fprintf('%s\n', repmat('-',1,88));
for n = 5:14
    A  = hilb(n);
    x0 = ones(n,1);
    [lam_max, ~, itM, ~] = pow_max(A, x0, tol4, maxit4);
    [lam_min, ~, itm, ~] = pow_min(A, x0, tol4, maxit4);
    c2    = abs(lam_max)/abs(lam_min);
    c2_ml = cond(A);
    fprintf('%-5d %-20.6e %-20.6e %-20.6e %-8d %-8d\n', ...
            n, c2, c2_ml, abs(c2-c2_ml)/c2_ml, itM, itm);
end
%% 4b-ii) Lehmer -------------------------------------------------------
fprintf('\n=== QUESTÃO 4b-ii – Lehmer Ln ===\n');
fprintf('%-6s %-20s %-20s %-20s %-8s %-8s\n', ...
        'n','cond2 (potências)','cond2 (MATLAB)','Erro relativo','it_max','it_min');
fprintf('%s\n', repmat('-',1,90));
for n = [10, 100, 200, 300, 400, 500]
    A  = lehmer(n);
    x0 = ones(n,1);
    [lam_max, ~, itM, ~] = pow_max(A, x0, tol4, maxit4);
    [lam_min, ~, itm, ~] = pow_min(A, x0, tol4, maxit4);
    c2 = abs(lam_max)/abs(lam_min);
    if n <= 200
        c2_ml = cond(A);
        fprintf('%-6d %-20.6e %-20.6e %-20.6e %-8d %-8d\n', ...
                n, c2, c2_ml, abs(c2-c2_ml)/c2_ml, itM, itm);
    else
        fprintf('%-6d %-20.6e %-20s %-20s %-8d %-8d\n', ...
                n, c2, '(omitido)', '(omitido)', itM, itm);
    end
end
%% 4c) Gráfico de convergência – exemplo com H_8 -----------------------
% Ilustra a convergência do método das potências (máximo e mínimo) para H_8.
fprintf('\n=== QUESTÃO 4 – Convergência para H_8 ===\n');
n_plt = 8;
A_plt = hilb(n_plt);
x0    = ones(n_plt,1);
[lam_max_plt, ~, itM_plt, lhist_max] = pow_max(A_plt, x0, tol4, maxit4);
[lam_min_plt, ~, itm_plt, lhist_min] = pow_min(A_plt, x0, tol4, maxit4);
figure('Name','Q4 – Convergence H_8');
subplot(1,2,1);
semilogy(1:itM_plt, abs(lhist_max(1:itM_plt) - lam_max_plt), 'b-o', 'MarkerSize', 4);
xlabel('Iteração'); ylabel('|\lambda^{(k)} - \lambda_{max}|');
title(sprintf('Método das potências – H_{%d}', n_plt)); grid on;
subplot(1,2,2);
semilogy(1:itm_plt, abs(lhist_min(1:itm_plt) - lam_min_plt), 'r-s', 'MarkerSize', 4);
xlabel('Iteração'); ylabel('|\lambda^{(k)} - \lambda_{min}|');
title(sprintf('Método das potências inversas – H_{%d}', n_plt)); grid on;
sgtitle('Convergência dos métodos das potências');
%% ================================================================
%  FUNÇÕES LOCAIS
%  (devem ficar após o fim do corpo do script)
%% ================================================================
% ----------------------------------------------------------------
function [X, iter, err_hist, res_hist] = ns_iterate(A, maxIter, tol)
% NS_ITERATE  Newton-Schulz iteration.  X_{k+1} = 2Xk - Xk*A*Xk
% Chute inicial: X0 = A' / ||A||_F^2   (garantia de convergência, ver 1b)
    n_   = size(A,1);
    I_   = eye(n_);
    X    = A.' / norm(A,'fro')^2;       % X0 da questão 1b
    err_hist = NaN(maxIter+1,1);
    res_hist = zeros(maxIter+1,1);
    res_hist(1) = norm(I_ - A*X,'fro');
    for k = 1:maxIter
        X_new = 2*X - X*A*X;
        step  = norm(X_new - X,'fro');
        res   = norm(I_ - A*X_new,'fro');
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
%
%   Algoritmo: para cada coluna k, constrói o vector de Householder
%   v = x + sign(x_1)*||x||*e_1 (o sinal evita cancelamento catastrófico),
%   normaliza v, e aplica H_k = I - 2vv' às submatrizes activas de R e Q^T.
    [m, n_] = size(A);
    R  = A;
    Qt = eye(m);      % acumula Q^T
    for k = 1:n_
        x   = R(k:m, k);
        s   = sign(x(1)); if s==0, s=1; end
        e1  = [1; zeros(m-k,1)];
        v   = x + s*norm(x)*e1;
        nv  = norm(v);
        if nv < eps*norm(x), continue; end
        v   = v/nv;
        R(k:m,  k:n_) = R(k:m,  k:n_) - 2*v*(v.'*R(k:m,  k:n_));
        Qt(k:m, :)    = Qt(k:m, :)     - 2*v*(v.'*Qt(k:m, :));
    end
    % anula ruído sub-diagonal
    for k=1:n_, R(k+1:m,k)=0; end
    Q = Qt.';
end
% ----------------------------------------------------------------
function [lam, x, iter, lam_hist] = pow_max(A, x0, tol, maxit)
% POW_MAX  Método das potências: valor próprio dominante de A.
%   Normalização por ||.||_inf; critério de paragem relativo em lambda.
    x        = x0 / norm(x0,inf);
    lam_prev = Inf;
    lam_hist = zeros(maxit,1);
    for k = 1:maxit
        z         = A*x;
        [~, i]    = max(abs(x));
        lam       = z(i)/x(i);
        lam_hist(k) = lam;
        x         = z / lam;        % normaliza: ||x||_inf = 1
        if abs(lam-lam_prev)/(abs(lam)+eps) < tol
            iter     = k;
            lam_hist = lam_hist(1:k);
            return;
        end
        lam_prev = lam;
    end
    iter     = maxit;
    lam_hist = lam_hist(1:maxit);
    warning('pow_max: max iterations reached.');
end
% ----------------------------------------------------------------
function [lam_min, x, iter, lam_hist] = pow_min(A, x0, tol, maxit)
% POW_MIN  Método das potências inversas: menor v.p. (em módulo) de A.
%   Resolve A*z = x via LU (sem calcular A^{-1} explicitamente).
    [L_, U_, P_] = lu(A);
    x         = x0 / norm(x0,inf);
    mu_prev   = Inf;
    lam_hist  = zeros(maxit,1);
    for k = 1:maxit
        z         = U_ \ (L_ \ (P_*x));
        [~, i]    = max(abs(x));
        mu        = z(i)/x(i);          % v.p. de A^{-1}
        lam_hist(k) = 1/mu;
        x         = z/mu;               % normaliza
        if abs(mu-mu_prev)/(abs(mu)+eps) < tol
            lam_min  = 1/mu;
            iter     = k;
            lam_hist = lam_hist(1:k);
            return;
        end
        mu_prev = mu;
    end
    lam_min  = 1/mu;
    iter     = maxit;
    lam_hist = lam_hist(1:maxit);
    warning('pow_min: max iterations reached.');
end
% ----------------------------------------------------------------
function print_ns_result(name, n, iter, errs, ress)
% PRINT_NS_RESULT  Imprime sumário de uma corrida do Newton-Schulz.
    fprintf('Matriz: %s  |  n=%d  |  iter=%d\n', name, n, iter);
    fprintf('  ||Xk+1-Xk||_F (final) = %.6e\n', errs(end));
    fprintf('  ||Rk||_F      (final) = %.6e\n', ress(end));
end