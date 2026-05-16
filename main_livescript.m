%% Projeto ALN 2025/26 – Livescript Principal
% Este livescript serve como ponto de entrada para todo o projeto.
% Executa as duas partes e agrega os resultados de forma organizada.
% Para gerar relatório: Export > PDF no MATLAB Live Editor.

format long;

%% ============================================================
%  PARTE 1 – Método de Newton-Schulz
%% ============================================================
%
% Iteração:  X_{k+1} = 2*X_k - X_k * A * X_k
%
% Escolha inicial:  X_0 = A' / ||A||_F^2
%
% Critério de paragem combinado:
%   ||X_{k+1} - X_k||_F < tol   E   ||I - A*X_{k+1}||_F < tol
%
% Parâmetros globais
maxIter = 200;
tol     = 1e-12;

%% --- 1a) Prova teórica (resumo no comentário) ---
%
%  Definindo  R_k = I - A*X_k, temos:
%    R_{k+1} = I - A*X_{k+1}
%            = I - A*(2*X_k - X_k*A*X_k)
%            = I - 2*A*X_k + (A*X_k)^2
%            = (I - A*X_k)^2
%            = R_k^2
%
%  Portanto  ||R_{k+1}|| <= ||R_k||^2.
%  Se rho(R_0) = rho(I - A*X_0) < 1, então ||R_k|| -> 0 quadraticamente.

%% --- 1b) Prova que X_0 = A'/||A||_F^2 satisfaz a condição ---
%
%  A*X_0 = A*A'/||A||_F^2
%
%  Os valores próprios de A*A' são os sigma_i^2 (quadrados dos valores
%  singulares de A). Assim os valores próprios de I - A*X_0 são:
%    lambda_i = 1 - sigma_i^2 / ||A||_F^2
%
%  Como ||A||_F^2 = sum(sigma_i^2) >= sigma_i^2 para cada i,
%  temos  0 <= sigma_i^2/||A||_F^2 <= 1,  logo |lambda_i| < 1
%  (a igualdade apenas ocorreria se A tivesse posto 1, mas como A é
%  não singular com n >= 2 a desigualdade é estrita).
%  Portanto rho(I - A*X_0) < 1.

%% --- 1c) Função newton_schulz ---
% Ver ficheiro:  newton_schulz.m
% Protótipo:
%   [X, iter, err_hist, res_hist] = newton_schulz(A, X0, maxIter, tol)

%% --- 1d) Teste nas quatro famílias de matrizes ---
fprintf('===== PARTE 1d – Resultados =====\n\n');
run('script_part1de.m');

%% --- 1e) Ordem de convergência computacional ---
%  (calculada dentro de script_part1de.m via convergence_order.m)
%  Fórmula usada:
%       p_k ≈ ln( ||R_{k+1}||_F / ||R_k||_F )
%             ------------------------------------
%              ln( ||R_k||_F   / ||R_{k-1}||_F )
%
%  Espera-se p ≈ 2 (convergência quadrática), conforme demonstrado em 1a.

%% ============================================================
%  PARTE 2 – Decomposição QR de Householder
%% ============================================================
%
% Algoritmo:
%   Para k = 1,...,n:
%     1. x  = R(k:m, k)
%     2. v  = x + sign(x_1)*||x||*e_1   (escolha de sinal evita cancelamento)
%     3. v  = v / ||v||
%     4. R(k:m, k:n) -= 2*v*(v'*R(k:m, k:n))    % aplica H_k a R
%     5. Q(k:m, :)   -= 2*v*(v'*Q(k:m, :))       % acumula Q^T
%   No fim: Q = Q^T transposta.
%
%  Ver ficheiro:  householder_qr.m
%  Protótipo:
%    [Q, R] = householder_qr(A)

fprintf('\n\n===== PARTE 2 – Resultados =====\n\n');
run('script_part2.m');

%% Notas finais
%  - Para matrizes de Hilbert de grande dimensão, o elevado número de
%    condicionamento (cond(H_n) cresce exponencialmente) traduz-se em
%    erros de reconstrução ||QR-H||_F maiores, apesar do algoritmo de
%    Householder ser backward-stable.
%  - Para matrizes de Lehmer, o condicionamento cresce mais lentamente;
%    os erros de ortogonalidade e de reconstrução mantêm-se próximos de
%    eps_maq * n.
