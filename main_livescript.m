%% Projeto ALN 2025/26 – Livescript Principal  (versão completa)
% Questões 1 a 4 do projeto de Álgebra Linear Numérica – LMAC 2025/26.
% Para gerar o relatório: Export > PDF no MATLAB Live Editor.

format long;

%% ============================================================
%  QUESTÃO 1 – Método de Newton-Schulz
%% ============================================================
%
% Iteração:  X_{k+1} = 2*X_k - X_k * A * X_k
%
% --- 1a) Rk+1 = Rk^2 ---------------------------------------------------
%
%  Rk+1 = I - A*Xk+1
%       = I - A*(2Xk - Xk*A*Xk)
%       = I - 2*A*Xk + A*Xk*A*Xk
%       = (I - A*Xk)^2  =  Rk^2
%
%  Portanto  ||Rk+1|| <= ||Rk||^2.  Se rho(R0) < 1, ||Rk|| -> 0
%  quadraticamente.
%
% --- 1b) X0 = A'/||A||_F^2 satisfaz rho(I - A*X0) < 1 -----------------
%
%  A*X0 = A*A'/||A||_F^2.  Os v.p. de I - A*X0 são 1 - sigma_i^2/||A||_F^2.
%  Como ||A||_F^2 = sum(sigma_i^2), cada sigma_i^2/||A||_F^2 in (0,1]
%  (A não singular => todos sigma_i > 0), logo cada v.p. in [0,1).
%  Portanto rho(I - A*X0) < 1.  QED.
%
% Funções: newton_schulz.m, ns_initial_guess.m, build_An.m,
%          build_Lehmer.m, convergence_order.m

fprintf('\n=== QUESTÃO 1 – Newton-Schulz ===\n');
run('script_part1de.m');

%% ============================================================
%  QUESTÃO 2 – Decomposição QR de Householder
%% ============================================================
%
% Algoritmo (coluna k, k=1..n):
%   x  = R(k:m, k)
%   v  = x + sign(x1)*||x||*e1     (sinal evita cancelamento)
%   v  = v / ||v||
%   R(k:m, k:n) -= 2*v*(v'*R(k:m, k:n))
%   Q(k:m, :)   -= 2*v*(v'*Q(k:m, :))   % acumula Q^T
%  No fim:  Q = Q^T transposta  =>  A = Q*R
%
% Função: householder_qr.m

fprintf('\n=== QUESTÃO 2 – QR de Householder ===\n');
run('script_part2.m');

%% ============================================================
%  QUESTÃO 3 – Compressão de Imagens via SVD
%% ============================================================
%
% Qualidade:  sigma = (sum sigma_i^2, i=1..p) / (sum sigma_i^2, i=1..r)
%
% 3a) Imagem a preto e branco
%   A = double(imread('image_gray.jpg'))
%   [U, S, V] = svd(A, 'econ')
%   Reconstrução com p valores singulares: U(:,1:p)*S(1:p,1:p)*V(:,1:p)'
%
% 3b) Imagem a cores: aplicar 3a) a cada canal R, G, B separadamente.
%
% Scripts: script_part3_gray.m, script_part3_color.m
%
% NOTA: Substitua 'image_gray.jpg' e 'image_color.jpg' pelos seus ficheiros.

fprintf('\n=== QUESTÃO 3 – Compressão SVD ===\n');
run('script_part3_gray.m');
run('script_part3_color.m');

%% ============================================================
%  QUESTÃO 4 – Número de condição via método das potências
%% ============================================================
%
% Para A simétrica:  cond2(A) = |lambda|_max / |lambda|_min
%
% Método das potências (dominante):
%   z(k) = A*x(k-1);  lambda(k) = z(k)_i / x(k-1)_i;  x(k) = z(k)/lambda(k)
%   Parar quando |lambda(k)-lambda(k-1)| / |lambda(k)| < tol
%
% Método das potências inversas (mínimo):
%   Igual mas aplicado a A^{-1}: resolve A*z = x por factorização LU.
%   NUNCA se calcula A^{-1} explicitamente.
%
% Funções: power_method.m, inverse_power_method.m, cond2_power.m

fprintf('\n=== QUESTÃO 4 – Número de condição ===\n');
run('script_part4.m');
