%% script_part4.m
% Question 4(b) – Condition number of Hilbert and Lehmer matrices via the
%                 power method and inverse power method.
%
% Results are compared with MATLAB's built-in cond(A) for validation.

format long;
tol     = 1e-10;
maxIter = 5000;

%% =========================================================
%  4(b)(i)  Hilbert matrices  H_n,  n = 5, 6, 7, ...
%% =========================================================
fprintf('========== Hilbert matrices H_n ==========\n');
fprintf('%-5s  %-20s  %-20s  %-20s  %-10s  %-10s\n', ...
        'n', 'cond2 (power method)', 'cond2 (MATLAB)', ...
        'Relative error', 'iter_max', 'iter_min');
fprintf('%s\n', repmat('-', 1, 92));

ns_hilb = 5 : 14;    % increase upper limit if needed; Hilbert ill-cond. grows fast

for n = ns_hilb
    A  = hilb(n);
    x0 = ones(n, 1);           % initial vector (will be normalised inside)

    [c2_pm, lam_max, lam_min, it_max, it_min] = cond2_power(A, x0, tol, maxIter);
    c2_ml = cond(A);           % MATLAB reference

    rel_err = abs(c2_pm - c2_ml) / c2_ml;

    fprintf('%-5d  %-20.6e  %-20.6e  %-20.6e  %-10d  %-10d\n', ...
            n, c2_pm, c2_ml, rel_err, it_max, it_min);
end

%% =========================================================
%  4(b)(ii)  Lehmer matrices  L_n,  n = 10, 100, 200, 300, ...
%% =========================================================
fprintf('\n========== Lehmer matrices L_n ==========\n');
fprintf('%-6s  %-20s  %-20s  %-20s  %-10s  %-10s\n', ...
        'n', 'cond2 (power method)', 'cond2 (MATLAB)', ...
        'Relative error', 'iter_max', 'iter_min');
fprintf('%s\n', repmat('-', 1, 94));

ns_lehmer = [10, 100, 200, 300, 400, 500];

for n = ns_lehmer
    A  = build_Lehmer(n);
    x0 = ones(n, 1);

    [c2_pm, lam_max, lam_min, it_max, it_min] = cond2_power(A, x0, tol, maxIter);

    % cond(A) is expensive for large n; skip for n > 200 if slow
    if n <= 200
        c2_ml   = cond(A);
        rel_err = abs(c2_pm - c2_ml) / c2_ml;
        fprintf('%-6d  %-20.6e  %-20.6e  %-20.6e  %-10d  %-10d\n', ...
                n, c2_pm, c2_ml, rel_err, it_max, it_min);
    else
        fprintf('%-6d  %-20.6e  %-20s  %-20s  %-10d  %-10d\n', ...
                n, c2_pm, '(skipped)', '(skipped)', it_max, it_min);
    end
end

%% ---- Convergence history plot (example: H_8) ----------------------------
n  = 8;
A  = hilb(n);
x0 = ones(n, 1);

[~, ~, ~, ~, lambda_hist_max] = power_method(A, x0, tol, maxIter);
% Note: power_method returns lambda_hist as 4th output
[lam_max, ~, it_max, lam_hist_max] = power_method(A, x0, tol, maxIter);
[lam_min, ~, it_min, lam_hist_min] = inverse_power_method(A, x0, tol, maxIter);

figure;
subplot(1,2,1);
semilogy(1:it_max, abs(lam_hist_max - lam_max), 'b-o');
xlabel('Iteration'); ylabel('|\lambda^{(k)} - \lambda_{max}|');
title(sprintf('Power method – H_{%d}', n)); grid on;

subplot(1,2,2);
semilogy(1:it_min, abs(lam_hist_min - lam_min), 'r-s');
xlabel('Iteration'); ylabel('|\lambda^{(k)} - \lambda_{min}|');
title(sprintf('Inverse power method – H_{%d}', n)); grid on;

sgtitle('Convergence of power methods');
