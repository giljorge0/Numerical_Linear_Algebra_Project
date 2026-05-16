function p = convergence_order(res_hist)
% CONVERGENCE_ORDER  Estimates the computational convergence order from a
%                    residual history using the formula
%
%       p_k ≈ log( ||R_{k+1}||_F / ||R_k||_F )
%             --------------------------------
%              log( ||R_k||_F / ||R_{k-1}||_F )
%
% Input:
%   res_hist - vector of residual norms ||R_k||_F, k = 0, 1, 2, ...
%              (at least 3 entries required)
%
% Output:
%   p - vector of estimated orders p_k for k = 1, ..., length-2

    m = length(res_hist);
    if m < 3
        p = NaN;
        return;
    end

    p = zeros(m - 2, 1);
    for k = 2 : m - 1
        num = log(res_hist(k+1) / res_hist(k));
        den = log(res_hist(k)   / res_hist(k-1));
        if abs(den) < eps
            p(k-1) = NaN;
        else
            p(k-1) = num / den;
        end
    end
end
