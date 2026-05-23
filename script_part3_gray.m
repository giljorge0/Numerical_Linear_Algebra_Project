%% script_part3_gray.m
% Question 3(a) – SVD-based compression of a grayscale image.
%
% The image file must be in the same folder (or on the MATLAB path).
% Replace 'image_gray.jpg' with your actual filename.
%
% The script:
%   1. Loads the image and converts it to a double matrix A.
%   2. Computes the full SVD of A.
%   3. For each target percentage of singular values kept (1%, 5%, 10%,
%      25%, 50%, 75%) reconstructs the compressed image and reports the
%      quality metric sigma.

format long;

%% --- Load image -----------------------------------------------------------
img_file = 'image_gray.jpg';          % <-- change to your file name
A_raw    = imread(img_file);

% If the image was loaded as RGB, convert to grayscale
if ndims(A_raw) == 3
    A_raw = rgb2gray(A_raw);
end

A = double(A_raw);                    % convert uint8 -> double
[m, n] = size(A);
fprintf('Image size: %d x %d\n', m, n);

%% --- Full SVD -------------------------------------------------------------
[U, S, V] = svd(A, 'econ');          % 'econ' gives min(m,n) singular values
sigma_all  = diag(S);                 % all singular values (descending order)
r          = length(sigma_all);
fprintf('Number of singular values: %d\n\n', r);

% Total energy (denominator of quality metric)
energy_total = sum(sigma_all.^2);

%% --- Compression loop -----------------------------------------------------
pct_keep = [1, 5, 10, 25, 50, 75];   % percentages of singular values to keep

fprintf('%-12s  %-10s  %-10s  %-20s\n', ...
        '% SV kept', 'p (SV)', 'Quality %', 'Storage ratio');
fprintf('%s\n', repmat('-', 1, 58));

figure('Name', 'Grayscale SVD Compression');
num_plots = length(pct_keep) + 1;
subplot(2, 4, 1);
imshow(uint8(A));
title('Original');

for idx = 1 : length(pct_keep)
    % Number of singular values to keep
    p = max(1, round(pct_keep(idx) / 100 * r));

    % Reconstruct: A_p = U(:,1:p) * S(1:p,1:p) * V(:,1:p)'
    A_compressed = U(:, 1:p) * diag(sigma_all(1:p)) * V(:, 1:p).';

    % Quality metric
    quality = sum(sigma_all(1:p).^2) / energy_total * 100;

    % Storage ratio: original needs m*n values; compressed needs (m+n+1)*p
    storage_ratio = (m + n + 1) * p / (m * n);

    fprintf('%-12.0f  %-10d  %-10.4f  %-20.6f\n', ...
            pct_keep(idx), p, quality, storage_ratio);

    % Show compressed image
    subplot(2, 4, idx + 1);
    imshow(uint8(A_compressed));
    title(sprintf('%d%% SV  (q=%.1f%%)', pct_keep(idx), quality));
end
sgtitle('Grayscale image – SVD compression');
