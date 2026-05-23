%% script_part3_color.m
% Question 3(b) – SVD-based compression of a colour (RGB) image.
%
% Each of the three channels (R, G, B) is compressed independently.
% Replace 'image_color.jpg' with your actual filename.

format long;

%% --- Load image -----------------------------------------------------------
img_file = 'image_color.jpg';         % <-- change to your file name
A_raw    = imread(img_file);

if ndims(A_raw) ~= 3 || size(A_raw, 3) ~= 3
    error('Expected an RGB image with 3 channels.');
end

[m, n, ~] = size(A_raw);
fprintf('Image size: %d x %d (colour)\n\n', m, n);

% Split into channels (convert to double)
R_ch = double(A_raw(:,:,1));
G_ch = double(A_raw(:,:,2));
B_ch = double(A_raw(:,:,3));

%% --- SVD for each channel -------------------------------------------------
[U_R, S_R, V_R] = svd(R_ch, 'econ');
[U_G, S_G, V_G] = svd(G_ch, 'econ');
[U_B, S_B, V_B] = svd(B_ch, 'econ');

sv_R = diag(S_R);   sv_G = diag(S_G);   sv_B = diag(S_B);
r    = length(sv_R);
fprintf('Number of singular values per channel: %d\n\n', r);

energy_R = sum(sv_R.^2);
energy_G = sum(sv_G.^2);
energy_B = sum(sv_B.^2);

%% --- Compression loop -----------------------------------------------------
pct_keep = [1, 5, 10, 25, 50, 75];

fprintf('%-12s  %-6s  %-12s  %-12s  %-12s  %-12s\n', ...
        '% SV kept', 'p', 'Quality R%', 'Quality G%', 'Quality B%', 'Mean Q%');
fprintf('%s\n', repmat('-', 1, 72));

figure('Name', 'Colour SVD Compression');
subplot(2, 4, 1);
imshow(A_raw);
title('Original');

for idx = 1 : length(pct_keep)
    p = max(1, round(pct_keep(idx) / 100 * r));

    % Reconstruct each channel
    R_c = U_R(:,1:p) * diag(sv_R(1:p)) * V_R(:,1:p).';
    G_c = U_G(:,1:p) * diag(sv_G(1:p)) * V_G(:,1:p).';
    B_c = U_B(:,1:p) * diag(sv_B(1:p)) * V_B(:,1:p).';

    % Clamp to [0,255] before casting
    R_c = min(max(R_c, 0), 255);
    G_c = min(max(G_c, 0), 255);
    B_c = min(max(B_c, 0), 255);

    % Assemble compressed image
    A_compressed          = zeros(m, n, 3, 'uint8');
    A_compressed(:,:,1)   = uint8(R_c);
    A_compressed(:,:,2)   = uint8(G_c);
    A_compressed(:,:,3)   = uint8(B_c);

    % Quality per channel
    q_R = sum(sv_R(1:p).^2) / energy_R * 100;
    q_G = sum(sv_G(1:p).^2) / energy_G * 100;
    q_B = sum(sv_B(1:p).^2) / energy_B * 100;
    q_m = (q_R + q_G + q_B) / 3;

    fprintf('%-12.0f  %-6d  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n', ...
            pct_keep(idx), p, q_R, q_G, q_B, q_m);

    subplot(2, 4, idx + 1);
    imshow(A_compressed);
    title(sprintf('%d%% SV  (Q≈%.1f%%)', pct_keep(idx), q_m));
end
sgtitle('Colour image – SVD compression (per RGB channel)');
