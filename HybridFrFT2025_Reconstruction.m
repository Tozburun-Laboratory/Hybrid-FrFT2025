%%% Article %% Hybrid fractional Fourier transform for filtering and imaging in a single digital holography workflow 
%%% File %% Hybrid Fractional Fourier Transform Reconstruction Framework
%%% Author %% Müge Topcu
%%% Supervisor %% Serhat Tozburun
%%% Upload Date %% 09 Aug 2025
%%% Version Date %% 09 Aug 2025
%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FrFT Fractions
a_o = 0.5; % optical FrFT - Default 0 (Depends on data)
a_d = 0.5; % digital FrFT - Default 1 (Depends on user)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sign of Diffraction Order
sgn = '-';
if strcmp(sgn, '-')
    multiplier = -1;
else
    multiplier = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 635*10^(-9); % Wavelength of the light in meters
dx = 5.3*10^(-6); % Pixel size of the camera sensor in meters
dy = 5.3*10^(-6); % Pixel size of the camera sensor in meters
k = 2*pi/lambda;
%% READ HOLOGRAM
holo_file = '/Users/.....'; %INPUT HOLOGRAM


imageArray = imread(holo_file);
imageArray = imageArray(:,:,1);
input_image = imageArray;
            targetSize = [1024, 1024];
            win1 = centerCropWindow2d(size(input_image),targetSize);
            input_image = imcrop(input_image,win1);
            input_image = double(input_image);
holo = input_image;
holo_image = (holo - min(holo(:))) / (max(holo(:)) - min(holo(:))); %for display
figure;imagesc(holo_image),colormap(gray),daspect([1 1 1]), axis off,  title('Hologram Image');
%% HOLOGRAM PREPROCESS
close all;

%%%% DISPLAY SPECTRUM DOMAIN OF RAW HOLOGRAM
[N,M] = size(holo);    
            [X,Y] = meshgrid(-M/2+1:M/2,-N/2+1:N/2); 
FT_holo = fftshift(fft2(fftshift(holo))); % Spectrum Domain of Hologram for display
figure;imagesc(log(abs(FT_holo))),colormap(gray),daspect([1 1 1]), axis off; title(sprintf('k space of hologram a_o=%.2f, a_d=%.2f', a_o, 1));

%%%% DC SUPPRESION
laplacian_kernel = [0 -1 0; -1 4 -1; 0 -1 0]; % Define the Laplacian kernel
lholo = conv2(holo, laplacian_kernel, 'same'); % Apply the Laplacian kernel to the hologram % Convolution with 'same' to maintain original size
FT_holo = fftshift(fft2(fftshift(lholo)));  %Spectrum Domain of Laplacian Hologram for display
figure; imagesc(abs(lholo)); colormap gray; axis image; axis off; title('Laplacian Hologram');
figure; imagesc(log(abs(FT_holo))); colormap gray; axis image; axis off; title(sprintf('k space of laplacian hologram a_o=%.2f, a_d=%.2f', a_o, 1));
holo = adapthisteq(holo);
holo = lholo;

%%%% ZERO PADDING
padded_size = 2048; % Customize based on the input image size 
image_size = padded_size;
pad_x = (padded_size - size(holo, 1)) / 2; % Calculate the padding sizes
pad_y = (padded_size - size(holo, 2)) / 2; % Calculate the padding sizes
holo = padarray(holo, [floor(pad_x), floor(pad_y)], 0, 'both'); % Perform zero-padding (centered)
figure; imagesc(abs(holo)); colormap gray; axis image; axis off; title(' Zero Padded Laplacian Hologram');
figure; imagesc(log(abs(FT_holo))); colormap gray; axis image; axis off; title('Fourier Domain of Z-Laplacian Hologram');

%%%% Digital FrFT 
angles_d = [a_d, a_d];
frft_holo = frft22d((single(holo)), angles_d); %Digital FrFT
figure; imagesc(log(abs(frft_holo))); colormap(gray); axis image; axis off; title('FrFT Domain of Hybrid FrFT Hologram');
frft_holo_image = log(abs(frft_holo));
%% ROI Isolation
% close all;

%%%% Hybrid FrFT 
a_h = a_d + (multiplier*a_o);
display(a_h);
display(sgn);
roi_size = [600, 600]; % Example

%%%% Choose ROI Interactively
figure('Position', [100, 100, 1024, 1024]);
imagesc(log(abs(frft_holo))); colormap gray; axis image; title('Click the ROI center');
[x_click, y_click] = ginput(1);  % user clicks once
center_points = [x_click, y_click]; % Store the selected center point
x_center = center_points(1);
y_center = center_points(2);
close(gcf);  
x1 = max(1, round(x_center - (roi_size(2) - 1) / 2)); % Left boundary
y1 = max(1, round(y_center - (roi_size(1) - 1) / 2)); % Top boundary
x2 = min(size(frft_holo, 2), round(x_center + (roi_size(2) - 1) / 2)); % Right boundary
y2 = min(size(frft_holo, 1), round(y_center + (roi_size(1) - 1) / 2)); % Bottom boundary
refined_ROI = false(size(frft_holo));
refined_ROI(y1:y2, x1:x2) = true;

%%%% Display the selected ROI 
figure('Position', [100, 100, 1024, 1024]);
imagesc(log(abs(frft_holo)));
colormap gray;
axis image;
hold on;
plot(x_center, y_center, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
rectangle('Position', [x1, y1, x2 - x1 + 1, y2 - y1 + 1], 'EdgeColor', 'r', 'LineWidth', 2);
title('Selected ROI');
hold off;

%%%% FRFT FILTERING
filtered_frft = frft_holo .* refined_ROI; 
ROI_area = sum(refined_ROI(:)); 
fprintf('Refined ROI area: %d pixels\n', ROI_area); 
figure, imagesc(abs(filtered_frft).^2), colormap gray;
title(['FT Hologram filter with area in pixel^2 ' + string(ROI_area)]);
daspect([1 1 1]), axis off;


%%%% SHIFT TO SPECTRUM DOMAIN 
rec = frft22d(filtered_frft, [1 - a_h, 1 - a_h]);
figure, imagesc(log(abs(rec))), colormap gray, daspect([1 1 1]), axis off,; title('Spectrum Domain after hybrid order FrFT');

%%%% MOVE TO CENTER
spectrum1 = abs(rec); 
maximum = max(max(spectrum1)); % Compute the absolute value of the spectrum and find its maximum
[x0, y0] = find(spectrum1 == maximum, 1); % Get the first occurrence of the maximum value
N = size(rec, 1); % Assume square matrix (N x N)
spectrum2 = zeros(N, N);
shiftX = round(N / 2) - round(x0); % Calculate shifts relative to the center
shiftY = round(N / 2) - round(y0);
% Ensure shifts are valid integers
if ~isfinite(shiftX) || ~isfinite(shiftY)
    error('Invalid shifts computed. Check the input data.');
end
% Perform the shifting
for ii = 1:N
    for jj = 1:N
        newX = ii + shiftX;
        newY = jj + shiftY;
        
        % Check if the new indices are within bounds
        if newX > 0 && newX <= N && newY > 0 && newY <= N
            spectrum2(newX, newY) = rec(ii, jj);
        end
    end
end
figure, imagesc(log(abs(spectrum2))), colormap gray, daspect([1 1 1]), axis off,; title('Centered Spectrum Domain after inverse hybrid order FrFT');

%%%% INVERSE FRFT
angles3 = [1, 1];
rec = frft22d(spectrum2, -angles3);
rec_amplitude = abs(rec);
figure, imagesc((abs(rec))), colormap gray, daspect([1 1 1]), axis off,; title('Amplitude Reconstruction after IFFT');
figure, imagesc(log(abs(rec))), colormap gray, daspect([1 1 1]), axis off,; title('Log Amplitude Reconstruction after IFFT');
figure, imagesc(angle(rec)), colormap gray, daspect([1 1 1]), axis off,; title('Phase Reconstruction IFFT');
rec_phase = angle(rec);
rec_amplitude = abs(rec);

%% PROPAGATION
% close all;

%%%% Propagation parameters
z_prop_final = multiplier*-80e-3;% z_prop = 0.0143; % Propagation distance (meters) (Adjusted based on the setup)
[M, N] = size(rec);
kx = 2 * pi * (-N/2:N/2-1) / (N * dx); % Spatial frequency in x
ky = 2 * pi * (-M/2:M/2-1) / (M * dy); % Spatial frequency in y
[Kx, Ky] = meshgrid(kx, ky);
H = exp(1i * z_prop_final * sqrt(k^2 - Kx.^2 - Ky.^2)); % Angular spectrum transfer function

%%%% Apply transfer function to spectrum
spectrum_centered = fftshift(fft2(fftshift(rec)));
propagated_spectrum = spectrum_centered .* H;

%%%% Inverse Fourier transform to get propagated field
propagated_field = fftshift(ifft2(fftshift(propagated_spectrum)));
figure, imagesc(abs(propagated_field)), colormap gray, daspect([1 1 1]), axis off,; %title(sprintf('Amplitude Reconstruction a_o = 0.%01d', a_o*10));
figure, imagesc(angle(propagated_field)), colormap gray, daspect([1 1 1]), axis off,; %title(sprintf('Phase Reconstruction a_o = 0.%01d', a_o*10));
rec_phase_prop = angle(propagated_field);
rec_amplitude_prop = abs(propagated_field);

%% MOVE TO CENTER
threshold = 1.3
[refined_ROI, bounding_box, center_point] = simple_rect_roi(propagated_field, threshold);
filtered_rec = propagated_field.*refined_ROI;
x0 = center_point(2);
y0 = center_point(1);
N = size(propagated_field, 1);
spectrum3 = zeros(N, N);
shiftX = round(N / 2) - round(x0);
shiftY = round(N / 2) - round(y0);
% Ensure shifts are valid integers
if ~isfinite(shiftX) || ~isfinite(shiftY)
    error('Invalid shifts computed. Check the input data.');
end
% Perform the shifting
for ii = 1:N
    for jj = 1:N
        newX = ii + shiftX;
        newY = jj + shiftY;

        % Check if the new indices are within bounds
        if newX > 0 && newX <= N && newY > 0 && newY <= N
            spectrum3(newX, newY) = filtered_rec(ii, jj);
        end
    end
end
spectrum_centered = spectrum3;
figure, imagesc(abs(spectrum_centered)), colormap gray, daspect([1 1 1]), axis off,; title(sprintf('Amplitude Reconstruction a_o=%.2f, a_d=%.2f', a_o, a_d));
figure, imagesc(angle(spectrum_centered)), colormap gray, daspect([1 1 1]), axis off,; title(sprintf('Phase Reconstruction a_o=%.2f, a_d=%.2f', a_o, a_d));
rec_phase_prop = angle(spectrum_centered);
rec_amplitude_prop = abs(spectrum_centered);
%% PHASE RECONSTRUCTION

%%%% PARAMETERS
targetSize = [1000, 1000];
N = targetSize(1);
M = targetSize(2);
phase_c = rec_phase_prop;
figure, imagesc(phase_c), colormap gray, daspect([1 1 1]), axis off; %title('Tilt Compensation (Phase)');

%%%% UNWRAPPING
[UP] = ls_unwrap(phase_c);
% figure; imagesc(UP); colormap gray; axis image; axis off; %title('Unwrapped Phase');
sx = 400;
cropped_phase_c = phase_c;
targetSize = [sx, sx];
win1 = centerCropWindow2d(size(cropped_phase_c),targetSize);
cropped_phase_c = imcrop(cropped_phase_c,win1);
% figure; imagesc(cropped_phase_c); colormap gray; axis square; axis off; %title('Cropped Initial Phase');
cropped_amp = rec_amplitude_prop;
targetSize = [sx, sx];
win1 = centerCropWindow2d(size(cropped_amp),targetSize);
cropped_amp = imcrop(cropped_amp,win1);
% figure; imagesc(cropped_amp); colormap gray; axis square; axis off; %title('Cropped Amplitude');
[CUP] = ls_unwrap(cropped_phase_c);
% figure; imagesc(CUP); colormap gray; axis image; axis off; %title('CUP');
% Optional: Remove any remaining tilt (if necessary)
[grad_x, grad_y] = gradient(UP);
avg_grad_x = mean(grad_x(:));
avg_grad_y = mean(grad_y(:));
[XX, YY] = meshgrid(1:min(size(UP)), 1:min(size(UP)));
tilt_correction = avg_grad_x * (XX - N/2) + avg_grad_y * (YY - M/2);
TUP = UP - tilt_correction;
% figure; imagesc(TUP); colormap gray; axis image; axis off; %title('TUP');
[K,L] = size(cropped_phase_c);    
[x1,Y1] = meshgrid(-L/2+1:L/2,-K/2+1:K/2); 
[grad_x, grad_y] = gradient(CUP);
avg_grad_x = mean(grad_x(:));
avg_grad_y = mean(grad_y(:));
[X1X1, Y1Y1] = meshgrid(1:K, 1:L);
tilt_correction = avg_grad_x * (X1X1 - K/2) + avg_grad_y * (Y1Y1 - L/2);
TCUP = CUP - tilt_correction;
% figure; imagesc(TCUP); colormap gray; axis image; axis off; %title('TCUP');
% Display 3D surface plots of the final unwrapped phases
figure;
subplot(2,2,1);
surf(UP);
colormap("jet")
shading(gca, 'interp');
axis square;
title('UP'); 
view(45,80);

subplot(2,2,2);
surf(TUP);
colormap("jet")
shading(gca, 'interp');
axis square;
title('TUP'); 
view(45,80);

subplot(2,2,3);
surf(CUP);
colormap("jet")
shading(gca, 'interp');
axis square;
title('CUP'); 
view(45,80);

subplot(2,2,4);
surf(TCUP);
colormap("jet")
shading(gca, 'interp');
axis square;
title('TCUP'); 
view(45,80);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%% RECTANGULAR ROI
function [rect_ROI, bounding_box, center_point] = simple_rect_roi(fft_holo,th_factor)
    % Compute the energy map from the FFT hologram.
    energy_map = log(abs(fft_holo).^2 + eps);
    energy_map = (energy_map - min(energy_map(:))) / (max(energy_map(:)) - min(energy_map(:)));
    threshold = graythresh(energy_map)*th_factor;
    % threshold = adaptthresh(energy_map,th_factor);
    % Create a binary map using the provided threshold.
    binary_map = energy_map > threshold;
    % binary_map = imbinarize(energy_map,threshold);
    % display(threshold);
    % Find connected components and select the largest one.
    cc = bwconncomp(binary_map);
    if cc.NumObjects < 1
        error('No region found above the threshold.');
    end
    num_pixels = cellfun(@numel, cc.PixelIdxList);
    [~, idx] = max(num_pixels);
    largest_component = false(size(binary_map));
    largest_component(cc.PixelIdxList{idx}) = true;
    
    % Get the bounding box of the largest connected component.
    stats = regionprops(largest_component, 'BoundingBox');
    bounding_box = stats.BoundingBox;  % Format: [x, y, width, height]
    
    % Create the rectangular ROI mask.
    rect_ROI = false(size(binary_map));
    x = round(bounding_box(1));
    y = round(bounding_box(2));
    width = round(bounding_box(3));
    height = round(bounding_box(4));
    rect_ROI(y:y+height-1, x:x+width-1) = true;
    
    % Compute the center of the ROI.
    center_point = [x + width/2, y + height/2];
    
    %%% Visualization
    figure('Name','Simple Rectangular ROI Detection','Position',[100, 100, 1200, 800]);
    
    % Panel 1: Energy Map
    subplot(2,2,1);
    imagesc(energy_map);
    colormap gray;
    axis image off;
    title('Energy Map');
    
    % Panel 2: Binary Map (Thresholded)
    subplot(2,2,2);
    imagesc(binary_map);
    colormap gray;
    axis image off;
    title(sprintf('Binary Map (Threshold = %.2f)', threshold));
    
    % Panel 3: Largest Component with Bounding Box
    subplot(2,2,3);
    imagesc(largest_component);
    colormap gray;
    axis image off;
    title('Largest Component');
    hold on;
    rectangle('Position', bounding_box, 'EdgeColor', 'r', 'LineWidth', 2);
    hold off;
    
    % Panel 4: Final ROI overlaid on Energy Map
    subplot(2,2,4);
    imagesc(energy_map);
    colormap gray;
    axis image off;
    title('Final ROI and Center');
    hold on;
    rectangle('Position', bounding_box, 'EdgeColor', 'g', 'LineWidth', 2);
    plot(center_point(1), center_point(2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
end

%%%% 2D PHASE UNWRAPPING
% Implements the “Goldstein‐Poisson” approach:
function phi = ls_unwrap(wrapped)
    [M,N] = size(wrapped);

    %--- a) compute forward wrapped gradients
    % pad a zero column/row so sizes match back‐diff below
    % horizontal gradient
    Wx = wrapToPi(diff(wrapped,1,2));
    Wx(:,N) = 0;
    % vertical gradient
    Wy = wrapToPi(diff(wrapped,1,1));
    Wy(M,:) = 0;

    %--- b) divergence of the wrapped gradient field
    % backward differences:
    dx = [Wx(:,1), diff(Wx,1,2)];    % ∂/∂x T
    dy = [Wy(1,:); diff(Wy,1,1)];    % ∂/∂y T
    div = dx + dy;

    %--- c) solve ∇² φ = div by DCT (Neumann BCs)
    % build cosine eigenvalues
    [u,v] = meshgrid(0:N-1,0:M-1);
    denom = (2*cos(pi*u/N)-2) + (2*cos(pi*v/M)-2);
    denom(1,1) = 1;  % avoid divide‐by‐zero for the DC term

    % DCT‐II of div
    div_hat = dct2(div);

    % spectral division
    phi_hat = div_hat ./ denom;

    % set the DC component to zero (mean‐free)
    phi_hat(1,1) = 0;

    % invert DCT
    phi = idct2(phi_hat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%