%%% Article %% Hybrid fractional Fourier transform for filtering and imaging in a single digital holography workflow 
%%% File %% Hybrid Fractional Order Hologram Simulation 
%%% Author %% Müge Topcu
%%% Supervisor %% Serhat Tozburun
%%% Upload Date %% 09 Aug 2025
%%% Version Date %% 09 Aug 2025
%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FRFT HOLOGRAM SIMULATION 

%%%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FrFT Fractions
a_o = 0.5; % optical FrFT - Default 0 
a_d = 0.5; % digital FrFT - Default 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters on interference
z_distance = 100e-2; % Distance for reference wave in meters    
z_o_distance = 100e-2; %  Distance for object wave in meters
tilt_angle1 = 0.25*(pi / 2); % Angle between reference and object waves
tilt_angle2 = 0.383*(pi / 2); % Angle between reference and object waves
ratio = 1; % Intensity ratio between object beam and reference beam
% Fix parameters
lambda = 635e-9; % Wavelength of the light in meters
k = 2 * pi / lambda; % Wavenumber
dx = 5.3*10^(-6); % Pixel size of the camera sensor in meters
dy = 5.3*10^(-6); % Pixel size of the camera sensor in meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% READING THE IMAGE
imageArray = imread('USAF1951_1024.jpg');
imageArray = imageArray(:,:,1);
image_size = min(size(imageArray)); % Size of the hologram
input_image = imcrop(imageArray, centerCropWindow2d(size(imageArray), [image_size, image_size]));
input_image = double(input_image);


%%%% ZERO PADDING
% Define the desired padded size
padded_size = 2560; % Customize based on the input image size 
image_size = padded_size;
% Calculate the padding sizes
pad_x = (padded_size - size(input_image, 1)) / 2;
pad_y = (padded_size - size(input_image, 2)) / 2;
% Perform zero-padding (centered)
input_image = padarray(input_image, [floor(pad_x), floor(pad_y)], 0, 'both');
% figure; imagesc(input_image); colormap gray; axis image; axis off; title('Zero-Padded Input Image');


%%%% CONSTRUCT OBJECT BEAM
amp_image = input_image;
object_amplitude = (amp_image / max(amp_image(:)));
object_phase = pi/4 * ones(size(amp_image)); %plane
% figure; imagesc(object_amplitude); colormap gray; axis image; axis off; title('Input Image Amplitude');
% figure; imagesc(object_phase); colormap gray; axis image; axis off; title('Input Image Phase');
object_wave = object_amplitude .* exp(1i * object_phase);
% figure; imagesc(abs(object_wave)); colormap gray; axis image; axis off; title('Object Wave amplitude');
% figure; imagesc(angle(object_wave)); colormap gray; axis image; axis off; title('Object Wave Phase');
%%%% PROPAGATE OBJECT WAVE
[x, y] = meshgrid(1:image_size, 1:image_size);
x = (x - image_size / 2) * dx;
y = (y - image_size / 2) * dy;
r = sqrt(x.^2 + y.^2 + z_distance^2);
r_o = sqrt(x.^2 + y.^2 + z_o_distance^2);
object_beam = object_wave;%
angles_o = [a_o, a_o];
object_beam = frft22d((single(object_beam)), angles_o); %Optical FrFT
object_beam = object_beam/(max(abs(object_beam(:))))*ratio;


%%%% CONSTRUCT REFERENCE BEAM
reference_beam = exp(1i * k * (x * cos(tilt_angle1) + y * sin(tilt_angle2))); % Linear reference wave with tilt angle


%%%% CONSTRUCT HOLOGRAM
hologram_intensity = abs(object_beam + reference_beam).^2;
figure; imagesc(hologram_intensity); colormap gray; axis image; axis off; title('Simulated  Hologram');
FT_holo = fftshift(fft2(fftshift(hologram_intensity))); %Spectrum Domain of Hologram for display
figure; imagesc(log(abs(FT_holo))); colormap(gray); axis image; axis off; title(sprintf('k space of hologram a_o=%.2f, a_d=%.2f', a_o, 1));
holo = hologram_intensity;


%%%% DC SUPPRESION
laplacian_kernel = [0 -1 0; -1 4 -1; 0 -1 0]; % Define the Laplacian kernel
lholo = conv2(holo, laplacian_kernel, 'same'); % Apply the Laplacian kernel to the hologram % Convolution with 'same' to maintain original size
FT_holo = fftshift(fft2(fftshift(lholo)));  %Spectrum Domain of Laplacian Hologram for display
figure; imagesc(abs(lholo)); colormap gray; axis image; axis off; title('Laplacian Hologram');
figure; imagesc(log(abs(FT_holo))); colormap gray; axis image; axis off; title(sprintf('k space of laplacian hologram a_o=%.2f, a_d=%.2f', a_o, 1));
holo = lholo;


%%%% Digital FrFT 
angles_d = [a_d, a_d]; 
frft_holo = frft22d((single(holo)), angles_d); %Digital FrFT
figure; imagesc(log(abs(frft_holo))); colormap(gray); axis image; axis off;  title(sprintf('k space of hybrid frft hologram a_o=%.2f, a_d=%.2f', a_o, a_d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% TOZBURUN LABORATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
