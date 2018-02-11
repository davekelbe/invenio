function [ output_args ] = sam_classification( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

options_delimiter = '_';
options_delimiter_wavelength = '+';
options_folder_structure  = 'mss_folio';
options_movetonewfolder = 1;

mask = imread('/Users/Kelbe/Desktop/stick_mask.tif');
mask = double(mask);
mask = mask./max(mask(:));
mask = logical(mask);
% figure; imagesc(mask);

path_ref = '/Volumes/Seagate Backup Plus Drive/CAPTURE_SpectraHutt/FLATTENED/Assorted_Papyrus/GC_PAP_000531_2_Y/';
D = dir(path_ref);
w_file = remove_hiddenfiles(D);
w_wavelength = w_file;
n_w = numel(w_wavelength);

for w = 1:n_w
    ix_delimiter = strfind(w_wavelength{w},options_delimiter_wavelength);
    w_wavelength{w} = w_wavelength{w}(ix_delimiter:end);
end

% Load spectralon mask
mask_spectralon = imread('/Users/Kelbe/Desktop/MOTB/Processed/GC_PAP_000531_2_Y/GC_PAP_000531_2_Y+tiff/GC_PAP_000531_2_Y_spectralon_mask.tif');
% Load each image, normalize to spectralon, and save average value of mask
path_save = '/Users/Kelbe/Desktop/MOTB/Processed/GC_PAP_000531_2_Y/GC_PAP_000531_2_Y+norm/';
for w = 1:n_w
    I = imread(sprintf('%s%s', path_ref, w_file{w}));
    LOW_HIGH = stretchlim(I, [0 0.999]);
    I = imadjust(I, LOW_HIGH, []);
    imwrite(uint8(I*120), sprintf('%s%sjpg',path_save, w_file{w}(1:end-3)));
end

chopsticks_spectra = zeros(5,1);
I = imread(sprintf('%s%s', path_ref, 'GC_PAP_000531_2_Y+MB365UV_011_F.tif'));
maxVal = median(I(mask_spectralon));
I = double(I)./double(maxVal);
chopsticks_spectra(1) = median(I(mask));
I = imread(sprintf('%s%s', path_ref, 'GC_PAP_000531_2_Y+MB530GN_005_F.tif'));
maxVal = median(I(mask_spectralon));
I = double(I)./double(maxVal);
chopsticks_spectra(2) = median(I(mask));
I = imread(sprintf('%s%s', path_ref, 'GC_PAP_000531_2_Y+MB850IR_012_F.tif'));
maxVal = median(I(mask_spectralon));
I = double(I)./double(maxVal);
chopsticks_spectra(3) = median(I(mask));
I = imread(sprintf('%s%s', path_ref, 'GC_PAP_000531_2_Y+W365B47_021_F.tif'));
maxVal = median(I(mask_spectralon));
I = double(I)./double(maxVal);
chopsticks_spectra(4) = median(I(mask));
I = imread(sprintf('%s%s', path_ref, 'GC_PAP_000531_2_Y+W365O22_014_F.tif'));
maxVal = median(I(mask_spectralon));
I = double(I)./double(maxVal);
chopsticks_spectra(5) = median(I(mask));
save('/Users/Kelbe/Desktop/chopsticks_spectra.mat', 'chopsticks_spectra')
end

