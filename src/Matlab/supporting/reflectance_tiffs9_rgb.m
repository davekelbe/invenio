function [  ] = reflectance_tiffs9_rgb( aux )
%NORMALIZED_ENVI_CUBE Create a normalized ENVI image cube
%
%   There is no input to this function. Typing reflectance_tiffs in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Reflectance TIFFS  Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%Greek
% V0.0 - Initial Version - January 4 2012
%
%
% Requirements:
%   *Commands are for UNIX and would need to be changed if used on a PC
%   *also requires these programs:
%       uipickfiles.m
%       binary_mask.m
%       combine_cube.m
%       enviwrite_bandnames.m
%
% Tips:
%   * Press ctrl+c to cancel execution and restart
%   *Set default paths in source code for efficiency
%% Preliminary setup
fprintf('\n***********************************************************\n');
fprintf('Truecolor RGB: \n');

m_name = aux.m_name;
n_m = aux.n_m;
subpath_tiff_dir = aux.subpath_tiff_dir;
subpath_jpg_dir = aux.subpath_jpg_dir;
subpath_matlab_dir = aux.subpath_matlab_dir;
w_wavelength = aux.w_wavelength;
m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;

%% Load shutter speed and aperture

is_red = cellfun(@(x) contains(x,'MB625Rd'), w_wavelength);
if sum(is_red) == 0
    is_red = cellfun(@(x) contains(x,'MB625RD'), w_wavelength);
end
if sum(is_red) == 0
    is_red = cellfun(@(x) contains(x,'MB630RD'), w_wavelength);
end
is_green = cellfun(@(x) contains(x,'MB535GN'), w_wavelength);
if sum(is_green) == 0
    is_green = cellfun(@(x) contains(x,'MB535GR'), w_wavelength);
end
if sum(is_green) == 0
    is_green = cellfun(@(x) contains(x,'MB535Gr'), w_wavelength);
end
if sum(is_green) == 0
    is_green = cellfun(@(x) contains(x,'MB530GN'), w_wavelength);
end

is_blue = cellfun(@(x) contains(x,'MB455RB'), w_wavelength);
if sum(is_blue) == 0
    is_blue = cellfun(@(x) contains(x,'MB450RB'), w_wavelength);
end

shutter_speed = zeros(n_m,3);
aperture = zeros(n_m,3);
for m = 1:n_m
    filepath_shutter_speed = sprintf('%s%s_shutter_speed.mat',subpath_matlab_dir{m},m_name{m});
    filepath_aperture = sprintf('%s%s_aperture.mat',subpath_matlab_dir{m},m_name{m});
    load(filepath_shutter_speed);
    load(filepath_aperture);
    
    shutter_speed(m,1) = w_shutter_speed(is_red);
    shutter_speed(m,2) = w_shutter_speed(is_green);
    shutter_speed(m,3) = w_shutter_speed(is_blue);
    aperture(m,1) = w_aperture(is_red);
    aperture(m,2) = w_aperture(is_green);
    aperture(m,3) = w_aperture(is_blue);
end
clear filepath_shutter_speed filepath_aperture m
% Output
% shutter_speed             - n_m x 3 (red, green, blue)
% aperture                  - n_m x 3 (red, green, blue)
%% Check if reference value exists for all folios
ref_exist = true;
filepath_reference = cell(n_m,1);
for m = 1:n_m
    filepath_reference{m} = sprintf('%s%s_rgb_reference.mat',subpath_matlab_dir{m},m_name{m});
    if ~exist(filepath_reference{m}, 'file')
        ref_exist = false;
    end
   filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    if ~exist(filepath_tiff, 'file')
        ref_exist = false;
    end
end
clear m
% Output
% ref_exist                 - true if all reference values already exists
%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});
for m = 1:n_m
    if ~ref_exist
    
    % Load first image for spectralon
    I_red = imread(m_wavelength_filepath{m}{is_red});
    I_red = double(I_red);
    h = figure('name','Please choose spectralon');
    
    %imagesc(imadjust(I_red,stretchlim(I_red),[]));
    imagesc(I_red);
    colormap(info_colormap);
    
    hFH = imfreehand();
    % Create a binary image ("mask") from the ROI object.
    mask = hFH.createMask();
    delete(h);
    reference = zeros(3,1);
    spectralon_DC = I_red(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(1) =  1*LOW_HIGH(2)*spectral_DCmax;
    
    I_green = imread(m_wavelength_filepath{m}{is_green});
    I_green = double(I_green);
    spectralon_DC = I_green(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(2) =  1*LOW_HIGH(2)*spectral_DCmax;
    
    I_blue = imread(m_wavelength_filepath{m}{is_blue});
    I_blue = double(I_blue);
    spectralon_DC = I_blue(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(3) =  1*LOW_HIGH(2)*spectral_DCmax;
    save(filepath_reference{m},'reference');
    
    
    %{
    ss_rep = repmat(shutter_speed(1,1),n_m,3);
    a_rep = repmat(aperture(1,1),n_m,3);
    exposure_factor = (shutter_speed./ss_rep).*(a_rep./aperture).^2;
    m_reference = exposure_factor*ref_val;
    for m = 1:n_m;
        reference = m_reference(m,:)';
        save(filepath_reference{m},'reference');
    end
    %}
    else
    end
    clear I h hfH mask spectralon_DC spectralon_DCmax LOW_HIGH
    clear ref_val ss_rep a_rep exposure_factor reference 
    % Output
    % m_reference               - reference value for reflectance calibration
    % For each folio, make truecolor RGB image
    
    fprintf('                 \t\t%s\n', m_name{m});
    
    filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg = sprintf('%s%s_DJK_true.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    if exist(filepath_tiff, 'file') && exist(filepath_jpg, 'file')
        continue
    end
    
    
    load(filepath_reference{m})
    %filepath_red = m_wavelength_filepath{m}{is_red};
    %filepath_green = m_wavelength_filepath{m}{is_green};
    %filepath_blue = m_wavelength_filepath{m}{is_blue};
    %I_red = double(imread(filepath_red));
    %I_green = double(imread(filepath_green));
    %I_blue = double(imread(filepath_blue));
    I_red = I_red./reference(1);
    I_green = I_green./reference(2);
    I_blue = I_blue./reference(3);
    
    RGB = zeros(size(I_red,1), size(I_red,2),3);
    RGB(:,:,1) = I_red;
    RGB(:,:,2) = I_green;
    RGB(:,:,3) = I_blue;
    
    RGB(RGB>1) = 1;
    
    % Get rotation
    filepath_rotation_angle = sprintf('%srotation.mat',...
        subpath_matlab_dir{m});
    load(filepath_rotation_angle);
                   
    RGB = imrotate(RGB,-rotation_angle);
    
    RGB_tiff = uint16(RGB*65536);
    RGB_jpg = imresize(RGB,.4);
    RGB_jpg = uint8(RGB_jpg*256);
    
    imwrite(RGB_tiff, filepath_tiff, 'tif');
    imwrite(RGB_jpg, filepath_jpg, 'jpg', 'quality', 50);
    
    fprintf('Writing %s\n', filepath_tiff);

end

%end




