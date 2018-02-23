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

m_path_upper = aux.m_path_upper;
m_folio = aux.m_folio;
m_mss = aux.m_mss;
m_name = aux.m_name;
is_band_subset = aux.is_band_subset;
bands = aux.bands;
info_rmcall = aux.info_rmcall;
info_slash = aux.info_slash;
info_user = aux.info_user;
n_m = aux.n_m;
options_delimiter = aux.options_delimiter;
options_delimiter_wavelength = aux.options_delimiter_wavelength;
options_folder_structure = aux.options_folder_structure;
options_movetonewfolder = aux.options_movetonewfolder;
path_source = aux.path_source;
path_target = aux.path_target;
subpath_tiff_dir = aux.path_tiff_dir;
subpath_jpg_dir = aux.path_jpg_dir;
subpath_tiff_mask_dir = aux.path_tiff_mask_dir;
subpath_jpg_mask_dir = aux.path_jpg_mask_dir;
subpath_matlab_dir = aux.path_matlab_dir;
subpath_envi_dir = aux.path_envi_dir;
m_wavelength = aux.m_wavelength;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength_file_new = aux.m_wavelength_file_new;
m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;

clear aux

exist_im = false(n_m,1);
for m = 1:n_m
    filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg = sprintf('%s%s_DJK_true.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    if exist(filepath_tiff, 'file') && exist(filepath_jpg, 'file')
        exist_im(m) = true;
    end
end

if ~any(~exist_im);
    %return
end


%% Load shutter speed and aperture



clear filepath_shutter_speed filepath_aperture m
% Output
% shutter_speed             - n_m x 3 (red, green, blue)
% aperture                  - n_m x 3 (red, green, blue)
%% Check if reference value exists for all folios
ref_exist = true;
filepath_reference = cell(n_m,1);
for m = 1:n_m;
    filepath_reference{m} = sprintf('%s%s_rgb_reference.mat',subpath_matlab_dir{m},m_name{m});
    if ~exist(filepath_reference{m}, 'file');
        ref_exist = false;
    end
    filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    if ~exist(filepath_tiff, 'file');
        ref_exist = false;
    end
end
clear m
% Output
% ref_exist                 - true if all reference values already exists
%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});
for m = 1:n_m;
    
    is_red1 = cellfun(@(x) ~isempty(strfind(x,'MB625Rd')), m_wavelength{m});
    is_red2 = cellfun(@(x) ~isempty(strfind(x,'MB630RD')), m_wavelength{m});
    is_red3 = cellfun(@(x) ~isempty(strfind(x,'MB625RD')), m_wavelength{m});
    is_red = (is_red1 | is_red2 | is_red3);
    ix_red = find(is_red);
    ix_red = ix_red(end);
    
    
    is_green1 = cellfun(@(x) ~isempty(strfind(x,'MB535GN')), m_wavelength{m});
    is_green2 = cellfun(@(x) ~isempty(strfind(x,'MB530GN')), m_wavelength{m});
    is_green3 = cellfun(@(x) ~isempty(strfind(x,'MB535GR')), m_wavelength{m});
    is_green4 = cellfun(@(x) ~isempty(strfind(x,'MB535Gr')), m_wavelength{m});
    is_green = (is_green1 | is_green2 | is_green3 | is_green4);
    ix_green = find(is_green);
    ix_green = ix_green(end);
    
    is_blue1 = cellfun(@(x) ~isempty(strfind(x,'MB455RB')), m_wavelength{m});
    is_blue2 = cellfun(@(x) ~isempty(strfind(x,'MB450RB')), m_wavelength{m});
    is_blue = (is_blue1 | is_blue2 );
    ix_blue = find(is_blue);
    ix_blue = ix_blue(end);
    
    filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg = sprintf('%s%s_DJK_true.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    filepath_tiff_mask = sprintf('%s%s_DJK_true_mask.tif',...
        subpath_tiff_mask_dir{m}, m_name{m});
    filepath_jpg_mask = sprintf('%s%s_DJK_true_mask.jpg',...
        subpath_jpg_mask_dir{m}, m_name{m});
    if ~exist(filepath_tiff, 'file') || ~exist(filepath_jpg, 'file') || ...
            ~exist(filepath_tiff_mask, 'file') || ~exist(filepath_jpg_mask, 'file')
        
        % Load first image for spectralon
        filepath_red = sprintf('%s%s', m_path_upper{m}, m_wavelength_file{m}{ix_red});
        I_red = imread(filepath_red);
        I_red = double(I_red);
        
        filepath_mask = sprintf('%s%s_spectralon_mask.tif',subpath_matlab_dir{m}, m_name{m});
        if ~exist(filepath_mask, 'file')
            h = figure('name','Please choose spectralon');
            
            %imagesc(imadjust(I_red,stretchlim(I_red),[]));
            imagesc(I_red);
            colormap(info_colormap);
            
            hFH = imfreehand();
            % Create a binary image ("mask") from the ROI object.
            mask = hFH.createMask();
            delete(h);
        else
            mask = imread(filepath_mask);
        end
        reference = zeros(3,1);
        spectralon_DC = I_red(mask);
        %spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
        %LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
        %reference(1) =  1*LOW_HIGH(2)*spectral_DCmax;
        reference(1) =  median(spectralon_DC);%1*LOW_HIGH(2)*spectral_DCmax;
        
        filepath_green = sprintf('%s%s', m_path_upper{m}, m_wavelength_file{m}{ix_green});
        I_green = imread(filepath_green);
        I_green = double(I_green);
        spectralon_DC = I_green(mask);
        %spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
        %LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
        %reference(2) =  1*LOW_HIGH(2)*spectral_DCmax;
        reference(2) =  median(spectralon_DC);%1*LOW_HIGH(2)*spectral_DCmax;
        
        filepath_blue = sprintf('%s%s', m_path_upper{m}, m_wavelength_file{m}{ix_blue});
        I_blue = imread(filepath_blue);
        I_blue = double(I_blue);
        spectralon_DC = I_blue(mask);
        %spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
        %LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
        %reference(3) =  1*LOW_HIGH(2)*spectral_DCmax;
        reference(3) =  median(spectralon_DC);%1*LOW_HIGH(2)*spectral_DCmax;
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
   % else
   % end
    clear I h hfH mask spectralon_DC spectralon_DCmax LOW_HIGH
    clear ref_val ss_rep a_rep exposure_factor
    % Output
    % m_reference               - reference value for reflectance calibration
    % For each folio, make truecolor RGB image
    
    
    
    
    %load(filepath_reference{m})
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
    %filepath_rotation_angle = sprintf('%srotation.mat',...
    %    subpath_matlab_dir{m});
    %load(filepath_rotation_angle);
    rotation_angle = 0;
    
    RGB = imrotate(RGB,-rotation_angle);
    
    RGB_tiff = uint16(RGB*65536);
    %RGB_jpg = imresize(RGB,.4);
    RGB_jpg = uint8(RGB*256);
    
    fprintf('                 \t\t%s\n', m_name{m});
    imwrite(RGB_tiff, filepath_tiff, 'tif');
    imwrite(RGB_jpg, filepath_jpg, 'jpg', 'quality', 50);
    
    
 
    filepath_mask = sprintf('%s%s_parchment_mask.tif',...
        subpath_matlab_dir{m}, m_name{m});
    mask = imread(filepath_mask);
    RGB_tiff1 = RGB_tiff(:,:,1);
    RGB_tiff2 = RGB_tiff(:,:,2);
    RGB_tiff3 = RGB_tiff(:,:,3);
    RGB_tiff1(~mask) = 65535;
    RGB_tiff2(~mask) = 65535;
    RGB_tiff3(~mask) = 65535;
    I = cat(3,RGB_tiff1,RGB_tiff2,RGB_tiff3);
    imwrite(I, filepath_tiff_mask);
    Jjpg = double(I)./65535;
    Jjpg = uint8(256*Jjpg);
    imwrite(Jjpg,filepath_jpg_mask,'jpeg','Quality', 50);
    
    
    
    end
end

%end




