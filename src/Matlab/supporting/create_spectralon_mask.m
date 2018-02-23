function [  ] = create_spectralon_mask( aux )
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
fprintf('Create spectralon mask: \n');

m_path_upper = aux.m_path_upper;
m_folio = aux.m_folio;
m_mss = aux.m_mss;
m_name = aux.m_name;
%m_wavelength_file_new = aux.m_wavelength_file_new;
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
subpath_matlab_dir = aux.path_matlab_dir;
subpath_envi_dir = aux.path_envi_dir;
%w_wavelength = aux.w_wavelength;
%w_wavelength = aux.w_wavelength;
%m_wavelength_file = aux.m_wavelength_file;
%m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;
m_wavelength_filepath = aux.m_wavelength_filepath;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength = aux.m_wavelength;
m_wavelength_file_new = aux.m_wavelength_file_new;

clear aux

exist_im = false(n_m,1);
for m = 1:n_m
    filepath_tiff = sprintf('%s%s_spectralon_mask.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg = sprintf('%s%s_spectralon_mask.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    filepath_matlab = sprintf('%s%s_spectralon_mask.tif',...
        subpath_matlab_dir{m}, m_name{m});
    if exist(filepath_tiff, 'file') && exist(filepath_jpg, 'file') && exist(filepath_matlab, 'file')
        exist_im(m) = true;
    end
end

if ~any(~exist_im)
    return %DJK REVERT
end

%% Check if reference value exists for all folios

%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});
for m = 1:n_m
    is_fluor = cellfun(@(x) ~isempty(strfind(x,'WBUVO22')), m_wavelength{m});
    if sum(is_fluor) == 0
            is_fluor = cellfun(@(x) ~isempty(strfind(x,'W365O22')), m_wavelength{m});
    end
    if sum(is_fluor) == 0
            is_fluor = cellfun(@(x) ~isempty(strfind(x,'TX450RB')), m_wavelength{m});
    end
    is_tx = cellfun(@(x) ~isempty(strfind(x,'TX500CN')), m_wavelength{m});
    if sum(is_tx) == 0
            is_tx = cellfun(@(x) ~isempty(strfind(x,'TX450RB')), m_wavelength{m});
    end
    if sum(is_tx) == 0
            is_tx = cellfun(@(x) ~isempty(strfind(x,'TX500CN')), m_wavelength{m});
    end
    is_refl = cellfun(@(x) ~isempty(strfind(x,'MB365UV')), m_wavelength{m});
    is_dr = cellfun(@(x) ~isempty(strfind(x,'MB655DR')), m_wavelength{m});
    if sum(is_dr) == 0
            is_dr = cellfun(@(x) ~isempty(strfind(x,'MB700IR')), m_wavelength{m});
    end
    ix_dr = find(is_dr);
    ix_dr = ix_dr(1);
    
    filepath_spectralon_mask_jpg =  sprintf('%s%s_spectralon_mask.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_spectralon_mask_matlab=  sprintf('%s%s_spectralon_mask.tif',subpath_matlab_dir{m},m_name{m});
    filepath_spectralon_mask_tif =  sprintf('%s%s_spectralon_mask.tif',subpath_tiff_dir{m},m_name{m});
    filepath_temp_jpg =  sprintf('%s%s_temp.jpg',subpath_jpg_dir{m},m_name{m});

    if true;%~exist(filepath_spectralon_mask_jpg, 'file') || ~exist(filepath_spectralon_mask_tif, 'file') || ~exist(filepath_spectralon_mask_matlab, 'file')
        %DJK REVERT
        
        % Load first image for spectralon
        filepath_fluor = sprintf('%s%s',m_path_upper{m},m_wavelength_file{m}{is_fluor});
        filepath_tx = sprintf('%s%s',m_path_upper{m},m_wavelength_file{m}{is_tx});
        ix = find(is_refl);
        %filepath_refl = sprintf('%s%s',m_path_upper{m},m_wavelength_file_new{m}{ix(1)});
        filepath_dr = sprintf('%s%s',m_path_upper{m},m_wavelength_file{m}{ix_dr});
        
        I_fluor = imread(filepath_fluor);
        I_tx = imread(filepath_tx);
        %I_refl = imread(filepath_refl);
        I_dr = imread(filepath_dr);
        
        LOW_HIGH = stretchlim(I_fluor, [0 0.999]);
        I_fluor = imadjust(I_fluor, LOW_HIGH, []);
        LOW_HIGH = stretchlim(I_tx, [0 0.999]);
        I_tx = imadjust(I_tx, LOW_HIGH, []);
        %figure; imagesc(I_fluor_stretch);
        LOW_HIGH = stretchlim(I_dr, [0 0.999]);
        I_dr = imadjust(I_dr, LOW_HIGH, []);
        
        %imwrite(uint8(255*double(I_fluor)./max(double(I_fluor(:)))), filepath_temp_jpg);
        
        mask_MegaVision = false(size(I_fluor));
        mask_MegaVision(I_fluor > 60000) = 1;
        mask_lightpanel = false(size(I_fluor));
        mask_lightpanel(I_tx > 10000) = 1;
        mask_spectralon = false(size(I_fluor));
        mask_spectralon(I_dr > 65000) = 1;
        mask_out = double(mask_spectralon & ~mask_MegaVision) - double(mask_lightpanel);
        mask_out = (mask_out >= 1);
        mask_out = imopen(mask_out, strel('disk', 30));
        
        path_temp = sprintf('%s%s%s', path_target, m_name{m}, info_slash);
        if ~exist(path_temp, 'dir');
            mkdir(path_temp);
        end
        path_temp = sprintf('%s%s%s%s+tiff%s', path_target, m_name{m}, info_slash,m_name{m}, info_slash);
        if ~exist(path_temp, 'dir');
            mkdir(path_temp);
        end
        path_temp = sprintf('%s%s%s%s+jpg%s', path_target, m_name{m},info_slash, m_name{m}, info_slash);
        if ~exist(path_temp, 'dir');
            mkdir(path_temp);
        end        
        path_temp = sprintf('%s%s%s%s+matlab%s', path_target, m_name{m}, info_slash, m_name{m}, info_slash);
        if ~exist(path_temp, 'dir');
            mkdir(path_temp);
        end        

        fprintf('                 \t\t%s\n', m_name{m});
        imwrite(mask_out, filepath_spectralon_mask_tif);
        imwrite(mask_out, filepath_spectralon_mask_matlab);
        imwrite(mask_out, filepath_spectralon_mask_jpg);
    end
end



%end




