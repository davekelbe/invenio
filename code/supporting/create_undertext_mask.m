function [  ] = create_undertext_mask( aux )
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
fprintf('Create overtext mask: \n');

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

info_min_pixels = 2000;
clear aux

%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});
for m = 1:n_m
    filepath_undertext_mask_jpg =  sprintf('%s%s_undertext_mask.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_undertext_mask_matlab=  sprintf('%s%s_undertext_mask.tif',subpath_matlab_dir{m},m_name{m});
    filepath_undertext_mask_tif =  sprintf('%s%s_undertext_mask.tif',subpath_tiff_dir{m},m_name{m});
    if (~exist(filepath_undertext_mask_jpg, 'file') || ~exist(filepath_undertext_mask_tif, 'file') || ~exist(filepath_undertext_mask_matlab, 'file'))
        
        is_dr = cellfun(@(x) contains(x,'W365R25'), m_wavelength{m});
        if sum(is_dr) == 0
            is_dr = cellfun(@(x) contains(x,'WBUVB47'), m_wavelength{m});
        end
        ix_dr = find(is_dr);
        ix_dr = ix_dr(end);
        
        % Load image and mask
        filepath_I = sprintf('%s%s', subpath_tiff_mask_dir{m}, m_wavelength_file_new{m}{ix_dr});
        I = imread(filepath_I);
        filepath_overtext_mask_tif =  sprintf('%s%s_overtext_mask.tif',subpath_tiff_dir{m},m_name{m});
        mask = imread(filepath_overtext_mask_tif);
        mask = imdilate(mask, strel('disk',8));
        I(mask) = 65535;
        
      %  filepath_parchment_mask=  sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
      %  mask_parchment = imread(filepath_parchment_mask);

        % Scale by parchment
        LOW_HIGH = stretchlim(I(I<65535), [0.01 0.99]);
        I_sc = imadjust(I, LOW_HIGH, []);
        mask_overtext = I_sc < 17000;
        mask_overtext = imopen(mask_overtext, strel('disk', 6));
        mask_overtext(mask) = 0;
        %I_tx = imadjust(I_tx ,stretchlim(I_tx), []);
        
        % Remove background 
        filepath_I = sprintf('%s%s_parchment_mask.tif', subpath_matlab_dir{m}, m_name{m});
        I = imread(filepath_I);   
        mask_overtext(~I) = 0;
        
        
        fprintf('                 \t\t%s\n', m_name{m});
        imwrite(mask_overtext, filepath_undertext_mask_tif);
        imwrite(mask_overtext, filepath_undertext_mask_matlab);
        imwrite(mask_overtext, filepath_undertext_mask_jpg);
    end
end



%end




