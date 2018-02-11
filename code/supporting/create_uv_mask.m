function [  ] = create_uv_mask( aux )
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
fprintf('Create uv mask: \n');

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
subpath_matlab_dir = aux.path_matlab_dir;
subpath_envi_dir = aux.path_envi_dir;
%w_wavelength = aux.w_wavelength;
%w_wavelength = aux.w_wavelength;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;

info_min_pixels = 2000;
clear aux


%% Check if reference value exists for all folios

%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});

            TOL = [0 .99];
            sat = 0.75;
            
for m = 1:n_m
    is_dr = cellfun(@(x) contains(x,'W450O22'), m_wavelength_file{m});
    ix_dr = find(is_dr);
    ix_dr = ix_dr(1);
    filepath_uv_mask_jpg =  sprintf('%s%s_uv_mask.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_uv_mask_matlab=  sprintf('%s%s_uv_mask.tif',subpath_matlab_dir{m},m_name{m});
    filepath_uv_mask_tif =  sprintf('%s%s_uv_mask.tif',subpath_tiff_dir{m},m_name{m});
    if (~exist(filepath_uv_mask_jpg, 'file') || ~exist(filepath_uv_mask_tif, 'file') || ~exist(filepath_uv_mask_matlab, 'file'))
        
    
    % Load first image for spectralon
    I_dr = imread(m_wavelength_filepath{m}{ix_dr});
    mask_spectralon = imread(sprintf('%s%s_spectralon_mask.tif',subpath_tiff_dir{m},m_name{m})); 
    parch1d = I_dr(mask_spectralon);
        max_parch = max(parch1d);
        Jparch1d = parch1d./max_parch;
        upper_stretch1 = stretchlim(Jparch1d,TOL);
        upper_stretch2(2) = upper_stretch1(2)*(1/sat);
        stretchval =  upper_stretch2(2) * max_parch;
        clear mask parch1d max_parch Jparch1d upper_stretch1 upper_stretch2
        
        
        J = double(I_dr)./double(stretchval);

    J(J<=.45) = 0; 
    J(J>.45) = 1;
    
    J = logical(J);
    figure; imagesc(J)

%     [~, ix] = max(numPixels);
%     %isvalid = true(numel(numPixels),1);
%     %isvalid(ix) = false;
% 
%     for i = 1:ix-1
%     I(CC.PixelIdxList{i}) = 1;
%     end
%     for i = ix+1:numel(numPixels)
%     I(CC.PixelIdxList{i}) = 1;
%     end
%     
    % Dilate
      J = imclose(J, strel('disk', 20));

    
    fprintf('                 \t\t%s\n', m_name{m});
    imwrite(~J, filepath_uv_mask_tif);
    imwrite(~J, filepath_uv_mask_matlab);
    imwrite(~J, filepath_uv_mask_jpg);
    end
end
    


%end




