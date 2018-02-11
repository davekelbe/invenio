function [  ] = apply_mask_rgb( aux )
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
    filepath_tiff_mask = sprintf('%s%s_DJK_true_mask.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg_mask = sprintf('%s%s_DJK_true_mask.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    if exist(filepath_tiff_mask, 'file') && exist(filepath_jpg_mask, 'file')
        continue
    end
    filepath_tiff = sprintf('%s%s_DJK_true.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_mask = sprintf('%s%s_parchment_mask.tif',...
        subpath_tiff_dir{m}, m_name{m});
    I = imread(filepath_tiff);
    mask = imread(filepath_mask);
    I1 = I(:,:,1);
    I2 = I(:,:,2);
    I3 = I(:,:,3);
    I1(~mask) = 65535;
    I2(~mask) = 65535;
    I3(~mask) = 65535;
    I = cat(3,I1,I2,I3);
    imwrite(I, filepath_tiff_mask);
    Jjpg = double(I)./65535;
    Jjpg = uint8(256*Jjpg);
    imwrite(Jjpg,filepath_jpg_mask,'jpeg','Quality', 50);

    
    
end

end

%end




