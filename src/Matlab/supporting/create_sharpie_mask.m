function [  ] = create_sharpie_mask( aux )
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
fprintf('Create sharpie mask: \n');

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
    filepath_sharpie_jpg =  sprintf('%s%s_sharpie.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_sharpie_tif =  sprintf('%s%s_sharpie.tif',subpath_tiff_dir{m},m_name{m});
    filepath_sharpiem_jpg =  sprintf('%s%s_sharpie.jpg',subpath_jpg_mask_dir{m},m_name{m});
    filepath_sharpiem_tif =  sprintf('%s%s_sharpie.tif',subpath_tiff_mask_dir{m},m_name{m});
    filepath_sharpie_mask_jpg =  sprintf('%s%s_sharpie_mask.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_sharpie_mask_tif =  sprintf('%s%s_sharpie_mask.tif',subpath_tiff_dir{m},m_name{m});
    if true %(    ~exist(filepath_sharpie_jpg, 'file') || ~exist(filepath_sharpie_tif, 'file') || ... 
           % ~exist(filepath_sharpiem_jpg, 'file') || ~exist(filepath_sharpiem_tif, 'file') || ...
           % ~exist(filepath_sharpie_mask_jpg, 'file') || ~exist(filepath_sharpie_mask_tif, 'file'))
        
        is_ir = cellfun(@(x) contains(x,'MB940IR'), m_wavelength{m});
        ix_ir = find(is_ir);
        if isempty(ix_ir)
              is_ir = cellfun(@(x) contains(x,'MB850IR'), m_wavelength{m});
              ix_ir = find(is_ir);
        end
        ix_ir = ix_ir(end);
        is_uv = cellfun(@(x) contains(x,'W365O22'), m_wavelength{m});
        if sum(is_uv) == 0
            is_uv = cellfun(@(x) contains(x,'WBRBO22'), m_wavelength{m});
        end
        ix_uv = find(is_uv);
        ix_uv = ix_uv(end);
        
        % Load image and mask
        filepath_IR = sprintf('%s%s', subpath_tiff_dir{m}, m_wavelength_file_new{m}{ix_ir});
        IR = imread(filepath_IR);
        filepath_UV = sprintf('%s%s', subpath_tiff_dir{m}, m_wavelength_file_new{m}{ix_uv});
        UV = imread(filepath_UV);
      %  filepath_parchment_mask=  sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
      %  mask_parchment = imread(filepath_parchment_mask);
        
        sharpie = double(IR)-double(UV);
        sharpie = localnormalize(double(sharpie),201,201);
        
        sharpie = sharpie - min(sharpie(:));
        sharpie = sharpie./max(sharpie(:));
        LOW_HIGH = stretchlim(sharpie, .0025);
        sharpie = imadjust(sharpie, LOW_HIGH);
        sharpie = abs(1 - sharpie);
       
        sharpie_mask = (sharpie < 0.2747);    
        sharpie_mask2 = imopen(sharpie_mask, strel('disk', 3));
        
        sharpiem = sharpie;

        % Remove background 
        filepath_parchment = sprintf('%s%s_parchment_mask.tif', subpath_tiff_dir{m}, m_name{m});
        parchment = imread(filepath_parchment);   

        sharpie_mask(~parchment) = 0;
        sharpiem(~parchment) = 1;

        fprintf('                 \t\t%s\n', m_name{m});
        imwrite(uint16(65535*sharpie), filepath_sharpie_tif);
        imwrite(uint8(255*sharpie), filepath_sharpie_jpg);
        imwrite(uint16(65535*sharpiem), filepath_sharpiem_tif);
        imwrite(uint8(255*sharpiem), filepath_sharpiem_jpg);
        imwrite(sharpie_mask, filepath_sharpie_mask_tif);
        imwrite(sharpie_mask, filepath_sharpie_mask_jpg);
    end
end



%end




