function [  ] = create_chopsticks2_mask( aux )
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
fprintf('Create chopsticks mask: \n');

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


chopsticks_spectra =   [0.0459
    0.3871
    0.7758
    0.1104
    1.4427];
chopsticks_wavelength = {'MB365UV', 'MB530GN', 'MB850IR', 'W365B47', 'W365O22'};

cd(m_path_upper{1});
D = dir('*_F.tif');
w_wavelength = remove_hiddenfiles(D);
n_w = numel(w_wavelength);
for w = 1:n_w
    ix_delimiter = strfind(w_wavelength{w},options_delimiter_wavelength);
    w_wavelength{w} = w_wavelength{w}(ix_delimiter:end);
end

                

%% Check if reference value exists for all folios
nref = numel(chopsticks_wavelength);
ref_val = zeros(nref,1);
for m = 1:n_m
    filepath_chopsticks_mask_tif =  sprintf('%s%s_chopsticks2_mask.tif',subpath_tiff_dir{m},m_name{m});
    filepath_chopsticks_mask_jpg =  sprintf('%s%s_chopsticks2_mask.jpg',subpath_jpg_dir{m},m_name{m});
    if exist(filepath_chopsticks_mask_tif, 'file') && exist(filepath_chopsticks_mask_jpg, 'file')
        continue 
    end
    mask_spectralon = imread(sprintf('%s%s_spectralon_mask.tif', subpath_tiff_dir{m}, m_name{m}));
    for r = 1:nref
        
        is_ref = cellfun(@(x) ~isempty(strfind(x,chopsticks_wavelength{r})), m_wavelength{m});
        ix_ref = find(is_ref);
        ix_ref = ix_ref(1);
        I = imread(m_wavelength_filepath{m}{ix_ref});
        if (m==1 && r==1)
            cube = zeros(size(I,1),size(I,2),nref);
            cube_ref = zeros(size(I,1),size(I,2),nref);
            sam = zeros(size(I,1),size(I,2));
        end
        maxVal = median(I(mask_spectralon));
        I = double(I)./double(maxVal);
        cube(:,:,r) = I; 
        cube_ref(:,:,r) = repmat(chopsticks_spectra(r),size(I));
    end
    sam = dot(cube,cube_ref,3);
    sam = sam./norm(chopsticks_spectra);
    cube_norm = cube(:,:,1).^2 + cube(:,:,2).^2 + cube(:,:,3).^2 + cube(:,:,4).^2 + cube(:,:,5).^2;
    cube_norm = power(cube_norm, 0.5);
    sam = acosd(sam./cube_norm);
    mask_sam = logical(sam < 0.5);
    mask_sam = imclose(mask_sam, strel('disk', 20));
    mask_sam = imdilate(mask_sam, strel('disk', 40));

    fprintf('                 \t\t%s\n', m_name{m});
    imwrite(mask_sam, filepath_chopsticks_mask_tif);
    imwrite(mask_sam, filepath_chopsticks_mask_jpg);

    %same = same.
   % for i = 1:size(I,1);
   %     for j = 1:size(I,2);
    %        same(i,j) = dot(cube(i,j), cube_ref(i,j)
end


    


%end




