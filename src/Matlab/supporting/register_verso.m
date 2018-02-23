function [  ] = register_verso( aux )
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
fprintf('Register reverse: \n');

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
                
%% Check if reference value exists for all folios
for m = 1:n_m
    filepath_reverse_tif =  sprintf('%s%s_DJK_reverse.tif',subpath_tiff_dir{m},m_name{m});
    filepath_reverse_jpg =  sprintf('%s%s_DJK_reverse.jpg',subpath_jpg_dir{m},m_name{m});
    if exist(filepath_reverse_tif, 'file') && exist(filepath_reverse_jpg, 'file')
       % continue 
    end
    filepath_I_front = sprintf('%s%s_DJK_true.tif',subpath_tiff_dir{m}, m_name{m});
    I_front = imread(filepath_I_front);
    suffix = m_name{m}(end-1:end);
    if strcmp(suffix,'_X')
        new_suffix = '_Y';
    elseif strcmp(suffix,'_Y')
        new_suffix = '_X';
    end
    str = sprintf('%s%s',m_name{m}(1:end-2),new_suffix);
    cd(path_target);
    D = dir();
    D = remove_hiddenfiles(D);
    is_reverse = cellfun(@(x) ~isempty(strfind(x,str)), D);
    if ~sum(is_reverse)
        return
    end
    name_reverse = D{is_reverse};
    filepath_I_reverse = sprintf('%s%s/%s+tiff/%s_DJK_true.tif', path_target, name_reverse, name_reverse, name_reverse);
    I_reverse = imread(filepath_I_reverse);
    
    filepath_mask_front = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m}, m_name{m});
    mask_front = imread(filepath_mask_front);
    filepath_mask_reverse = sprintf('%s%s/%s+tiff/%s_parchment_mask.tif', path_target, name_reverse, name_reverse, name_reverse);
    mask_reverse = imread(filepath_mask_reverse);
    
    I_reverse = cat(3, fliplr(I_reverse(:,:,1)), fliplr(I_reverse(:,:,2)), fliplr(I_reverse(:,:,3)));
    mask_reverse = fliplr(mask_reverse);
    
    points1 = detectSURFFeatures(mask_front);
    points2 = detectSURFFeatures(mask_reverse);
    
    [f1,vpts1] = extractFeatures(mask_front,points1);
    [f2,vpts2] = extractFeatures(mask_reverse,points2);
    
    indexPairs = matchFeatures(f1,f2) ;
    matchedPoints1 = vpts1(indexPairs(:,1));
    matchedPoints2 = vpts2(indexPairs(:,2));
    
    %figure; showMatchedFeatures(mask_front,mask_reverse,matchedPoints1,matchedPoints2);
    %legend('matched points 1','matched points 2');

    [tform,inlierPtsDistorted,inlierPtsOriginal] = ...
    estimateGeometricTransform(matchedPoints2,matchedPoints1,...
    'similarity');


    if rcond(tform.T) < 1
       % fprintf('%s', 'Warning: Matrix is close to singular. Aborting')   
       % register_verso_flipud( aux );
      %  continue
    end
    
%     figure;
%     showMatchedFeatures(mask_front,mask_reverse,...
%     inlierPtsOriginal,inlierPtsDistorted);
%     title('Matched inlier points - lr');
    
    outputView = imref2d(size(mask_front));
    Ir = imwarp(I_reverse,tform,'OutputView',outputView);
   % figure; imshow(Ir);
   % title('Recovered image');
    
    imwrite(Ir, filepath_reverse_tif);
    Ir_jpg = uint8(255*double(Ir)./65535);

    imwrite(Ir_jpg, filepath_reverse_jpg);

end


    


%end




