function [  ] = register_verso_flipud( aux )
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

%% Check if reference value exists for all folios
for m = 1:n_m
    
    cd(path_target);
    D = dir();
    D = remove_hiddenfiles(D);
    
    filepath_reverse_tif =  sprintf('%s%s_DJK_reverse_true.tif',subpath_tiff_dir{m},m_name{m});
    filepath_reverse_jpg =  sprintf('%s%s_DJK_reverse_true.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_reverse_gray_tif =  sprintf('%s%s_DJK_reverse_gray.tif',subpath_tiff_dir{m},m_name{m});
    filepath_reverse_gray_jpg =  sprintf('%s%s_DJK_reverse_gray.jpg',subpath_jpg_dir{m},m_name{m});
    
    filepath_reverse_mask_tif =  sprintf('%s%s_DJK_reverse_true.tif',subpath_tiff_mask_dir{m},m_name{m});
    filepath_reverse_mask_jpg =  sprintf('%s%s_DJK_reverse_true.jpg',subpath_jpg_mask_dir{m},m_name{m});
    filepath_reverse_mask_gray_tif =  sprintf('%s%s_DJK_reverse_gray.tif',subpath_tiff_mask_dir{m},m_name{m});
    filepath_reverse_mask_gray_jpg =  sprintf('%s%s_DJK_reverse_gray.jpg',subpath_jpg_mask_dir{m},m_name{m});
    if exist(filepath_reverse_tif, 'file') && exist(filepath_reverse_jpg, 'file') && ...
             exist(filepath_reverse_mask_jpg, 'file') && exist(filepath_reverse_mask_tif, 'file') && ...
             exist(filepath_reverse_mask_gray_jpg, 'file') && exist(filepath_reverse_mask_gray_tif, 'file')
        % continue
    end
    %filepath_I_front = sprintf('%s%s_DJK_true.tif',subpath_tiff_dir{m}, m_name{m});
    %I_front = imread(filepath_I_front);
    suffix = m_name{m}(end);
    if strcmp(suffix,'X')
        new_suffix = 'Y';
    elseif strcmp(suffix,'Y')
        new_suffix = 'X';
    elseif contains(suffix, 'r') && strcmp(suffix(1),'_')
        new_suffix = sprintf('%sv', suffix(1));
    elseif contains(suffix, 'v') && strcmp(suffix(1),'_')
        new_suffix = sprintf('%sr', suffix(1));
    elseif contains(m_name{m}(end), 'v')
        new_suffix = 'r';
    elseif contains(m_name{m}(end), 'r')
        new_suffix = 'v';
    else
        cd(path_source);
        prompt_str = sprintf('%s%s%s','Please choose verso for ', m_name{m}, ' or press cancel');
        [verso] = uipickfiles('Prompt',prompt_str);
        new_suffix = verso{1}(end-1:end);
    end    
    
    nend = numel(new_suffix);
   
    str = sprintf('%s%s',m_name{m}(1:end-nend),new_suffix);
    if (strcmp(m_name{m}(end-2), 'A') || strcmp(m_name{m}(end-2), 'B')) && ...
            (strcmp(suffix, '_X') || strcmp(suffix, '_Y'))
        str_wild = sprintf('%s.%s',m_name{m}(1:end-3),new_suffix);
        ix_valid = cellfun(@(x) regexp(x,str_wild), D, 'uniformOutput', false);
        ix_valid = cellfun(@any, ix_valid);
        str = D{ix_valid};
    end
    
    new_suffix = strrep(m_name{m}, suffix, new_suffix);
    suffix = m_name{m};
    
%     is_reverse = cellfun(@(x) contains(x,str), D);
%     if ~sum(is_reverse)
%         return
%     end
    

    
    name_reverse = str;
    
    filepath_mask_front = sprintf('%s%s_parchment_mask.tif',subpath_matlab_dir{m}, m_name{m});
    mask_front = imread(filepath_mask_front);
    filepath_mask_reverse = sprintf('%s%s/%s+matlab/%s_parchment_mask.tif', path_target, name_reverse, name_reverse, name_reverse);
    filepath_RGB_reverse = sprintf('%s%s/%s+tiff/%s_DJK_true.tif', path_target, name_reverse, name_reverse, name_reverse);
    path_reverse = sprintf('%s%s/', path_source, name_reverse);
    D = dir(path_reverse);
    D = remove_hiddenfiles(D);
    is_ir = cellfun(@(x) contains(x,'MB655DR'), D);
    if sum(is_ir) == 0
        is_ir = cellfun(@(x) contains(x,'MB700IR'), m_wavelength{m});
    end
    ix_ir = find(is_ir);
    ix_ir = ix_ir(1);

    filepath_gray_reverse = sprintf('%s%s%s%s', path_target, name_reverse, info_slash, D{ix_ir});
    if ~exist(filepath_mask_reverse, 'file') || ~exist(filepath_RGB_reverse, 'file') || ~exist(filepath_gray_reverse, 'file');
        aux.m_name = {name_reverse};
        aux.m_path_upper = {strrep(m_path_upper{m}, suffix, new_suffix)};
        ix_delim = strfind(aux.m_name{m}, aux.options_delimiter);
        aux.m_mss = {aux.m_name{m}(1:ix_delim(end)-1)};
        aux.m_folio = {aux.m_name{m}(ix_delim(end)+1:end)};
        aux.path_tiff_dir = {strrep(aux.path_tiff_dir{m}, suffix, new_suffix)};
        aux.path_jpg_dir = {strrep(aux.path_jpg_dir{m}, suffix, new_suffix)};
        aux.path_tiff_mask_dir = {strrep(aux.path_tiff_mask_dir{m}, suffix, new_suffix)};
        aux.path_jpg_mask_dir = {strrep(aux.path_jpg_mask_dir{m}, suffix, new_suffix)};
        aux.path_matlab_dir = {strrep(aux.path_matlab_dir{m}, suffix, new_suffix)};
        aux.path_envi_dir = {strrep(aux.path_envi_dir{m}, suffix, new_suffix)};
        aux.m_wavelength_filepath = [];
        aux.m_wavelength_file = [];
        aux.m_wavelength = [];
        aux.m_wavelength_file_new = [];
        dir_local = sprintf('%s%s%s',path_source, new_suffix, info_slash);
        D = dir(dir_local);
        D = remove_hiddenfiles(D);
        
        for w = 1:numel(m_wavelength_filepath{m})
            aux.m_wavelength_filepath{1}{w} = strrep(m_wavelength_filepath{m}{w}, suffix, new_suffix);
            aux.m_wavelength_file{1}{w} = strrep(m_wavelength_file{m}{w}, suffix, new_suffix);
            aux.m_wavelength{1}{w} = strrep(m_wavelength{m}{w}, suffix, new_suffix);
            aux.m_wavelength_file_new{1}{w} = strrep(m_wavelength_file_new{m}{w}, suffix, new_suffix);
        end
        is_valid = false(numel(m_wavelength_file{1}),1);
        for d = 1:numel(D)
            is_valid_temp = cellfun(@(x) contains(x,D{d}), aux.m_wavelength_file{1})';
            is_valid = is_valid | is_valid_temp;
        end
        aux.m_wavelength_filepath{1} = aux.m_wavelength_filepath{1}(is_valid);
        aux.m_wavelength_file{1} = aux.m_wavelength_file{1}(is_valid);
        aux.m_wavelength{1} = aux.m_wavelength{1}(is_valid);
        aux.m_wavelength_file_new{1} = aux.m_wavelength_file_new{1}(is_valid);

        
        path_up = sprintf('%s%s%s',path_target, aux.m_name{1}, info_slash);
        if ~exist(path_up, 'dir')
            mkdir(path_up);
        end
        path_temp = sprintf('%s%s%s%s+tiff%s',path_target, aux.m_name{1}, info_slash, aux.m_name{1}, info_slash);
        if ~exist(path_temp, 'dir')
            mkdir(path_temp);
        end
        path_temp = sprintf('%s%s%s%s+tiffm%s',path_target, aux.m_name{1}, info_slash, aux.m_name{1}, info_slash);
        if ~exist(path_temp, 'dir')
            mkdir(path_temp);
        end
        path_temp = sprintf('%s%s%s%s+jpg%s',path_target, aux.m_name{1}, info_slash, aux.m_name{1}, info_slash);
        if ~exist(path_temp, 'dir')
            mkdir(path_temp);
        end
        path_temp = sprintf('%s%s%s%s+jpgm%s',path_target, aux.m_name{1}, info_slash, aux.m_name{1}, info_slash);
        if ~exist(path_temp, 'dir')
            mkdir(path_temp);
        end        
        aux.n_m = 1;
        create_spectralon_mask(aux);
        %create_chopsticks_mask(aux);
        %create_chopsticks2_mask(aux);
        create_overtext_mask(aux);
        create_parchment_mask(aux);
        reflectance_tiffs9_rgb(aux);
        reflectance_tiffs11(aux,'MB655DR');
    end
        
        
    mask_reverse = imread(filepath_mask_reverse);
    
    mask_reverse_flip = flipud(mask_reverse);
    
    points1 = detectSURFFeatures(mask_front);
    points2 = detectSURFFeatures(mask_reverse_flip);
    
    [f1,vpts1] = extractFeatures(mask_front,points1);
    [f2,vpts2] = extractFeatures(mask_reverse_flip,points2);
    
    indexPairs = matchFeatures(f1,f2) ;
    matchedPoints1 = vpts1(indexPairs(:,1));
    matchedPoints2 = vpts2(indexPairs(:,2));
    
    %  figure; showMatchedFeatures(mask_front,mask_reverse,matchedPoints1,matchedPoints2);
    %  legend('matched points 1','matched points 2');
    
    [tform_ud,inlierPtsDistorted_ud,inlierPtsOriginal_ud] = ...
        estimateGeometricTransform(matchedPoints2,matchedPoints1,...
        'similarity');
    %% Repeat  lr
    mask_reverse_flip = fliplr(mask_reverse);
    
    points1 = detectSURFFeatures(mask_front);
    points2 = detectSURFFeatures(mask_reverse_flip);
    
    [f1,vpts1] = extractFeatures(mask_front,points1);
    [f2,vpts2] = extractFeatures(mask_reverse_flip,points2);
    
    indexPairs = matchFeatures(f1,f2) ;
    matchedPoints1 = vpts1(indexPairs(:,1));
    matchedPoints2 = vpts2(indexPairs(:,2));
    
    %  figure; showMatchedFeatures(mask_front,mask_reverse,matchedPoints1,matchedPoints2);
    %  legend('matched points 1','matched points 2');
    
    [tform_lr,inlierPtsDistorted_lr,inlierPtsOriginal_lr] = ...
        estimateGeometricTransform(matchedPoints2,matchedPoints1,...
        'similarity');
    
    %% Compare reciprocal condition
    
    singularity = zeros(2,1);
    singularity(1) = rcond(tform_ud.T);
    singularity(2) = rcond(tform_lr.T);
    
    if singularity(1) > singularity(2)
        tform = tform_ud;
        inlierPtsDistorted = inlierPtsDistorted_ud;
        inlierPtsOriginal = inlierPtsOriginal_ud;
        flip = 'ud';
    else
        tform = tform_lr;
        inlierPtsDistorted = inlierPtsDistorted_lr;
        inlierPtsOriginal = inlierPtsOriginal_lr;
        flip = 'lr';
    end
    %%
    filepath_parchment_reverse = sprintf('%s%s%s%s+tiff%s%s_parchment_mask_reverse.tif', path_target, m_name{m}, info_slash, m_name{m},info_slash, m_name{m});
    filepath_true_reverse = sprintf('%s%s/%s+tiff/%s_DJK_true.tif', path_target, name_reverse, name_reverse, name_reverse);
    I_true_reverse = imread(filepath_true_reverse);
    switch flip
        case 'ud'
            I_true_reverse = cat(3, flipud(I_true_reverse(:,:,1)), flipud(I_true_reverse(:,:,2)), flipud(I_true_reverse(:,:,3)));
            mask_reverse = flipud(mask_reverse);

    case 'lr'
            I_true_reverse = cat(3, fliplr(I_true_reverse(:,:,1)), fliplr(I_true_reverse(:,:,2)), fliplr(I_true_reverse(:,:,3)));
            mask_reverse = fliplr(mask_reverse);
    end
    
    
    outputView = imref2d(size(mask_front));
    Ir = imwarp(I_true_reverse,tform,'OutputView',outputView);
    %  figure; imshow(Ir);
    %  title('Recovered image');
    
    outputViewMask = imref2d(size(mask_front));
    mask = imwarp(mask_reverse,tform,'OutputView',outputViewMask);
    
    imwrite(Ir, filepath_reverse_tif);
    Ir_jpg = uint8(255*double(Ir)./65535);
    imwrite(Ir_jpg, filepath_reverse_jpg);
    
    imwrite(mask, filepath_parchment_reverse);
    
    % Repeat with mask
    Ir1 = Ir(:,:,1);
    Ir2 = Ir(:,:,2);
    Ir3 = Ir(:,:,3);
    Ir1(~mask) = 65535;
    Ir2(~mask) = 65535;
    Ir3(~mask) = 65535;
    I = cat(3,Ir1,Ir2,Ir3);
    imwrite(I, filepath_reverse_mask_tif);
    Jjpg = double(I)./65535;
    Jjpg = uint8(256*Jjpg);
    imwrite(Jjpg,filepath_reverse_mask_jpg,'jpeg','Quality', 50);

   % return
    %% Repeat for IR band 
    path_out = sprintf('%s%s/%s+tiff/',path_target, name_reverse, name_reverse');
    cd(path_out);
    D2 = dir();
    D2 = remove_hiddenfiles(D2);
    is_IR = cellfun(@(x) contains(x,'MB655DR'), D2);
    if sum(is_IR) == 0
        is_IR = cellfun(@(x) contains(x,'MB700IR'), m_wavelength{m});
    end
    ix_IR = find(is_IR);
    if isempty(ix_IR);
        reflectance_tiffs11(aux,'MB655DR');
    end
    ix_IR = ix_IR(end);
    
    filepath_reverse = sprintf('%s%s', path_out,D2{ix_IR});
    if ~exist(filepath_reverse, 'file');
        foo = 1;
    end
    
    I_reverse = imread(filepath_reverse);
    switch flip
        case 'ud'
            I_reverse_flip = flipud(I_reverse);
    case 'lr'
            I_reverse_flip = fliplr(I_reverse);
    end
    
    outputView = imref2d(size(mask_front));
    Ir = imwarp(I_reverse_flip,tform,'OutputView',outputView);
    %  figure; imshow(Ir);
    %  title('Recovered image');
    
    
    cd(subpath_tiff_dir{m});
    D2 = dir();
    D2 = remove_hiddenfiles(D2);
    is_IR = cellfun(@(x) contains(x,'MB655DR'), D2);
    if sum(is_IR) == 0
        is_IR = cellfun(@(x) contains(x,'MB700IR'), m_wavelength{m});
    end
    ix_IR = find(is_IR);
    ix_IR = ix_IR(end);
    
    
    filepath_gray_front = sprintf('%s%s',subpath_tiff_dir{m}, D2{ix_IR});
    gray_front = imread(filepath_gray_front);
    filepath_gray_front = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
    if ~exist(filepath_gray_front, 'file');
            filepath_gray_front = sprintf('%s%s_parchment_mask2.tif',subpath_tiff_dir{m},m_name{m});
    end
    gray_front = imread(filepath_gray_front);
    filepath_gray_front = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
    gray_front = imread(filepath_gray_front);
    %I_reverse_flip = fliplr(imread('/Users/Kelbe/Desktop/MOTB/Processed/MOTB_MS_000566_r/MOTB_MS_000566_r+tiff/MOTB_MS_000566_r_parchment_mask.tif'));
    figure;
    showMatchedFeatures(gray_front,I_reverse_flip,...
        inlierPtsOriginal,inlierPtsDistorted);
    title('Matched inlier points - ud');
    saveas(gcf, sprintf('%s%s',subpath_matlab_dir{m},'verso_registration4.jpg'));
    close(gcf);
    
    imwrite(Ir, filepath_reverse_gray_tif);
    Ir_jpg = uint8(255*double(Ir)./65535);
    
    imwrite(Ir_jpg, filepath_reverse_gray_jpg);
    
    % Repeat with mask
    
    Ir(~mask) = 65535;
    imwrite(Ir, filepath_reverse_mask_gray_tif);
    Jjpg = double(Ir)./65535;
    Jjpg = uint8(256*Jjpg);
    imwrite(Jjpg,filepath_reverse_mask_gray_jpg,'jpeg','Quality', 50);
    
end
end





%end




