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

m_name = aux.m_name;
%m_wavelength_file_new = aux.m_wavelength_file_new;
info_slash = aux.info_slash;
n_m = aux.n_m;
path_source = aux.path_source;
path_target = aux.path_target;
subpath_tiff_dir = aux.subpath_tiff_dir;
subpath_jpg_dir = aux.subpath_jpg_dir;
subpath_matlab_dir = aux.subpath_matlab_dir;
m_wavelength_filepath = aux.m_wavelength_filepath;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength = repmat(aux.w_wavelength, 1, n_m); % hack Cambridge
m_wavelength_file_new = aux.m_wavelength_file; % cambridge hack

%%
% Make additional directories 
% for m = 1:aux.n_m
%     aux.subpath_tiff_mask_dir{m} = sprintf('%s%s_%s%s%s_%s+tiffm%s',...
%         aux.path_target,...
%         aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
%     aux.subpath_jpg_mask_dir{m} = sprintf('%s%s_%s%s%s_%s+jpgm%s',...
%         aux.path_target,...
%         aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
%     if ~exist(aux.subpath_tiff_mask_dir{m},'dir')
%         mkdir(aux.subpath_tiff_mask_dir{m})
%     end
%     if ~exist(aux.subpath_jpg_mask_dir{m},'dir')
%         mkdir(aux.subpath_jpg_mask_dir{m})
%     end
% end
% 
% subpath_tiff_mask_dir = aux.subpath_tiff_mask_dir;
% subpath_jpg_mask_dir = aux.subpath_jpg_mask_dir;

%% Check if reference value exists for all folios
for m = 1:n_m
    
    cd(path_target);
    D = dir();
    D = remove_hiddenfiles(D);
     
    filepath_reverse_tif =  sprintf('%s%s_DJK_reverse_true.tif',subpath_tiff_dir{m},m_name{m});
    filepath_reverse_jpg =  sprintf('%s%s_DJK_reverse_true.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_gray_reverse_tif =  sprintf('%s%s_DJK_reverse_gray.tif',subpath_tiff_dir{m},m_name{m});
    filepath_gray_reverse_jpg =  sprintf('%s%s_DJK_reverse_gray.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_reverse_tx_tif =  sprintf('%s%s_DJK_reverse_tx.tif',subpath_tiff_dir{m},m_name{m});

    filepath_true_reverse_reg = sprintf('%s%s/%s+tiff/%s_DJK_true_reverse_reg.tif', path_target, m_name{m}, m_name{m}, m_name{m});
    filepath_true_reverse_reg_jpg = sprintf('%s%s/%s+jpg/%s_DJK_true_reverse_reg.jpg', path_target, m_name{m}, m_name{m}, m_name{m});

    filepath_reverse = sprintf('%s%s/%s+tiff/%s_DJK_reverse_gray.tif',path_target, m_name{m}, m_name{m}, m_name{m}');
    filepath_out_reg_tiff = sprintf('%s%s/%s+tiff/%s_DJK_reverse_gray_reg.tif',path_target, m_name{m}, m_name{m}, m_name{m}');
    filepath_out_reg_jpg = sprintf('%s%s/%s+jpg/%s_DJK_reverse_gray_reg.jpg',path_target, m_name{m}, m_name{m}, m_name{m}');

    
    %filepath_reverse_mask_tif =  sprintf('%s%s_DJK_reverse_true.tif',subpath_tiff_mask_dir{m},m_name{m});
    %filepath_reverse_mask_jpg =  sprintf('%s%s_DJK_reverse_true.jpg',subpath_jpg_mask_dir{m},m_name{m});
    %filepath_reverse_mask_gray_tif =  sprintf('%s%s_DJK_reverse_gray.tif',subpath_tiff_mask_dir{m},m_name{m});
    %filepath_reverse_mask_gray_jpg =  sprintf('%s%s_DJK_reverse_gray.jpg',subpath_jpg_dir{m},m_name{m});
    if exist(filepath_true_reverse_reg, 'file') && exist(filepath_true_reverse_reg_jpg, 'file') &&...
             exist(filepath_gray_reverse_jpg, 'file') && exist(filepath_gray_reverse_tif, 'file')
        continue
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
    
    filepath_mask_front = sprintf('%s%s_DJK_mask.tif',subpath_matlab_dir{m}, m_name{m});
    
    if ~exist(filepath_mask_front, 'file')
        error('Please manually make a mask for parchment');
        %create_parchment_mask(aux);
    end 
    
    mask_front = imread(filepath_mask_front);
    filepath_mask_reverse = sprintf('%s%s/%s+matlab/%s_DJK_reverse_mask.tif', path_target, m_name{m}, m_name{m}, m_name{m});
    filepath_RGB_reverse = sprintf('%s%s/%s+tiff/%s_DJK_true.tif', path_target, m_name{m}, m_name{m}, m_name{m});
    path_reverse = sprintf('%s%s/', path_source, name_reverse);
    D = dir(path_reverse);
    D = remove_hiddenfiles(D);
    is_ir = cellfun(@(x) contains(x,'MB655DR'), D);
    if sum(is_ir) == 0
        is_ir = cellfun(@(x) contains(x,'MB700IR'), m_wavelength{m});
    end
    ix_ir = find(is_ir);
    ix_ir = ix_ir(1);
    
    is_red = cellfun(@(x) contains(x,'MB630RD'), D);
    ix_red = find(is_red);
    ix_red = ix_red(1);
    is_green = cellfun(@(x) contains(x,'MB530GN'), D);
    ix_green = find(is_green);
    ix_green = ix_green(1); 
    is_blue = cellfun(@(x) contains(x,'MB450RB'), D);
    ix_blue = find(is_blue);
    ix_blue = ix_blue(1);
    is_tx = cellfun(@(x) contains(x,'TX500CN'), D);
    ix_tx = find(is_tx);
    ix_tx = ix_tx(1);
    
    filepath_ir = sprintf('%s%s', path_reverse, D{ix_ir});
    filepath_red = sprintf('%s%s', path_reverse, D{ix_red});
    filepath_green = sprintf('%s%s', path_reverse, D{ix_green});
    filepath_blue = sprintf('%s%s', path_reverse, D{ix_blue});
    filepath_tx = sprintf('%s%s', path_reverse, D{ix_tx});

    I_ir = double(imread(filepath_ir));
    I_red = double(imread(filepath_red));
    I_green = double(imread(filepath_green));
    I_blue = double(imread(filepath_blue));
    I_tx = double(imread(filepath_tx));

    h = figure('name','Please choose spectralon');
    
    %imagesc(imadjust(I_red,stretchlim(I_red),[]));
    imagesc(I_red);
    
    hFH = imfreehand();
    % Create a binary image ("mask") from the ROI object.
    mask = hFH.createMask();
    delete(h);
    reference = zeros(5,1);
    spectralon_DC = I_red(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(1) =  1*LOW_HIGH(2)*spectral_DCmax;
    
    spectralon_DC = I_green(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(2) =  1*LOW_HIGH(2)*spectral_DCmax;
   
    spectralon_DC = I_blue(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(3) =  1*LOW_HIGH(2)*spectral_DCmax;
   % save(filepath_reference{m},'reference');
   
    spectralon_DC = I_ir(mask);
    spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    reference(4) =  1*LOW_HIGH(2)*spectral_DCmax;
   % save(filepath_reference{m},'reference');
    
   I_tx = I_tx./max(I_tx(:));
    LOW_HIGH = stretchlim(I_tx,[0 .99]);
    reference(5) =  LOW_HIGH(2);%1*LOW_HIGH(2)*spectral_DCmax;
   % save(filepath_reference{m},'reference');
    
   
    I_red = I_red./reference(1);
    I_green = I_green./reference(2);
    I_blue = I_blue./reference(3);
    I_ir = I_ir./reference(4);
    I_tx = I_tx./reference(5);

    % Get rotation
    got_rotation = false;
    % Get rotation
        if ~got_rotation
            got_rotation = true;
            filepath_nospace = strrep( filepath_red, ' ', '\ ');
            command = sprintf('%s -Orientation %s', aux.exiftoolcall, filepath_nospace);
            [~, exifout] = system(command);
             k = strfind(exifout, ': ');
            if numel(exifout)>0
                exifrotation = strtrim(exifout(k+2:end));
            end
            
            if numel(exifout)>0
                switch exifrotation
                    case 'Horizontal (normal)'
                        m_rotation_angle(m) = 0;
                    case 'Rotate 180'
                        m_rotation_angle(m) = 180;
                    case 'Rotate 90 CW'
                        m_rotation_angle(m) = 90;
                    case 'Rotate 270 CW'
                        m_rotation_angle(m) = 270;
                end
            end
            %filepath_rotation = sprintf('%srotation.mat', aux.subpath_matlab_dir{m});
            %rotation_angle = m_rotation_angle(m);
            %save(filepath_rotation, 'rotation_angle'); 
        end
        I_red = imrotate(I_red,-m_rotation_angle(m));
        I_green = imrotate(I_green,-m_rotation_angle(m));
        I_blue = imrotate(I_blue,-m_rotation_angle(m));
        I_ir = imrotate(I_ir,-m_rotation_angle(m));
        I_tx = imrotate(I_tx,-m_rotation_angle(m));

    RGB = zeros(size(I_red,1), size(I_red,2),3);
    RGB(:,:,1) = I_red;
    RGB(:,:,2) = I_green;
    RGB(:,:,3) = I_blue;
    
    RGB(RGB>1) = 1;
    I_ir(I_ir > 1) = 1;
    I_tx(I_tx > 1) = 1;

    RGB_tiff = uint16(RGB*65536);
    RGB_jpg = imresize(RGB,.4);
    RGB_jpg = uint8(RGB_jpg*256);
    IR_tiff = uint16(I_ir*65536);
    IR_jpg = imresize(I_ir,.4);
    IR_jpg = uint8(I_ir*256);
    Tx_tiff = uint16(I_tx*65536);
    %Tx_jpg = imresize(I_tx,.4);
    %Tx_jpg = uint8(I_tx*256);
    
    imwrite(RGB_tiff, filepath_reverse_tif, 'tif');
    imwrite(RGB_jpg, filepath_reverse_jpg, 'jpg', 'quality', 50);
    imwrite(IR_tiff, filepath_gray_reverse_tif, 'tif');
    imwrite(IR_jpg, filepath_gray_reverse_jpg, 'jpg', 'quality', 50);      
    imwrite(Tx_tiff, filepath_reverse_tx_tif, 'tif');
    %imwrite(Tx_jpg, filepath_reverse_tx_jpg, 'jpg', 'quality', 50);   
    
    if ~exist(filepath_mask_front, 'file')
        error('Please manually make a mask for reverse parchment');
      %create_parchment_mask(aux);
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
    
      figure; showMatchedFeatures(mask_front,mask_reverse_flip,matchedPoints1,matchedPoints2);
      legend('matched points 1','matched points 2');
    
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
    
      figure; showMatchedFeatures(mask_front,mask_reverse_flip,matchedPoints1,matchedPoints2);
      legend('matched points 1','matched points 2');
    
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
    %filepath_parchment_reverse = sprintf('%s%s%s%s+tiff%s%s_parchment_mask_reverse.tif', path_target, m_name{m}, info_slash, m_name{m},info_slash, m_name{m});
    filepath_true_reverse = sprintf('%s%s/%s+tiff/%s_DJK_reverse_true.tif', path_target, m_name{m}, m_name{m}, m_name{m});
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
    
    imwrite(Ir, filepath_true_reverse_reg);
    Ir_jpg = uint8(255*double(Ir)./65535);
    imwrite(Ir_jpg, filepath_true_reverse_reg_jpg);
    
    %imwrite(mask, filepath_parchment_reverse);
    
    % Repeat with mask
    Ir1 = Ir(:,:,1);
    Ir2 = Ir(:,:,2);
    Ir3 = Ir(:,:,3);
    Ir1(~mask) = 65535;
    Ir2(~mask) = 65535;
    Ir3(~mask) = 65535;
    I = cat(3,Ir1,Ir2,Ir3);
    imwrite(I, filepath_mask_reverse_tif);
    Jjpg = double(I)./65535;
    Jjpg = uint8(256*Jjpg);
    imwrite(Jjpg,filepath_reverse_mask_jpg,'jpeg','Quality', 50);

   % return
    %% Repeat for IR band 

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
    
    
    
%     filepath_gray_front = sprintf('%s%s',subpath_tiff_dir{m}, D2{ix_IR});
%     gray_front = imread(filepath_gray_front);
%     filepath_gray_front = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
%     if ~exist(filepath_gray_front, 'file');
%             filepath_gray_front = sprintf('%s%s_parchment_mask2.tif',subpath_tiff_dir{m},m_name{m});
%     end
%     gray_front = imread(filepath_gray_front);
%     filepath_gray_front = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
%     gray_front = imread(filepath_gray_front);
%     %I_reverse_flip = fliplr(imread('/Users/Kelbe/Desktop/MOTB/Processed/MOTB_MS_000566_r/MOTB_MS_000566_r+tiff/MOTB_MS_000566_r_parchment_mask.tif'));
%     figure;
%     showMatchedFeatures(gray_front,I_reverse_flip,...
%         inlierPtsOriginal,inlierPtsDistorted);
%     title('Matched inlier points - ud');
%     saveas(gcf, sprintf('%s%s',subpath_matlab_dir{m},'verso_registration4.jpg'));
%     close(gcf);
    
    imwrite(Ir, filepath_out_reg_tiff);
    Ir_jpg = uint8(255*double(Ir)./65535);
    
    imwrite(Ir_jpg, filepath_out_reg_jpg);
    
    % Repeat with mask
    
%     Ir(~mask) = 65535;
%     imwrite(Ir, filepath_reverse_mask_gray_tif);
%     Jjpg = double(Ir)./65535;
%     Jjpg = uint8(256*Jjpg);
%     imwrite(Jjpg,filepath_reverse_mask_gray_jpg,'jpeg','Quality', 50);
    
end
end





%end




