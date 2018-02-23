function [  ] = create_parchment_mask( aux )
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
fprintf('Create parchment mask: \n');

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
m_wavelength = aux.m_wavelength;
%m_wavelength_file = aux.m_wavelength_file;
%m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;

info_min_pixels = 2000;
clear aux

cd(m_path_upper{1});
D = dir('*_F.tif');
w_wavelength = remove_hiddenfiles(D);
n_w = numel(w_wavelength);
for w = 1:n_w
    ix_delimiter = strfind(w_wavelength{w},options_delimiter_wavelength);
    w_wavelength{w} = w_wavelength{w}(ix_delimiter:end);
end

m_wavelength_filepath = cell(n_m,1);
m_wavelength_file = cell(n_m,1);
for m = 1:n_m
    
    % Determine wavelength filepath
    m_wavelength_filepath{m} = cell(n_w,1);
    m_wavelength_file{m} = cell(n_w,1);
    switch options_folder_structure
        case 'mss_folio'
            for w = 1:n_w;
                m_wavelength_filepath{m}{w} = sprintf('%s%s%s', ...
                    m_path_upper{m},m_name{m},w_wavelength{w});
                m_wavelength_file{m}{w} = sprintf('%s%s', ...
                    m_name{m},w_wavelength{w});
            end
    end
end
                
exist_im = false(n_m,1);
for m = 1:n_m
    filepath_tiff = sprintf('%s%s_parchment_mask.tif',...
        subpath_tiff_dir{m}, m_name{m});
    filepath_jpg = sprintf('%s%s_parchment_mask.jpg',...
        subpath_jpg_dir{m}, m_name{m});
    if exist(filepath_tiff, 'file') && exist(filepath_jpg, 'file')
        exist_im(m) = true;
    end
end

if ~any(~exist_im)
   % return
end

%% Check if reference value exists for all folios
is_tx = cellfun(@(x) contains(x,'TX500CN'), w_wavelength);
if sum(is_tx) == 0
    is_tx = cellfun(@(x) contains(x,'TX450RB'), m_wavelength{m});
end
is_uv = cellfun(@(x) contains(x,'W450O22'), w_wavelength);
if sum(is_uv) == 0
    is_uv = cellfun(@(x) contains(x,'WBUVO22'), m_wavelength{m});
end
is_ir = cellfun(@(x) contains(x,'TX940IR'), w_wavelength);
is_gr = cellfun(@(x) contains(x,'MB530GN'), w_wavelength);
if sum(is_gr) == 0
    is_gr = cellfun(@(x) contains(x,'MB535Gr'), m_wavelength{m});
end
ix_gr = find(is_gr);
ix_gr = ix_gr(1);
%% Find spectralon and make reference
%filepath_reference = sprintf('%srgb_reference.txt',subpath_matlab_dir{m});
for m = 1:n_m
    filepath_parchment_mask_jpg =  sprintf('%s%s_parchment_mask.jpg',subpath_jpg_dir{m},m_name{m});
    filepath_parchment_mask_matlab=  sprintf('%s%s_parchment_mask.tif',subpath_matlab_dir{m},m_name{m});
    filepath_parchment_mask_tif =  sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m},m_name{m});
    I_tx = imread(m_wavelength_filepath{m}{is_tx});

    if ( exist(filepath_parchment_mask_tif, 'file') && exist(filepath_parchment_mask_matlab, 'file'))% || ~exist(filepath_parchment_mask_jpg, 'file'))
        continue
    end
    
    fn_old = m_wavelength_filepath{m}{is_tx};
        
    % Load first image for spectralon
    LOW_HIGH = stretchlim(I_tx, [0 0.99]);
    I_tx = imadjust(I_tx, LOW_HIGH, []);
    
    fn_new =  sprintf('%s%s',subpath_tiff_dir{m},m_wavelength_file{m}{is_tx});
    imwrite(I_tx, fn_new);
    
    
    %I_tx = imadjust(I_tx ,stretchlim(I_tx), []);
    
    I = I_tx;
    %I(I_tx<=41768) = 0; 
    %I(I_tx>41768) = 1;
    I(I_tx<=54068) = 0; 
    I(I_tx>540968) = 1;
    
    I = logical(I);
    
    % More strict criteria for edges 
     edges = stdfilt(I);
     edges = edges > 0;
     se = strel('disk', 11);
     edges = imdilate(edges, se);
%     I(edges) = 0;
    
    %foo = I_tx(edges); 
    %x = linspace(double(min(foo)), double(max(foo)),500);
    %y = hist(foo, x);
    I_tx_edges = I_tx;
    %I_tx_edges(~edges) = 0;
    ix_step = 1:1000:size(I,1);
    jx_step = 1:1000:size(I,2);
    %I_tx_local = I_tx;
    for i = 1:numel(ix_step)-1
        for j = 1:numel(jx_step)-1
            I_local = false(size(I));
            I_local(ix_step(i):ix_step(i+1)-1,jx_step(j):jx_step(j+1)-1) = true;
            valid = I_local & edges;
            tx_valid = I_tx(valid);
            if ~isempty(tx_valid)
            foo = I_tx_edges(ix_step(i):ix_step(i+1)-1,jx_step(j):jx_step(j+1)-1);
            LOW_HIGH = stretchlim(foo, [0 0.99]);
            foo2 = imadjust(foo, LOW_HIGH, []);
            I_tx_edges(ix_step(i):ix_step(i+1)-1,jx_step(j):jx_step(j+1)-1) = foo2;
            I_tx_edges(~edges) = 0;
            end
        end
    end

    
   % I(edges & I_tx > 47000) = 1;
    I(I_tx_edges > 53000) = 1;

    filepath_chopsticks_mask_tif =  sprintf('%s%s_chopsticks_mask.tif',subpath_tiff_dir{m},m_name{m});

%     mask_chopsticks = imread(filepath_chopsticks_mask_tif);
%     
%     filepath_chopsticks2_mask_tif =  sprintf('%s%s_chopsticks2_mask.tif',subpath_tiff_dir{m},m_name{m});
%     mask_chopsticks2 = imread(filepath_chopsticks2_mask_tif);
%     I = I | mask_chopsticks | mask_chopsticks2;
%     
    %{
    filepath_felt_mask_tif =  sprintf('%s%s_felt_mask.tif',subpath_tiff_dir{m},m_name{m});
    mask_felt = imread(filepath_felt_mask_tif);
    I = I | mask_felt;
    
    filepath_uv_mask_tif =  sprintf('%s%s_uv_mask.tif',subpath_tiff_dir{m},m_name{m});
    mask_uv = imread(filepath_uv_mask_tif);
    I = I | mask_uv;
    %}
    
%     border = false(size(I));
%     border(1,:) = true;
%     border(end,:) = true;
%     border(:,1) = true;
%     border(:,end) = true;
%     
%     I(border) = false;

    % Fill edges
    locations = [1 1];%; 1 size(I,2); size(I,1) 1; size(I,1) size(I,2)];
    J = imfill(I, locations);
    
    if max(J(:)) == 0;
       didntwork = true;
       fprintf('                 \t\t%s\n', m_name{m});
        imwrite(~I, filepath_parchment_mask_tif);
        imwrite(~I, filepath_parchment_mask_matlab);
        imwrite(~I, filepath_parchment_mask_jpg);
        error('Parchment is not disconnected. Please fix in Photoshop and continue');
    end
    
    % Hole filling: Assume parchment is connected
    CC = bwconncomp(~J);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    
    for i = 1:numel(numPixels)
        if numPixels(i) < info_min_pixels
              J(CC.PixelIdxList{i}) = 1;
        end
    end

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
    %se= strel('disk',1);
    %I = imerode(I,se);
    
    fprintf('                 \t\t%s\n', m_name{m});
    imwrite(~J, filepath_parchment_mask_tif);
    imwrite(~J, filepath_parchment_mask_matlab);
    imwrite(~J, filepath_parchment_mask_jpg);
    end
end
    


%end




