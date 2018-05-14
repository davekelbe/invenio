function [ ] = gui_batch_select_ROI_wrapper2( aux )  
%UNTITLED Summary of this function goes here 
%   Detailed explanation goes here       

%% Preliminary setup
fprintf('\n***********************************************************\n');
fprintf('Choose Statistics ROIs: \n');

%m_path_upper = aux.m_path_upper;
%m_folio = aux.m_folio;
%m_mss = aux.m_mss;
%m_name = aux.m_name;
%is_band_subset = aux.is_band_subset;
%bands = aux.bands;
%info_rmcall = aux.info_rmcall;
info_slash = aux.info_slash;
%info_user = aux.info_user;
n_m = aux.n_m;
%options_delimiter = aux.options_delimiter;
%options_delimiter_wavelength = aux.options_delimiter_wavelength;
%options_folder_structure = aux.options_folder_structure;
%options_movetonewfolder = aux.options_movetonewfolder;
%path_source = aux.path_source;
%path_target = aux.path_target;
%subpath_tiff_dir = aux.path_tiff_dir;
subpath_jpg_dir = aux.subpath_jpg_dir;
%subpath_matlab_dir = aux.path_matlab_dir;
%subpath_envi_dir = aux.path_envi_dir;
%w_wavelength = aux.w_wavelength;
%w_wavelength = aux.w_wavelength;
%m_wavelength_file = aux.m_wavelength_file;
%m_wavelength_filepath = aux.m_wavelength_filepath;
%m_rotation_angle = aux.m_rotation_angle;

clear aux
%%

for m = 1:n_m;
    
cd(subpath_jpg_dir{m});

%[file.images] = uipickfiles('Prompt','Please choose images files', 'FilterSpec', '*.jpg');

counter = 1;
D = dir('*DJK_true*');
D = remove_hiddenfiles(D);
if numel(D)
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{1});
    counter = counter + 1;
end


D = dir('*KTK*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*MP*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVO22_019_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVG58_020_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end
 
D = dir('*WBUVB47_021_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVUVP_022_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVUVb_023_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end
 
D = dir('*TX940IR_031_F_stretch.jpg*'); 
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

% New lights system 

D = dir('*WBUVO22_020_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVG58_021_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end
 
D = dir('*WBUVB47_022_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVUVP_023_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

D = dir('*WBUVUVb_024_F_stretch.jpg*');
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end
 
D = dir('*TX940IR_033_F_stretch.jpg*'); 
D = remove_hiddenfiles(D);
for d = 1:numel(D);
    file.images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
    counter = counter + 1;
end

%


n_f = size(file.images,2); 
I = cell(n_f,1);
  
counter = 1;
for b = 1:n_f; 
    B = imread(file.images{b});
    B = imresize(B,0.4);
    B = (double(B)./255);
    %[~,~,n_b] = size(B);
    I{counter} = double(B); 
    counter = counter + 1;
end      

bandnames = cell(n_f);
for b = 1:n_f;
    filepath_full = file.images{b};
    slash_ix = strfind(filepath_full,info_slash);
    filename = filepath_full(slash_ix(end)+1:end-4);
    bandnames{b} = filename;
end

target_dir = filepath_full(1:slash_ix(end));
gui_select_roi2(I,target_dir,bandnames) 
end

end

