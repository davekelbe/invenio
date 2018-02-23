function [ ] = package_data_tiffs( aux  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

%% Setup  
% Clear without breakpoints
%tmp = dbstatus;
%save('tmp.mat','tmp');
%clear all; close all; clc;
%load('tmp.mat');
%dbstop(tmp);
%clear tmp;
%delete('tmp.mat');
% Set User
info_initials = 'DJK';%inputdlg('Please type initials, e.g., DJK');
info_initials = char(info_initials);

slashix = strfind(path_target,info_slash);
%{
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
tstr = datestr(t);
tstr = strrep(tstr, ' ', '_');
tstr = strrep(tstr, ':', '-');
%}
tstr = date;

path_deliver = sprintf('%s%s_%s_%s%s',path_target(1:slashix(end-1)),'Deliver',info_initials, tstr,info_slash);
command = sprintf('mkdir %s', strrep(path_deliver, ' ', '\ '));
[~,~] = system(command);

% New stuff to rectify rotation
%path_flattened = uigetdir(path_target_previous_upper,'Please choose deliver path, e.g., Flattened_Images');
%path_flattened = sprintf('%s%s', path_flattened, info_slash);


for m = 1:n_m;
    
    path_sub = sprintf('%s%s%s', path_target, m_name{m}, info_slash);
    path_envi = sprintf('%s%s+envi%s',path_sub,m_name{m},info_slash);
    path_matlab = sprintf('%s%s+matlab%s',path_sub,m_name{m},info_slash);
    path_tiff = sprintf('%s%s+tiff%s',path_sub,m_name{m},info_slash);

    %{
    % Check rotation
    path_tif = sprintf('%s%s+tiff%s',path_sub,m_name{m},info_slash);
    filepath_I_work = sprintf('%s%s+MB365UV_001_F_stretch.tif',path_tif,m_name{m});
    I_work = imread(filepath_I_work);
    filepath_I_flat = sprintf('%s%s%s%s+MB365UV_001_F.tif',path_flattened,m_name{m},info_slash,m_name{m});
    I_flat = imread(filepath_I_flat);
    filepath_stretchval = sprintf('%s%s_stretchval.mat',path_matlab, m_name{m});
    load(filepath_stretchval);
    rotation_score = inf(2,1);
    for i = 1:2;
        rotation_score(i) = 
    filepath_I_flat = sprintf('%s%s',path_flattened,m_name{m},
    %}
    
    path_sub_deliver_tiff = sprintf('%s%s%s%s+tiff%s', path_deliver, m_name{m}, info_slash, m_name{m},info_slash );
    command = sprintf('mkdir %s', strrep(path_sub_deliver_tiff, ' ', '\ '));
    [~,~] = system(command); 
    
%     path_sub_deliver = sprintf('%s%s%s', path_deliver, m_name{m}, info_slash);
%     command = sprintf('mkdir %s', path_sub_deliver);
%     [~,~] = system(command);
%     path_envi_deliver = sprintf('%senvi%s', path_sub_deliver, info_slash);
%     command = sprintf('mkdir %s', path_envi_deliver);
%     [~,~] = system(command);
%     path_matlab_deliver = sprintf('%smatlab%s', path_sub_deliver, info_slash);
%     command = sprintf('mkdir %s', path_matlab_deliver);
%     [~,~] = system(command);
    %filepath_cubelist = sprintf('%s%s_cube.txt',path_matlab,m_name{m});
    ix_slash = strfind(path_matlab, info_slash);
    info_cube = path_matlab(ix_slash(end-2)+1:ix_slash(end-1)-1);
    
    filepath_cubelist = sprintf('%s%s_cube.txt', path_matlab,info_cube);
    if ~exist(filepath_cubelist, 'file')
        continue;
    end
    cube = parse_cubelist( filepath_cubelist );
    n_cube = numel(cube);
    all_bands = [];
    for c = 1:n_cube;
        all_bands = [all_bands; cube{c}];
    end
    all_bands = unique(all_bands);
    n_bands = numel(all_bands);
    for b = 1:n_bands;
         filepath_source = sprintf('%s%s%s.tif', path_tiff, m_name{m},all_bands{b});
         filepath_target = sprintf('%s%s%s.tif', path_sub_deliver_tiff,m_name{m},all_bands{b});
         if exist(filepath_target, 'file');
             continue;
         end
          if ispc();
             % command = sprintf('xcopy %s %s', filepath_source, filepath_target);
             %[a,b] = system(command);
             copyfile(filepath_source, filepath_target);
          else
             command = sprintf('cp %s %s',strrep(filepath_source, ' ', '\ '), strrep(filepath_target, ' ', '\ '));
             [~,~] = system(command);
          end
    end
    
%     if ispc();
%         command = sprintf('xcopy %s %s', path_envi(1:end-1), path_envi_deliver(1:end-1));
%         [~,~] = system(command);
%         command = sprintf('xcopy %s %s', path_matlab(1:end-1), path_matlab_deliver(1:end-1));
%         [~,~] = system(command);
%     else
%         command = sprintf('cp -r %s %s', path_envi,path_envi_deliver );
%         [~,~] = system(command);
%         command = sprintf('cp -r %s %s', path_matlab,path_matlab_deliver );
%         [~,~] = system(command);
%     end
end
    
% fprintf('\n Completed Successfully \n');
% fprintf('\n Please email the folder \n');
% fprintf('\n \t%s \n', path_deliver) 
% fprintf('\n to dave.kelbe@gmail.com for processing \n');

end

