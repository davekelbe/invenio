function [ ] = package_data_tiffs(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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
info_initials = inputdlg('Please type initials, e.g., DJK');
info_initials = char(info_initials);
%info_initials = 'GR';
info_user = 'emel';
switch info_user
    case 'emel'
        options_delimiter = '_';
        options_delimiter_wavelength = '+';
        options_folder_structure  = 'mss_folio';
        options_movetonewfolder = 1;
end

% Define OS-specific parameters 
if ispc();
    info_slash = '\';
    info_root = 'C:\'; 
    info_rmcall = 'del';
    exiftoolcall = 'exiftool.pl';
    % May have to define full path if not in $PATH, e.g., below 
    % exiftoolcall = 'C:\Users\KevinSacca\Documents\MATLAB\exiftool.pl';
else isunix();
    info_slash = '/';
    info_rmcall = 'rm';
    exiftoolcall = 'exiftool';
end

if ispc();
    command = sprintf('%s',exiftoolcall);
else 
    command = sprintf('which %s', exiftoolcall);
end
[exiftf,~] = system(command);
if exiftf;
    error('Please install Exiftool');
end

% Add Matlab directory to path (?)
% addpath

% OS-independent root directory 
filepath_matlab = matlabroot;
ix_slash = strfind(filepath_matlab,info_slash);
path_matlab = filepath_matlab(1:ix_slash(end));

filepath_target_previous = sprintf('%spath_target_previous.txt', ...
    path_matlab);
if exist(filepath_target_previous, 'file');
    fid = fopen(filepath_target_previous);
    path_target_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_target_previous = char(path_target_previous{1});
    ix_slash = strfind(path_target_previous, info_slash);
    path_target_previous_upper = path_target_previous(1:ix_slash(end-1));
else
    path_target_previous = info_slash;
    path_target_previous_upper = info_slash;
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_target_previous, 'dir');
    path_target_previous = info_slash;
end
path_target = uigetdir(path_target_previous_upper,'Please choose working path, e.g., Processed-0065');
path_target = char(path_target);
if ~strcmp(path_target(end), info_slash);
    path_target = sprintf('%s%s',path_target, info_slash);
end

slashix = strfind(path_target,info_slash);
%{
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
tstr = datestr(t);
tstr = strrep(tstr, ' ', '_');
tstr = strrep(tstr, ':', '-');
%}
tstr = date;

path_deliver = sprintf('%s%s_%s_%s%s',path_target(1:slashix(end-1)),'Deliver',info_initials, tstr,info_slash);
path_deliver_tiff = sprintf('%s%s_%s_%s%s',path_target(1:slashix(end-1)),'Deliver-tiff',info_initials, tstr,info_slash);
command = sprintf('mkdir %s', path_deliver_tiff);
[~,~] = system(command);
command = sprintf('mkdir %s', path_deliver);
[~,~] = system(command);

D = dir(path_target);
D = remove_hiddenfiles(D);
n_D = numel(D);

% New stuff to rectify rotation
%path_flattened = uigetdir(path_target_previous_upper,'Please choose deliver path, e.g., Flattened_Images');
%path_flattened = sprintf('%s%s', path_flattened, info_slash);


for d = 1:n_D;
    
    path_sub = sprintf('%s%s%s', path_target, D{d}, info_slash);
    path_envi = sprintf('%s%s+envi%s',path_sub,D{d},info_slash);
    path_matlab = sprintf('%s%s+matlab%s',path_sub,D{d},info_slash);
    path_tiff = sprintf('%s%s+tiff%s',path_sub,D{d},info_slash);

    %{
    % Check rotation
    path_tif = sprintf('%s%s+tiff%s',path_sub,D{d},info_slash);
    filepath_I_work = sprintf('%s%s+MB365UV_001_F_stretch.tif',path_tif,D{d});
    I_work = imread(filepath_I_work);
    filepath_I_flat = sprintf('%s%s%s%s+MB365UV_001_F.tif',path_flattened,D{d},info_slash,D{d});
    I_flat = imread(filepath_I_flat);
    filepath_stretchval = sprintf('%s%s_stretchval.mat',path_matlab, D{d});
    load(filepath_stretchval);
    rotation_score = inf(2,1);
    for i = 1:2;
        rotation_score(i) = 
    filepath_I_flat = sprintf('%s%s',path_flattened,D{d},
    %}
    
    path_sub_deliver_tiff = sprintf('%s%s%s', path_deliver_tiff, D{d}, info_slash);
    command = sprintf('mkdir %s', path_sub_deliver_tiff);
    [~,~] = system(command); 
    
%     path_sub_deliver = sprintf('%s%s%s', path_deliver, D{d}, info_slash);
%     command = sprintf('mkdir %s', path_sub_deliver);
%     [~,~] = system(command);
%     path_envi_deliver = sprintf('%senvi%s', path_sub_deliver, info_slash);
%     command = sprintf('mkdir %s', path_envi_deliver);
%     [~,~] = system(command);
%     path_matlab_deliver = sprintf('%smatlab%s', path_sub_deliver, info_slash);
%     command = sprintf('mkdir %s', path_matlab_deliver);
%     [~,~] = system(command);
    filepath_cubelist = sprintf('%s%s_cube.txt',path_matlab,D{d});
    cube = parse_cubelist( filepath_cubelist );
    n_cube = numel(cube);
    all_bands = [];
    for c = 1:n_cube;
        all_bands = [all_bands; cube{c}];
    end
    all_bands = unique(all_bands);
    n_bands = numel(all_bands);
    for b = 1:n_bands;
         filepath_source = sprintf('%s%s%s.tif', path_tiff, D{d},all_bands{b});
         filepath_target = sprintf('%s%s%s.tif', path_sub_deliver_tiff, D{d},all_bands{b});
          if ispc();
              command = sprintf('xcopy %s %s', filepath_source, filepath_target);
             [~,~] = system(command);
          else
             command = sprintf('cp %s %s',filepath_source, filepath_target);
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
    
fprintf('\n Completed Successfully \n');
fprintf('\n Please email the folder \n');
fprintf('\n \t%s \n', path_deliver) 
fprintf('\n to dave.kelbe@gmail.com for processing \n');

end

