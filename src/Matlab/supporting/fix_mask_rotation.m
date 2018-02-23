% Fix masks to make sure they match the flattened data 
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
%info_initials = inputdlg('Please type initials, e.g., DJK');
%info_initials = char(info_initials);
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
    exiftoolcall = '/usr/bin/exiftool';
end

% if ispc();
%     command = sprintf('%s',exiftoolcall);
% else 
%     command = sprintf('which %s', exiftoolcall);
% end
% [exiftf,~] = system(command);
% if exiftf;
%     error('Please install Exiftool');
% end

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

path_delivery = uigetdir(path_target_previous_upper,'Please choose deliver path, e.g., Deliver_DJK_05-Oct-2015');
path_flattened = uigetdir(path_target_previous_upper,'Please choose deliver path, e.g., Flattened_Images');

D = dir(path_target);
D = remove_hiddenfiles(D);

n_D = numel(D);
for d = 1:n_D;
    path_d = sprintf('



