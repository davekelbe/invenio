function [ aux ] = processing_setup(aux )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Set User
switch aux.info_user
    case 'emel'
        aux.options_delimiter = '_';
        aux.options_delimiter_wavelength = '+';
        aux.options_folder_structure  = 'mss_folio';
        aux.options_movetonewfolder = 1;
    case 'generic'
        aux.options_delimiter = '_';
        aux.options_delimiter_wavelength = '+';
        aux.options_folder_structure  = 'mss_folio';
        aux.options_movetonewfolder = 1;    
end

% Define OS-specific parameters 
if ispc()
    aux.info_slash = '\';
    aux.info_root = 'C:\';
    aux.info_rmcall = 'del';
    %aux.exiftoolcall = 'exiftool.pl'; %Roger
    aux.exiftoolcall = 'exiftool'; %Damian 
    % May have to define full path if not in $PATH, e.g., below
elseif isunix()
    aux.info_slash = '/';
    aux.info_root = '/';
    aux.info_rmcall = 'rm';
    aux.exiftoolcall = '/usr/local/bin/exiftool';
end

% Set up exiftool 
if ispc()
    command = sprintf(aux.exiftoolcall);
elseif isunix()
    command = sprintf('which %s', aux.exiftoolcall);
end
[exiftf,~] = system(command);
if exiftf
%    error('Please install Exiftool'); % Cambridge - not sure why not
%    working? 
end

% Add Matlab directory to path (?)
% addpath

% OS-independent root directory 
filepath_matlab = matlabroot;
ix_slash = strfind(filepath_matlab,aux.info_slash);
path_matlab = filepath_matlab(1:ix_slash(end));
% May have to define new root directory if no write permissions, e.g.,
%    path_matlab = 'C:\Users\KevinSacca\Documents\';
% Determine previous directory for source data 
filepath_source_previous = sprintf('%spath_source_previous.txt', ...
    path_matlab);
if exist(filepath_source_previous, 'file')
    fid = fopen(filepath_source_previous);
    path_source_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_source_previous = char(path_source_previous{1});
else
    path_source_previous = aux.info_root; 
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_source_previous, 'dir')
    path_source_previous = aux.info_root;
end
    
% Choose parent files
%cd();
[aux.m_path_upper] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose folders to process');
aux.m_path_upper = aux.m_path_upper';
aux.n_m = numel(aux.m_path_upper); 
for m = 1:aux.n_m
    aux.m_path_upper{m} = sprintf('%s%s',aux.m_path_upper{m},aux.info_slash);
end
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

% Find upper level directory
ix_slash = strfind(aux.m_path_upper{1},aux.info_slash);
aux.path_source = aux.m_path_upper{1}(1:ix_slash(end-1));

% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', aux.path_source); 
fclose(fid);

% Determine manuscript name 
aux.m_name = cell(aux.n_m,1);
for m = 1:aux.n_m
    ix_slash = strfind(aux.m_path_upper{m},aux.info_slash);
    aux.m_name{m} = aux.m_path_upper{m}(ix_slash(end-1)+1:end-1);
end

% Determine mss and folio
aux.m_mss = cell(aux.n_m,1);
aux.m_folio = cell(aux.n_m,1);
switch aux.options_folder_structure
    case 'mss_folio'
        for m = 1:aux.n_m
        ix_delimiter = strfind(aux.m_name{m}, aux.options_delimiter);
        aux.m_mss{m} = aux.m_name{m}(1:ix_delimiter(end)-1);
        aux.m_folio{m} = aux.m_name{m}(ix_delimiter(end)+1:end);
        end
end

% Determine target directory
filepath_target_previous = sprintf('%spath_target_previous.txt', ...
    path_matlab);
if exist(filepath_target_previous, 'file')
    fid = fopen(filepath_target_previous);
    path_target_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_target_previous = char(path_target_previous{1});
    ix_slash = strfind(path_target_previous, aux.info_slash);
    path_target_previous_upper = path_target_previous(1:ix_slash(end-1));
else
    path_target_previous = aux.info_root;
    path_target_previous_upper = aux.info_root;
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_target_previous, 'dir')
    path_target_previous = aux.info_root;
end
aux.path_target = uigetdir(path_target_previous_upper,'Please choose output path');
aux.path_target = char(aux.path_target);
if ~strcmp(aux.path_target(end), aux.info_slash)
    aux.path_target = sprintf('%s%s',aux.path_target, aux.info_slash);
end
k = strfind(aux.path_target, 'Processed');
if isempty(k)
    aux.path_target = sprintf('%sProcessed-%s%s', aux.path_target, aux.m_mss{1}, aux.info_slash);
end

% Update target directory 
aux.path_target = '/Volumes/Tyndale/Tyndale/Processed-CCR/'; % Hack Cambridge
fid = fopen(filepath_target_previous, 'w+');
fprintf(fid, '%s', aux.path_target); 
fclose(fid);
if ~exist(aux.path_target, 'dir')
    mkdir(aux.path_target);
end

aux.subpath_tiff_dir = cell(aux.n_m,1);
aux.subpath_jpg_dir = cell(aux.n_m,1);
aux.subpath_matlab_dir = cell(aux.n_m,1);
aux.subpath_envi_dir = cell(aux.n_m,1);
% Make additional directories 
for m = 1:aux.n_m
    aux.subpath_tiff_dir{m} = sprintf('%s%s_%s%s%s_%s+tiff%s',...
        aux.path_target,...
        aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
    aux.subpath_jpg_dir{m} = sprintf('%s%s_%s%s%s_%s+jpg%s',...
        aux.path_target,...
        aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
    aux.subpath_matlab_dir{m} = sprintf('%s%s_%s%s%s_%s+matlab%s',...
        aux.path_target,...
        aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
    aux.subpath_envi_dir{m} = sprintf('%s%s_%s%s%s_%s+envi%s',...
        aux.path_target,...
        aux.m_mss{m},aux.m_folio{m},aux.info_slash,aux.m_mss{m},aux.m_folio{m},aux.info_slash);
    if ~exist(aux.subpath_tiff_dir{m},'dir')
        mkdir(aux.subpath_tiff_dir{m})
    end
    if ~exist(aux.subpath_jpg_dir{m},'dir')
        mkdir(aux.subpath_jpg_dir{m})
    end
    if ~exist(aux.subpath_envi_dir{m},'dir')
        mkdir(aux.subpath_envi_dir{m})
    end
    if ~exist(aux.subpath_matlab_dir{m},'dir')
        mkdir(aux.subpath_matlab_dir{m})
    end
end

% Print outputs 
fprintf('\n');
fprintf('Source Directory:\t\t%s\n', aux.path_source);
fprintf('Target Directory:\t\t%s\n', aux.path_target);
fprintf('Files to process:\n');
for m = 1:aux.n_m
fprintf('                 \t\t%s\n', aux.m_name{m});
end

aux.is_band_subset = true;
aux.bands = true;
aux.info_colormap = 'parula';

clear fid filepath_matlab filepath_source_previous path_matlab slashindex
clear path_source_previous m ans filepath_target_previous ix_slash 
clear ix_delimiter k path_target_previous_upper path_target_previous
clear command exiftf exiftoolcall 


% Output 
% path_source              - Path to source (flattened) data
% aux.m_path_upper                - Concatenated path and name 
% m_name                   - Cell array of manuscript names 
% info                     - predefined variables   

end

