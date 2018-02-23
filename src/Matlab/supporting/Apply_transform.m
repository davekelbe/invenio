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
isKS = false; 
info_user = 'generic';
switch info_user
    case 'emel'
        options_delimiter = '_';
        options_delimiter_wavelength = '+';
        options_folder_structure  = 'mss_folio';
        options_movetonewfolder = 1;
    case 'generic'
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
    if isKS;
        exiftoolcall = 'C:\Users\KevinSacca\Documents\MATLAB\exiftool.pl';
    end
else isunix();
    info_slash = '/';
    info_root = '/';
    info_rmcall = 'rm';
    exiftoolcall = '/usr/bin/exiftool';
end

% if ispc();
%     command = sprintf(exiftoolcall);
% else isunix();
%     command = sprintf('which %s', exiftoolcall);
% end
% 
% [exiftf,~] = system(command);
% if exiftf;
%     error('Please install Exiftool');
% end

% Add Matlab directory to path (?)
% addpath

% OS-independent root directory 
filepath_matlab = matlabroot;
ix_slash = strfind(filepath_matlab,info_slash);
path_matlab = filepath_matlab(1:ix_slash(end));
% May have to define new root directory if no write permissions, e.g.,
if isKS;
    path_matlab = 'C:\Users\KevinSacca\Documents\';
end
% Determine previous directory for source data 
filepath_source_previous = sprintf('%spath_source_previous.txt', ...
    path_matlab);
if exist(filepath_source_previous, 'file');
    fid = fopen(filepath_source_previous);
    path_source_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_source_previous = char(path_source_previous{1});
else
    path_source_previous = info_root; 
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_source_previous, 'dir');
    path_source_previous = info_root;
end
    
% Choose parent files
%cd();
[m_path_upper] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose folders to process');
m_path_upper = m_path_upper';
n_m = numel(m_path_upper); 
for m = 1:n_m;
    m_path_upper{m} = sprintf('%s%s',m_path_upper{m},info_slash);
end
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

% Find upper level directory
ix_slash = strfind(m_path_upper{1},info_slash);
path_source = m_path_upper{1}(1:ix_slash(end-1));


% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', path_source); 
fclose(fid);

% Determine manuscript name 
m_name = cell(n_m,1);
for m = 1:n_m;
    ix_slash = strfind(m_path_upper{m},info_slash);
    m_name{m} = m_path_upper{m}(ix_slash(end-1)+1:end-1);
end

% Determine mss and folio
m_mss = cell(n_m,1);
m_folio = cell(n_m,1);
switch options_folder_structure
    case 'mss_folio';
        for m = 1:n_m;
        ix_delimiter = strfind(m_name{m}, options_delimiter);
        m_mss{m} = m_name{m}(1:ix_delimiter(end)-1);
        m_folio{m} = m_name{m}(ix_delimiter(end)+1:end);
        end
end

% Determine transformation directory
filepath_trans_previous = sprintf('%spath_trans_previous.txt', ...
    path_matlab);
if exist(filepath_trans_previous, 'file');
    fid = fopen(filepath_trans_previous);
    path_trans_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_trans_previous = char(path_trans_previous{1});
    ix_slash = strfind(path_trans_previous, info_slash);
    path_target_previous_upper = path_trans_previous(1:ix_slash(end-1));
else
    path_trans_previous = info_root;
    path_target_previous_upper = info_root;
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_trans_previous, 'dir');
    path_trans_previous = info_root;
end

% Determine source transformations 
[m_trans] = uipickfiles('filterspec',path_trans_previous,'Prompt','Please choose transformations to replicate');
m_trans = m_trans';
n_m = numel(m_trans); 
path_trans = m_trans{1}(1:ix_slash(end));
ix_slash = strfind(m_trans{1}, info_slash);
for m = 1:n_m;
    m_trans{m} = sprintf('%s%s',m_trans{m}(ix_slash(end)+1:end));
end

path_trans = char(path_trans);
if ~strcmp(path_trans(end), info_slash);
    path_trans = sprintf('%s%s',path_trans, info_slash);
end

% Update target directory 
fid = fopen(filepath_trans_previous, 'w+');
fprintf(fid, '%s', path_trans); 
fclose(fid);
if ~exist(path_trans, 'dir')
    mkdir(path_trans);
end




subpath_tiff_dir = cell(n_m,1);
subpath_jpg_dir = cell(n_m,1);
subpath_matlab_dir = cell(n_m,1);
subpath_envi_dir = cell(n_m,1);
% Make additional directories 
for m = 1:n_m;
    subpath_tiff_dir{m} = sprintf('%s%s_%s%s%s_%s+tiff%s',...
        path_trans,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_jpg_dir{m} = sprintf('%s%s_%s%s%s_%s+jpg%s',...
        path_trans,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_matlab_dir{m} = sprintf('%s%s_%s%s%s_%s+matlab%s',...
        path_trans,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_envi_dir{m} = sprintf('%s%s_%s%s%s_%s+envi%s',...
        path_trans,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    if ~exist(subpath_tiff_dir{m},'dir');
        mkdir(subpath_tiff_dir{m})
    end
    if ~exist(subpath_jpg_dir{m},'dir');
        mkdir(subpath_jpg_dir{m})
    end
    if ~exist(subpath_envi_dir{m},'dir');
        mkdir(subpath_envi_dir{m})
    end
    if ~exist(subpath_matlab_dir{m},'dir');
        mkdir(subpath_matlab_dir{m})
    end
end

% Print outputs 
fprintf('\n');
fprintf('Source Directory:\t\t%s\n', path_source);
fprintf('Files to process:\n');
for m = 1:n_m;
fprintf('                 \t\t%s\n', m_name{m});
end


clear fid filepath_matlab filepath_source_previous path_matlab slashindex
clear path_source_previous m ans filepath_target_previous ix_slash 
clear ix_delimiter k path_target_previous_upper path_target_previous
% Output 
% path_source              - Path to source (flattened) data
% m_path_upper                - Concatenated path and name 
% m_name                   - Cell array of manuscript names 
% info                     - predefined variables   
%%



% Choose images to replicate 

% Choose folios to apply transform to 