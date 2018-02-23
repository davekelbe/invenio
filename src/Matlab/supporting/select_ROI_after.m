fprintf('\n***********************************************************\n');
fprintf('Tips\n');
fprintf('            Press ctrl+c to cancel execution and restart\n');
fprintf('            *Change default paths in source code (line 57-60)\n');
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
[m_path_upper] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose [processed] folders to process');
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

subpath_jpg_dir = cell(n_m,1);
for m = 1:n_m;
  subpath_jpg_dir{m} = sprintf('%s%s%s_%s+jpg%s',...
        m_path_upper{m},...
        info_slash,m_mss{m},m_folio{m},info_slash);
end

aux.n_m = n_m;
aux.path_jpg_dir = subpath_jpg_dir;
aux.info_slash = info_slash;

gui_batch_select_ROI_wrapper2(aux)
