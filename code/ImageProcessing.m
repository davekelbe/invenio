function [  ] = ImageProcessing(  )
%IMAGEPROCESSING Wrapper for image processing steps
%
%   There is no input to this function. Typing ImageProcessing in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Image Processing   Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%
% V0.0 - Initial Version - February 6 2015
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
cd(path_source_previous);
[m_path_upper] = uipickfiles('Prompt','Please choose folders to process');
m_path_upper = m_path_upper';
n_m = numel(m_path_upper); 
for m = 1:n_m;
    m_path_upper{m} = sprintf('%s%s',m_path_upper{m},info_slash);
  %  m_path_upper{m} = strrep(m_path_upper{m}, ' ', '\ ');
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

% Determine target directory
filepath_target_previous = sprintf('%spath_target_previous.txt', ...
    path_matlab);
if exist(filepath_target_previous, 'file');
    fid = fopen(filepath_target_previous);
    path_target_previous = textscan(fid, '%s', 'delimiter', '\t');
    path_target_previous = char(path_target_previous{1});
    ix_slash = strfind(path_target_previous, info_slash);
    path_target_previous_upper = path_target_previous(1:ix_slash(end-1));
else
    path_target_previous = info_root;
    path_target_previous_upper = info_root;
end
% Change source directory if no longer exists (e.g. drive removed) 
if ~exist(path_target_previous, 'dir');
    path_target_previous = info_root;
end
path_target = uigetdir(path_target_previous_upper,'Please choose output path');
path_target = char(path_target);
if ~strcmp(path_target(end), info_slash);
    path_target = sprintf('%s%s',path_target, info_slash);
end
%path_target = strrep(path_target, ' ', '\ ');
k = strfind(path_target, 'Processed');
if isempty(k);
    %path_target = sprintf('%sProcessed-%s%s', path_target, m_mss{1}, info_slash);
    path_target = sprintf('%sProcessed%s', path_target, info_slash);
end

% Update target directory 
fid = fopen(filepath_target_previous, 'w+');
fprintf(fid, '%s', path_target); 
fclose(fid);
%path_target_nospace = strrep(path_target, '\ ', ' ');
if ~exist(path_target, 'dir')
    mkdir(path_target);
end

subpath_tiff_dir = cell(n_m,1);
subpath_jpg_dir = cell(n_m,1);
subpath_tiff_mask_dir = cell(n_m,1);
subpath_jpg_mask_dir = cell(n_m,1);
subpath_matlab_dir = cell(n_m,1);
subpath_envi_dir = cell(n_m,1);
% Make additional directories 
for m = 1:n_m;
    subpath_tiff_dir{m} = sprintf('%s%s_%s%s%s_%s+tiff%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_jpg_dir{m} = sprintf('%s%s_%s%s%s_%s+jpg%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_tiff_mask_dir{m} = sprintf('%s%s_%s%s%s_%s+tiffm%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_jpg_mask_dir{m} = sprintf('%s%s_%s%s%s_%s+jpgm%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_matlab_dir{m} = sprintf('%s%s_%s%s%s_%s+matlab%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_envi_dir{m} = sprintf('%s%s_%s%s%s_%s+envi%s',...
        path_target,...
        m_mss{m},m_folio{m},info_slash,m_mss{m},m_folio{m},info_slash);
    subpath_tiff_dir_nospace = strrep(subpath_tiff_dir{m}, '\ ', ' ');
    if ~exist(subpath_tiff_dir_nospace,'dir');
        mkdir(subpath_tiff_dir{m})
    end
    subpath_jpg_dir_nospace = strrep(subpath_jpg_dir{m}, '\ ', ' ');
    if ~exist(subpath_jpg_dir_nospace,'dir');
        mkdir(subpath_jpg_dir{m})
    end
    subpath_tiff_mask_dir_nospace = strrep(subpath_tiff_mask_dir{m}, '\ ', ' ');
    if ~exist(subpath_tiff_mask_dir_nospace,'dir');
        mkdir(subpath_tiff_mask_dir{m})
    end
    subpath_jpg_mask_dir_nospace = strrep(subpath_jpg_mask_dir{m}, '\ ', ' ');
    if ~exist(subpath_jpg_mask_dir_nospace,'dir');
        mkdir(subpath_jpg_mask_dir{m})
    end
    subpath_envi_dir_nospace = strrep(subpath_envi_dir{m}, '\ ', ' ');
    if ~exist(subpath_envi_dir_nospace,'dir');
        mkdir(subpath_envi_dir{m})
    end
    subpath_matlab_dir_nospace = strrep(subpath_matlab_dir{m}, '\ ', ' ');
    if ~exist(subpath_matlab_dir_nospace,'dir');
        mkdir(subpath_matlab_dir{m})
    end
end

% Print outputs 
fprintf('\n');
fprintf('Source Directory:\t\t%s\n', path_source);
fprintf('Target Directory:\t\t%s\n', path_target);
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

%% Set up aux 
aux.m_path_upper = m_path_upper;
aux.m_folio = m_folio;
aux.m_mss = m_mss;
aux.m_name = m_name; 
aux.is_band_subset = true;
aux.bands = true;
aux.info_rmcall = info_rmcall;
aux.info_slash = info_slash; 
aux.info_user  = info_user;
aux.info_colormap = 'jet';
if isKS; 
    aux.info_colormap = 'jet';
end
aux.n_m = n_m;
aux.options_delimiter = options_delimiter;
aux.options_delimiter_wavelength = options_delimiter_wavelength;
aux.options_folder_structure = options_folder_structure;
aux.options_movetonewfolder = options_movetonewfolder;
aux.path_source = path_source; 
aux.path_target = path_target; 
aux.path_tiff_dir = subpath_tiff_dir;
aux.path_jpg_dir = subpath_jpg_dir;
aux.path_tiff_mask_dir = subpath_tiff_mask_dir;
aux.path_jpg_mask_dir = subpath_jpg_mask_dir;
aux.path_matlab_dir = subpath_matlab_dir;
aux.path_envi_dir = subpath_envi_dir;
aux.exiftoolcall = exiftoolcall;

%% Setup aux
%aux.m_wavelength_file_new = m_wavelength_file_new;

aux = setup_aux(aux);

%% Create spectralon mask 
%aux.m_wavelength_file_new = m_wavelength_file_new;

%create_spectralon_mask(aux);

%% Create chopstick mask 
%create_chopsticks_mask(aux);

%create_chopsticks2_mask(aux);
%create_uv_mask(aux);

%create_felt_mask(aux);


%% Create parchment mask 
%create_parchment_mask(aux);
%% Create reflectance tiffs 

%reflectance_tiffs11(aux);

%% Create ink mask 
create_sharpie_mask(aux);
%create_overtext_mask(aux);
%create_undertext_mask(aux);

%% Create truecolor images 
%aux.m_rotation_angle = 0;%m_rotation_angle; 
%reflectance_tiffs9_rgb(aux);



 
%filepath_rotation = sprintf('%srotation.mat', subpath_matlab_dir);
%load(filepath_rotation);

%% Register verso 

%register_verso_flipud(aux);




%% Resize auxiliary files 

%resize_auxiliary_files(aux)

%% Select ROI for statistics 
%gui_batch_select_ROI_wrapper2(aux)
 

%% Select band subsets 

%choose_cube_bands_aux(aux);

%% Package data for IDL/ENVI Processing 
package_data(aux);
package_data_tiffs(aux);

end

