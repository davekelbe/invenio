function [  ] = resize_KTK_jpgs(  )
%MAKE_STRETCHED_FLATS Saves stretched single band images of flattened
%image for display
%
%   There is no input to this function. Typing MAKE_STRETCHED_FLATS in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Make Stretched Flats Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%
% V0.0 - Initial Version - February 23 2013
%
%
% Requirements:
%
%
% Tips:
%   * Press ctrl+c to cancel execution and restart
%   *Set default paths in source code for efficiency
fprintf('\n***********************************************************\n');
fprintf('Tips\n');
fprintf('            Press ctrl+c to cancel execution and restart\n');

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
% m_name = cell(n_m,1);
% for m = 1:n_m;
%     ix_slash = strfind(m_path_upper{m},info_slash);
%     m_name{m} = m_path_upper{m}(ix_slash(end-1)+1:end-1);
% end

for m = 1:n_m;
    D = dir(m_path_upper{m});
    D = remove_hiddenfiles(D);
    n_D = numel(D);
    for d = 1:n_D;
        dir_tiff = sprintf('%s%s%s%s+temp%s',m_path_upper{m},D{d},info_slash,D{d},info_slash);
        F = dir(dir_tiff);
        F = remove_hiddenfiles(F);
        n_F = numel(F);
        for f = 1:n_F;
        filepath_F_jpg = sprintf('%s%s%s%s+jpg%s%sjpg',m_path_upper{m},D{d},info_slash,D{d},info_slash,F{f}(1:end-3));
        filepath_F_tiff = sprintf('%s%s%s%s+temp%s%s',m_path_upper{m},D{d},info_slash,D{d},info_slash,F{f});

        if ~exist(filepath_F_jpg, 'file');
            I = imread(filepath_F_tiff);
            J = imresize(I,.4);
            iminfo = imfinfo(filepath_F_tiff);
            bitdepth = iminfo.BitsPerSample(1);
            maxval = 2^bitdepth;
            J = uint8(256*double(J)/maxval);
            imwrite(J,filepath_F_jpg,'jpeg','Quality', 50);
        end
        end
    end
end
    
    
