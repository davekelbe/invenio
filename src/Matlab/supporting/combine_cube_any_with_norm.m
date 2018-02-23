function [] = combine_cube_any_with_norm( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%COMBINE_CUBE2 Create a normalized ENVI image cube
%
%   There is no input to this function. Typing combine_envi_cube in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Combine ENVI Cube Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%
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
    
% %% Set User
% user = 'dave';
% %user = 'roger';
% 
% %% Set default paths
% 
% %current_path = cd;
% %current_path = sprintf('%s/',current_path);
% % Change default paths here
% 
% switch user
%     case 'roger'
%         slash = '\';
%         rmcall = 'del';
%         movetonewfolder = 0;
%         default.source_path = 'g:\research\StC\';
%         default.processed_dir = default.source_path;
%     case 'dave'
%         slash = '/';
%         rmcall = 'rm';
%         movetonewfolder = 1;
%         default.source_path = '/Volumes/';
%         default.processed_dir = default.source_path;
% end

fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
%fprintf('Source:     %s\n',default.source_path);
%fprintf('Save:       %s\n',default.processed_dir);

%% Update paths

options.source_dir = path_source_previous;

%% Choose source file, parent files, and output directory
% Choose parent files

%% Choose source file, parent files, and output directory
% Choose parent files
cd(options.source_dir);
%[SourceFile, source_path] = uigetfile('*.tif','Please choose a file for submission', 'MultiSelect', 'off');

[bands] = uipickfiles('Prompt','Please choose tif images to combine into ENVI cube', 'FilterSpec', '*.tif','REFilter', '^[^\.]', 'REDirs', true);
fprintf('\n***********************************************************\n');

% ix_slash = strfind(bands, slash);
% n.band = size(bands,2);
% for i = 1:n.band;
%    filename =  bands{i}(

% Find upper level directory
ix_slash = strfind(bands{1},info_slash);
path_source = bands{1}(1:ix_slash(end-1));

% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', path_source); 
fclose(fid);



options.processed_dir = uigetdir(path_source_previous, 'Please choose output directory');%default.processed_dir; %GUI to customize
options.processed_dir = sprintf('%s%s',options.processed_dir, info_slash);
%% Get band list
 n.band = size(bands,2);
bandnames = cell(n.band,1);
for i = 1:n.band
    file = bands{i};
    I = double(imread(file));
    I = I(:,:,1);
    % load(file);
    if i==1;
        [n.row, n.col] = size(I);
        J = zeros(n.row,n.col,n.band);
        % Normalize 
        I = I./max(I(:));
        Iview = imadjust(I,stretchlim(I), []);
        h = figure('name','Please choose spectralon for normalization');
        imagesc(Iview);
        hFH = imfreehand();
        % Create a binary image ("mask") from the ROI object.
        mask = hFH.createMask();
        delete(h);
    end
    Ispec = I(mask);
    Ispecmean = mean(Ispec);
    Inorm = I./Ispecmean;
   % Inorm(Inorm>1) = 1;
   Inorm_mean = mean(Inorm(:));
   Inorm_std = std(Inorm(:));
   maxlim = Inorm_mean+(4*Inorm_std);
   Inorm(Inorm>maxlim) = maxlim;
    J(:,:,i) = Inorm;
end


    
ix_slash = strfind(file, info_slash);
filename = file(ix_slash(end)+1:end);
ix_plus = strfind(filename, '+');
ix_KTK = strfind(filename, 'KTK');
if isempty(ix_plus);
    name = filename(1:ix_KTK-2);
else%if isempty(ix_plus);
    name = filename(1:ix_plus-1);
end

%options.processed_dir = sprintf('%s%s+envi%s',file(1:ix_slash(end-1)),name,slash);
%options.processed_dir = '/Volumes/Corinth/MOTB/Cubes/';
file_filename = file(ix_slash(end)+1:end);
ix_underscore = strfind(file_filename, '_');
ix_plus = strfind(file_filename, '+');
options.shoot_list = file_filename(1:ix_underscore(1)-1);
options.shot_seq = file_filename(ix_underscore(1)+1:ix_plus-1);
if isempty(options.shot_seq);
    options.shot_seq = file_filename(ix_underscore(1)+1:ix_underscore(2)-1);
end

counter = 1;
output_filename = sprintf('%s%s_%s_cube_norm_%02.0f.img',options.processed_dir,...
    options.shoot_list, options.shot_seq,counter);
while exist(output_filename);
    counter = counter + 1;
    output_filename = sprintf('%s%s_%s_cube_norm_%02.0f.img',options.processed_dir,...
    options.shoot_list, options.shot_seq,counter);
end

bands2 = bands;
for b = 1:n.band
    bands2{b} = bands{b}(ix_slash(end)+1:end);
end

text_filename = sprintf('%stxt',output_filename(1:end-3));
fid = fopen(text_filename, 'w+');
for b  = 1:n.band;
fprintf(fid, sprintf('\n%s\n',bands2{b}));
end


bandnames = char(bands2);

% Change to correct orientation
J = imrotate(J,90);
for b = 1:n.band;
    J(:,:,b) = flipud(J(:,:,b));
end

enviwrite_bandnames(J,output_filename, bandnames);
fprintf('     %s\n',output_filename);

%cd(default.source_path);


end

