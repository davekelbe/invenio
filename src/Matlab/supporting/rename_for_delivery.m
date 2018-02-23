% Rename files 

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
    
path_source = uigetdir(path_source_previous,'Please choose finished path, e.g., Finished');
D = dir(path_source);
D = remove_hiddenfiles(D);
if ~strcmp(path_source(end), info_slash);
    path_source = sprintf('%s%s', path_source, info_slash);
end

ix_slash = strfind(path_source, info_slash);
path_up = path_source(1:ix_slash(end-1));
path_target = '/Users/Kelbe/Desktop/Kelbe-NA/';
%path_target = sprintf('%s%s%s',path_up, 'ShortFnms',info_slash);
if ~exist(path_target, 'dir');
    mkdir(path_target);
end

n_D = numel(D);
for d = 1:n_D;
    increment = true;
    start = 1;
    ix_delim = strfind(D{d},'_');
    name = D{d}(1:ix_delim(2)-1);
    while increment       
        filepath_cand = sprintf('%s%s+DJK_%02.0f.tif',path_target,name,start);
        if exist(filepath_cand,'file');
            start = start + 1;
        else
        increment = false;
        end
    end
    filepath_I = sprintf('%s%s',path_source, D{d});
    I = imread(filepath_I);
    info = imfinfo(filepath_I);
    BitsPerSample = info.BitsPerSample(1);
    if BitsPerSample ~= 8;
        I = double(I);
        I = uint8(I./max(I(:))*256);
    end 
    imwrite(I,filepath_cand,'tif','Compression', 'none')
    filepath_jpg = sprintf('%sjpg',filepath_cand(1:end-3));
    imwrite(I,filepath_jpg,'jpeg','Quality', 100)
end
foo = 1;




