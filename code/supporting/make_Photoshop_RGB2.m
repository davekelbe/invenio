function [  ] = make_Photoshop_RGB2(  )
%MAKE_PHOTOSHOP_RGB Creates a 2% Contrast Adjusted, 16-bit Photoshop RGB
%TIF Image for fine-tuning contrast, hue, etc.
%
% 
%   There is no input to this function. Typing make_Photoshop_RGB in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Make Photoshop RGB Tool
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
% Tips:
%   * Press ctrl+c to cancel execution and restart
%   *Set default paths in source code for efficiency
fprintf('\n***********************************************************\n');
fprintf('Tips\n');
fprintf('            Press ctrl+c to cancel execution and restart\n');
fprintf('            *Change default paths in source code (line 57-60)\n');
%%
%% Setup  
isGR = false;
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
%filepath_matlab = strrep(filepath_matlab, ' ' , '\ ');
ix_slash = strfind(filepath_matlab,info_slash);
path_matlab = filepath_matlab(1:ix_slash(end));
% May have to define new root directory if no write permissions, e.g.,
if isKS;
    path_matlab = 'C:\Users\KevinSacca\Documents\';
end
if isGR;
    path_matlab = 'F:\';
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
     path_source = path_source_previous;


% 2% Contrast Stretch
lowfrac = .0002;
TOL = [lowfrac 1-lowfrac];
fprintf('\n%g%% Contrast Stretch\n', 100*(lowfrac*2));

%% Choose source file, parent files, and output directory
% Choose parent files
cd(path_source);
%[SourceFile, source_path] = uigetfile('*.tif','Please choose a file for submission', 'MultiSelect', 'off');

[files] = uipickfiles('Prompt','Please choose 2-3 images to combine into Photoshop RGB', 'FilterSpec', '*.tif','REFilter', '^[^\.]', 'REDirs', true);
fprintf('\n***********************************************************\n');
ix_slash = strfind(files{1},info_slash);
path_source = files{1}(1:ix_slash(end-2));
% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', path_source); 
fclose(fid);

n.files = size(files,2);

%info = imfinfo(files{1});
%n.r = info.Height; 
%n.c = info.Width;
%n.c = 7216;
%n.r = 5412;
%I = uint16(zeros(n.r,n.c,3));

if n.files > 3;
    error('Please supply 2 or 3 images');
end

bandname = cell(3,1);
for f = 1:n.files
        A = double(imread(files{f}));
        if size(A,3)~=1;
            A = A(:,:,1);
        end
        if f ==1; 
            [n.r,n.c] = size(A);
            I = uint16(zeros(n.r,n.c,3));
        end
        A = A./max(A(:));
        A = imadjust(A,stretchlim(A,TOL),[]);
        A = double(A);
        A = A-min(A(:));
        A = A./max(A(:));
        A = uint16(65536*A);
        if f == 2 && n.files ==2;
        I(:,:,2) = A;
        I(:,:,3) = A;
        bandname{2} = files{f}(end-16:end-10);
        bandname{3} = files{f}(end-16:end-10);
        else
        I(:,:,f) = A;
        bandname{f} = files{f}(end-16:end-10);
        end
end
f=1;
outfilename = sprintf('%s',files{1}(1:end-4));
k = strfind(outfilename,info_slash);
%outfilename = outfilename(1:k(end-1));
%outfilename = sprintf('%sRGB%s%s',outfilename,info_slash,temp);
for f = 2:n.files;
    temp_outfilename = sprintf('%s',files{f}(1:end-4));
    k = strfind(temp_outfilename,info_slash);
    temp = temp_outfilename(k(end)+13:end);
    outfilename = sprintf('%s+%s',outfilename,temp);
end
outfilename = sprintf('%s_DJK_RGB.tif',outfilename);

imwrite(I,outfilename,'tif')


fprintf('\n***********************************************************\n');
fprintf('\nCompleted successfully\n');


end









