% Add to cube text file 

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
    exiftoolcall = 'exiftool';
end

if ispc();
    command = sprintf(exiftoolcall);
else isunix();
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
    
[m_filepath_band] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose bands to add to cube');

m_filepath_band = m_filepath_band';
n_m = numel(m_filepath_band); 

% Find upper level directory
ix_slash = strfind(m_filepath_band{1},info_slash);
path_source = m_filepath_band{1}(1:ix_slash(end-1));

% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', path_source); 
fclose(fid);

% Determine manuscript name 
bandnames = cell(n_m,1);
for m = 1:n_m;
        filepath_full = m_filepath_band{m};
        slash_ix = strfind(filepath_full,'/');
        name = filepath_full(slash_ix(end-2)+1:slash_ix(end-1)-1);
        ix_name = strfind(filepath_full,name);
        bandname = filepath_full(ix_name(end)+numel(name):end-4);
        bandnames{m} = bandname;        
end

% cube file 
path_matlab = sprintf('%s%s+matlab%s',path_source, name,info_slash);
filepath_cube = sprintf('%s%s_cube.txt',path_matlab, name);

if exist(filepath_cube, 'file');
   % cubeno = 0;
   fid = fopen(filepath_cube, 'r');
   text = textscan(fid,'%s');
   text = text{1};
   iscube = ~cellfun(@isempty,strfind(text,'Cube'));
   ix = find(iscube);
   cubeno = str2double(text{ix(end)}(5:6))+1;
%    for l = 1:numel(iscube);
%        if ~iscube(l);
%            continue
%        end
%        cubeno_curr = str2double(text{l}(5:6)) 
%        cubeno = max(cubeno,cubeno_curr);
%    end
   fclose(fid);
else
    cubeno = 1;
end


fid = fopen(filepath_cube, 'a+');
    fprintf(fid, '\n');
fprintf(fid, 'Cube%02.0f\n', cubeno);
for m = 1:n_m
    fprintf(fid, '%s\n', bandnames{m});
end
%    fprintf(fid, '\n');



