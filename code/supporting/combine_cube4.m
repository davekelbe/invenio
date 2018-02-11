function [] = combine_cube4( bands )
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

%% Set User
user = 'dave';
%user = 'roger';

%% Set default paths

%current_path = cd;
%current_path = sprintf('%s/',current_path);
% Change default paths here

switch user
    case 'roger'
        slash = '\';
        rmcall = 'del';
        movetonewfolder = 0;
        default.source_path = 'g:\research\StC\';
        default.processed_dir = default.source_path;
    case 'dave'
        slash = '/';
        rmcall = 'rm';
        movetonewfolder = 1;
        default.source_path = '/Volumes/';
        default.processed_dir = default.source_path;
end
 
fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
fprintf('Save:       %s\n',default.processed_dir);

%% Update paths

options.source_dir = default.source_path;
options.processed_dir = default.processed_dir; %GUI to customize

%% Choose source file, parent files, and output directory
% Choose parent files
cd(options.source_dir);
[folders] = uipickfiles('Prompt','Please choose folders to process','REFilter', '^[^\.]', 'REDirs', true);
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

slashindex = strfind(folders{1},slash);
options.source_updir = folders{1}(1:slashindex(end));
fprintf('Path:%s \n',options.source_updir);

for i = 1:size(folders,2)
    fprintf('      %s\n',folders{i}(slashindex(end)+1:end));
end



%% Get band list 
  
n.band = size(bands,2);
n.folders = size(folders,2); 
for f = 1:n.folders;
    bandnames = cell(n.band,1);
    for i = 1:n.band
        options.source_dir = sprintf('%s/tiff/',folders{f});
        cd(options.source_dir);
        D = dir(sprintf('*_%03d_*',bands(i)));
        if numel(D) == 2;
            file = D(2).name; 
        else
            file = D.name;
        end
        for c = 1:numel(D);
            n_char = numel(D(c).name);
            if n_char==37;
                file = D(c).name;
            end
        end
        I = double(imread(file));
       % load(file);
        if i==1;
            [n.row, n.col] = size(I);
            J = zeros(n.row,n.col,n.band);
        end
        J(:,:,i) = I;
        
        bandnames{i} = file(13:23);
    end
    options.shoot_list = file(1:4);
    options.shot_seq = file(6:11);
    
    band_str = [];
    for b = 1:n.band;
        band_str = sprintf('%s_%s',band_str,sprintf('%03d',bands(b)));
    end
    bandnames = char(bandnames);
    output_filename = sprintf('%s/envi%s%s_%s_cube_norm_b%s.img',...
        folders{f},...
        slash, options.shoot_list, options.shot_seq,band_str);
    
    % Change to correct orientation
    J = imrotate(J,90);
    for b = 1:n.band;
        J(:,:,b) = flipud(J(:,:,b));
    end
    
    enviwrite_bandnames(J,output_filename, bandnames);
    fprintf('     %s\n',output_filename);
end

cd(default.source_path);


end

