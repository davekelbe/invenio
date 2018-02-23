function [  ] = batch_resized_jpgs(  )
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
        default.source_path = '/Users/Kelbe/Desktop/';
end
 
fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
  
%% Update paths
 
options.source_path = default.source_path;

%% Choose source file, parent files, and output directory
% Choose parent files
cd(default.source_path);
[folders] = uipickfiles('Prompt','Please choose folders to process');
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

     
    %% Load images and process
    % Load Images
        
      
     
    n.f = size(folders,2);
    
    for f = 1:n.f;
        fprintf('\nWorking on folder %g of %g: folio %s \n', f, n.f, folders{f}(end-5:end));
        dir_tiff = sprintf('%s/tiff/',folders{f});
        dir_jpg = sprintf('%s/jpg/',folders{f});
        if ~exist(dir_jpg);
            mkdir(dir_jpg);
        end
        cd(dir_tiff); 
        D = dir();
        n_d = numel(D);
        for d = 1:n_d;
            filename_I = D(d).name;
            if strcmp(filename_I(1), '.')
                continue
            end 
            filepath_I = sprintf('%s%s',dir_tiff,filename_I);
            filename_J = sprintf('%sjpg', filename_I(1:end-3));
            filepath_J = sprintf('%s%s',dir_jpg,filename_J);
            if exist(filepath_J,'file');
                continue
            end
            I = imread(filepath_I);
            J = imresize(I,.4);
            iminfo = imfinfo(filepath_I);
            bitdepth = iminfo.BitsPerSample(1);
            maxval = 2^bitdepth;
            J = uint8(256*double(J)/maxval);
            imwrite(J,filepath_J,'jpeg','Quality', 50);
        end
     end
    
    fprintf('\n***********************************************************\n');
fprintf('\nCompleted successfully\n');

end
