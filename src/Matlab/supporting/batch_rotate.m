function [  ] = batch_rotate(  )
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
        default.source_path = '/';
        default.processed_dir = default.source_path;
end

fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
fprintf('Save:       %s\n',default.processed_dir);
 
%% Update paths
 
options.source_path = default.source_path;
options.processed_dir = default.processed_dir; %GUI to customize

%% Choose source file, parent files, and output directory
% Choose parent files
cd(default.source_path);
[files] = uipickfiles('Prompt','Please choose files to process');
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

     
    %% Load images and process
    % Load Images
    
    
    
    n.b = size(files,2);
    
    for b = 1:n.b
        fprintf('\nWorking on %g of %g\n', b, n.b);
        I = imread(sprintf('%s%s',files{b}));
        I = imrotate(I,-90); %90 = rotate left (counter clockwise) 
        outfilename.jpeg = sprintf('%s.jpg',files{b}(1:end-4));
        outfilename.tif = sprintf('%s.tif',files{b}(1:end-4));
        if strcmp(files{b}(end-2:end),'jpg')
            imwrite(I,outfilename.jpeg,'jpeg','Quality', 50);
        elseif strcmp(files{b}(end-2:end),'tif')
           % imwrite(I,outfilename.tif,'tif');
            imwrite2tif(I, [], outfilename.tif, 'uint16');
        end
    end
    
    fprintf('\n***********************************************************\n');
fprintf('\nCompleted successfully\n');

end
