function [  ] = separate_ENVI_bands2(  )
%SEPARATE_ENVI_BANDS Writes 2% Contrast Adjusted, 8-bit or 16-bit TIF and JPEG files for
%Photoshop
%
%   There is no input to this function. Typing separate_ENVI_bands in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Separate ENVI Bands Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of TechnologyProce
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
        default.source_path = pwd;
        default.processed_dir = default.source_path;
end

fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
fprintf('Save:       %s\n',default.processed_dir);

%% Update paths

options.source_path = default.source_path;
options.processed_dir = default.processed_dir; %GUI to customize
%% Set Contrast Stretch Level
 
% 2% Contrast Stretch
lowfrac = .005; 
TOL = [lowfrac 1-lowfrac];
fprintf('\n%g%% Contrast Stretch\n', 100*(lowfrac*2));

%% Choose source file, parent files, and output directory
% Choose parent files
cd(default.source_path);
%[SourceFile, source_path] = uigetfile('*.tif','Please choose a file for submission', 'MultiSelect', 'off');
 
[files] = uipickfiles('Prompt','Please choose envi images to process', 'FilterSpec', '*.img');
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');
 
n.files = size(files,2);
for f = 1:n.files
    I = enviread(files{f});
    I.z = double(I.z);
    n.bands  = I.info.bands;
    if n.bands >1; 
        for b = 1:n.bands;
            fprintf('Working on band %g of %g for image %g of %g\n',b, n.bands, f, n.files);
            slashindex = strfind(files{f},slash);
            options.source_updir = files{f}(1:slashindex(end-1));
            ix_slash = strfind(options.source_updir, '/');
            name = options.source_updir(ix_slash(end-1)+1:end-1);
            file_subset = files{f}(slashindex(end):end-4);
            outfilenametiff = sprintf('%s%s+tiff%s_b%g.tif',options.source_updir,name,file_subset,b);
            outfilenamejpeg = sprintf('%s%s+jpg%s_b%g.jpg',options.source_updir,name,file_subset,b);      
            fprintf('       %s\n',outfilenametiff);
         %   fprintf('       %s\n',outfilenamejpeg);
            I.z(:,:,b) = I.z(:,:,b) + abs(min(min(I.z(:,:,b))));
            I.z(:,:,b) = I.z(:,:,b)./max(max(I.z(:,:,b)));
            A = imadjust(I.z(:,:,b),stretchlim(I.z(:,:,b),TOL),[]);
            A = double(A);
            A = A-min(A(:));
            A = A./max(A(:));
            imwrite(uint16(65536*A),outfilenametiff,'tif')   
            Aresize = imresize(A, 0.4);
            imwrite(uint8(256*Aresize),outfilenamejpeg,'jpeg','Quality', 50)          
        end
    end
    
    
end

fprintf('\n***********************************************************\n');
fprintf('\nCompleted successfully\n');


end













