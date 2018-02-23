function [  ] = batch_make_truecolorRGB( varargin )
%NORMALIZED_ENVI_CUBE Create a normalized ENVI image cube
% 
%   There is no input to this function. Typing reflectance_tiffs in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Reflectance TIFFS  Tool
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
%%
file1 = '*MB625RD*';
file2 = '*MB535GN*';
file3 = '*MB450RB*';

n.files = 3;

if numel(varargin)>0;
    bandsubset = true;
    bands = varargin{1};
else
    bandsubset = false;
end

%% Set User
user = 'dave';
%user = 'roger';

%% Set default paths

%current_path = cd;
%current_path = sprintf('%s/',current_path);
% Change default paths here
%% Set default paths
%current_path = cd;
%current_path = sprintf('%s/',current_path);
% Change default paths here
switch user 
    case 'roger'
        slash = '\';
        exiftoolcall = 'exiftool.pl';
        rmcall = 'del';
        default.processing_type = 'PCA';
        info.Creator = 'Roger Easton';
        info.Contributor = 'Dave Kelbe';
        default.source_path = 'g:\research\StC\';
        default.parentpath = 'g:\research\StC\';
        default.outdir = 'C:\Users\rlepci\Documents\Research\StC\Submission\Images\';
        default.excelpath = 'C:\Users\rlepci\Documents\Research\StC\Submission\Log\';
    case 'dave'
        slash = '/';
        exiftoolcall = '/usr/bin/exiftool';
        rmcall = 'rm';
        cpcall = 'cp';
        info.Creator = 'Dave Kelbe';
        % info.Contributor = 'Dave Kelbe';
        %default.source_path = '/Volumes/Tarsus/Processed/';
        default.source_path = '/Volumes/';
        default.outdir =  default.source_path;
end

fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
fprintf('Save:       %s\n',default.outdir);

if ~exist(default.outdir,'dir');
    mkdir(default.outdir);
end


%% Update paths

%options.processed_dir = default.processed_dir; %GUI to customize

%% Choose source file, parent files, and output directory
% Choose parent files
cd(default.source_path);
[folders] = uipickfiles('Prompt','Please choose folders to process');
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');

slashindex = strfind(folders{1},slash);
options.source_updir = folders{1}(1:slashindex(end));
%fprintf('Path:%s \n',options.source_updir);

%for i = 1:size(folders,2)
%    fprintf('      %s\n',folders{i}(slashindex(end)+1:end));
%end

% Make new folders for reflectance images
%options.processed_updir = '/Volumes/Vienna/Processed/';

for f = 1:size(folders,2)
    fprintf('\n***********************************************************\n');
    fprintf('***********************************************************\n');
    fprintf('####################### %s #######################\n',folders{f}(slashindex(end)+1:end));
    fprintf('***********************************************************\n');
    fprintf('***********************************************************\n');
    
    
    source_path = folders{f};% sprintf('%s%s%s%s',folders{f},slash,'tiff',slash);
    if ~strcmp(source_path(end), slash);
        source_path = sprintf('%s%s',source_path, slash);
    end
    
    slash_ix = strfind(source_path, slash);
    cube = source_path(slash_ix(end-1)+1:end-1);
    
    cd(source_path);
    D = dir(sprintf('*%s',file1));
    
    if numel(D) ==1
        files{1} = sprintf('%s%s',source_path,D.name);
    else
        error('Multiple files matching string');
    end
    
    D = dir(sprintf('*%s',file2));
    if numel(D) ==1 
        files{2} = sprintf('%s%s',source_path, D.name);
    else
        error('Multiple files matching string');
    end
    
    
    if n.files ==3;
        D = dir(sprintf('*%s',file3));
        if numel(D) ==1
            files{3} = sprintf('%s%s',source_path,D.name);
        else
            error('Multiple files matching string');
        end
    end
    
  %  options.shoot_list = SourceFile(1:4);
  %  options.shot_seq = SourceFile(6:11);
    
    
    
    
    if n.files > 3;
        error('Please supply 2 or 3 images');
    end
    
    I1 = double(imread(files{1}));
      % Load first image for spectralon
    h = figure('name','Please choose spectralon');
    
    imagesc(I1);
    colormap jet 
    hFH = imrect();
    % Create a binary image ("mask") from the ROI object.
    mask = hFH.createMask();
    delete(h);
    
    TOL = 1.01;
    spectralon_DC1 = mean(double(I1(mask)));
    I1cal = TOL*I1./spectralon_DC1;
    I2 = double(imread(files{2}));
    spectralon_DC2 = mean(double(I2(mask)));
    I2cal = TOL*I2./spectralon_DC2;
    I3 = double(imread(files{3}));
    spectralon_DC3 = mean(double(I3(mask)));
    I3cal = TOL*I3./spectralon_DC3;
    
    I1cal(I1cal > 1) = 1;
    I2cal(I2cal > 1) = 1;
    I3cal(I3cal > 1) = 1;

       [n.r,n.c] = size(I1);
       RGB = uint16(zeros(n.r,n.c,3));
        RGB(:,:,1) = 65535*I1cal;
        RGB(:,:,2) = 65535*I2cal;
        RGB(:,:,3) = 65535*I3cal;

        RGB = imadjust(RGB,[0; 1],[.15; .85]);
       
    %spectral_DCmax = mean(spectralon_DC(:))+2*std(spectralon_DC(:));
    %LOW_HIGH = stretchlim(spectralon_DC./spectral_DCmax,[0 .99]);
    
    %ref_val =  1*LOW_HIGH(2)*spectral_DCmax;
    
    %lowfrac = .0002;
    %TOL = [lowfrac 1-lowfrac];
     
%     bandname = cell(3,1);
%     for b = 1:n.files
%         A = double(imread(files{b}));
%         if b ==1;
%             [n.r,n.c] = size(A);
%             I = uint16(zeros(n.r,n.c,3));
%         end
%         A = A./max(A(:));
%         A = imadjust(A,stretchlim(A,TOL),[]);
%         A = double(A);
%         A = A-min(A(:));
%         A = A./max(A(:));
%         A = uint16(65536*A);
%         if b == 2 && n.files ==2;
%             I(:,:,2) = A;
%             I(:,:,3) = A;
%             bandname{2} = files{b}(end-16:end-10);
%             bandname{3} = files{b}(end-16:end-10);
%         else
%             I(:,:,b) = A;
%             bandname{b} = files{b}(end-16:end-10);
%         end
%     end
    %outfilename = outfilename(1:k(end-1));
    %outfilename = sprintf('%sRGB%s%s',outfilename,slash,temp);
    
    ix_slash = strfind(folders{f}, slash);
    name = folders{f}(ix_slash(end)+1:end);
    ix_underscore = strfind(name, '_');
    mss = name(1:ix_underscore(end)-1);
    
    ix_slash = strfind(folders{f}, slash);
    dir_top = folders{f}(1:ix_slash(end-1));
    dir_current = sprintf('%sProcessed-%s/%s/',dir_top, mss, name);
    dir_jpg = sprintf('%s%s+jpg/', dir_current,name );
    dir_tif = sprintf('%s%s+tiff/', dir_current,name );
    filepath_jpg = sprintf('%s%s_DJK_true.jpg', dir_jpg, name);
    filepath_tif = sprintf('%s%s_DJK_true.tif', dir_tif, name);

    
    %outpath = '/Volumes/NCM-Pilot/Data/Preliminary Processed Images/DJKColor/';
    %outfilename = sprintf('%s%s_DJK_true.tif',outpath,cube);
    
    imwrite(RGB,filepath_tif,'tif')
    
    RGBjpg = uint8((double(RGB)./65536)*256);
    RGBjpg = imresize(RGBjpg, 0.4);
    imwrite(RGBjpg,filepath_jpg,'jpg', 'quality', 50)

  

end
fprintf('Completed %s successfully\n',folders{f}(slashindex(end)+1:end));

end




