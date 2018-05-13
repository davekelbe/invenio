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
        default.source_path = pwd;
        default.processed_dir = default.source_path;
    case 'dave'
        slash = '/';
        rmcall = 'rm';
        movetonewfolder = 0;
        default.source_path = '/Users/Kelbe/';
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
lowfrac = .0002;
TOL = [lowfrac 1-lowfrac];
fprintf('\n%g%% Contrast Stretch\n', 100*(lowfrac*2));

%% Choose source file, parent files, and output directory
% Choose parent files
cd(default.source_path);
%[SourceFile, source_path] = uigetfile('*.tif','Please choose a file for submission', 'MultiSelect', 'off');

[files] = uipickfiles('Prompt','Please choose 2-3 images to combine into Photoshop RGB', 'FilterSpec', '*.tif','REFilter', '^[^\.]', 'REDirs', true);
fprintf('\n***********************************************************\n');

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
k = strfind(outfilename,slash);
%outfilename = outfilename(1:k(end-1));
%outfilename = sprintf('%sRGB%s%s',outfilename,slash,temp);
for f = 2:n.files;
    temp_outfilename = sprintf('%s',files{f}(1:end-4));
    k = strfind(temp_outfilename,slash);
    temp = temp_outfilename(k(end)+13:end);
    outfilename = sprintf('%s+%s',outfilename,temp);
end
outfilename = sprintf('%s_DJK_RGB.tif',outfilename);

imwrite(I,outfilename,'tif')


fprintf('\n***********************************************************\n');
fprintf('\nCompleted successfully\n');


end









