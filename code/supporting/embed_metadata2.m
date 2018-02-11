function [  ] = embed_metadata2( varargin )
%EMBED_METADATA Embed metadata into a tif image according to <ap> standard
%   This function prepares a processed image file for subission into the
%   repository by saving a correctly rotated JPEG and TIF version of a
%   user-selected input file into a directory of choice. Metadata is
%   embedded into the output file using exiftool and json. 
% 
%   There is no input to this function. Typing embed_metadata in the
%   command line brings up a series of user interfaces which allow the user
%   to select an input file for submission, select the parent files, and
%   confirm and/or edit the metadata fields. Finally, the user chooses an
%   output directory, to which the prepared file is saved as both a JPEG
%   and TIFF version. An option also allows the saving of individual bands
%   of an RGB image with metadata automatically attached. Additionally, a
%   log file is created/updated which stores pertinent information
%   regarding the preparation of processed files. It is recommended that
%   the user change the source code directly to adjust default paths and
%   values which will go into the various metadata fields. 
%
%           
% Metadata Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%
% V0.0 - Initial Version - December 29 2012
% V0.1 - Updated field for free text to HostComputer; Fixed saving band 1
% V0.2 - Added user case for switching between Dave/Roger
%      - Added ability to load and convert 16 bit image to 8 bit
%      - Fixed < sign in Get Rotation information
%      - 
%
% This function adds metadata to a processed image file and saves a
% correctly-rotated TIFF and JPEG image to a folder specified by the user.
%
% Requirements:
%   *<ap> custom EXIF namespace requires .ExifTool_config file be placed in
%   home directory
%   *Input images should be TIFF format
%
%
% Tips:
%   * Press ctrl+c to cancel execution and restart
%   *Set default paths in source code for efficiency (Line 56-59)
%   *Set <ap> and Exif metadata defaults in source code for efficiency
%   (Line 158-177)
%   *Input image should be unrotated, i.e., capture rotation. Rotation is
%    automatically determined from EXIF Orientation tag
% 
fprintf('\n***********************************************************\n');
fprintf('Tips\n');
fprintf('            Press ctrl+c to cancel execution and restart\n');
fprintf('            *PLEASE change default paths in source code (line 57-60)\n');
fprintf('            Change default metadata parameters in source code (line 158-174)\n');
fprintf('            Close log.csv file before continuing \n');
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
        exiftoolcall = 'exiftool';
        rmcall = 'rm';
        cpcall = 'cp';
        default.processing_type = 'ICA';
        info.Creator = 'Dave Kelbe';
        %info.Contributor = 'Dave Kelbe';
        %default.source_path = '/Volumes/Tarsus/Processed/';
        default.source_path = '/Volumes/Timothy/';

        default.parentpath = '/Volumes/Vatican/Day2raw/Flattened/';
        %default.parentpath = '/Volumes/Babel/2013-05/';
        %'/Volumes/cyclone.cis.rit.edu - -dirs-grad-djk2312/other/data/Flattened_Images/2011-12';
        default.outdir =  '/Volumes/Vatican/Delivery/';
        %default.outdir =  '/Volumes/Philippi/Delivery/DJK/2014-01-09/data/';
        default.excelpath = '/Users/Kelbe/Documents/EMEL/Sinai/Processed/Log/'; 
%        default.exifpath = '/Volumes/Babel/2013-05/';
        default.exifpath = '/Volumes/Vatican/Day2raw/Flattened/';
 
end
 
fprintf('\n***********************************************************\n');
fprintf('Setting Default Paths \n');
fprintf('Source:     %s\n',default.source_path);
fprintf('Parent:     %s\n',default.parentpath);
fprintf('Save:       %s\n',default.outdir);

%% Choose source file, parent files, and output directory

% Choose source file
cd(default.source_path);
[SourceFile, source_path] = uigetfile('*.tif','Please choose a file for submission', 'MultiSelect', 'off');

%% PARAMETERS

type = varargin{1};
switch type
    case '3band'
        info.DAT_File_Processing = sprintf('RGB image created from ICA bands. Hue and levels adjustment in Photoshop.\n %s',SourceFile);%sprintf('%s',SourceFile);
        processing_type = 'ICA_3band';
    case '2band'
        info.DAT_File_Processing = sprintf('Two-band RGB image created from ICA bands (gives less neon colors). Hue and levels adjustment in Photoshop.\n %s',SourceFile);%sprintf('%s',SourceFile);
        processing_type = 'ICA_2band';
    case 'gray'
        info.DAT_File_Processing = sprintf('Grayscale image created from ICA bands.\n %s',SourceFile);%sprintf('%s',SourceFile);
        processing_type = 'ICA_gray';
    case 'color'
        info.DAT_File_Processing = sprintf('Natural color image with contrast enhancement via Spectralon.\n %s',SourceFile);%sprintf('%s',SourceFile);
        processing_type = 'color';
    case 'uv'
        info.DAT_File_Processing = sprintf('UV Fluourescence image made by combining orange, green, and blue fluorescence images under UV illumination.\n %s',SourceFile);%sprintf('%s',SourceFile);
        processing_type = 'UV_fluorescence';
end
%%

% Choose parent files 
%cd(source_path);
%{
cd(default.parentpath); %Roger
filterspec{1} = '*.tif';
filterspec{2} = sprintf('*%s*',SourceFile(1:11));
[parentfile, parentpath] = uigetfile(filterspec,'Please choose parent files', 'MultiSelect', 'on');
%{
parentfile{1} = sprintf('%s+MB365UV_001.tif',SourceFile(1:11));   % Roger
parentfile{2} = sprintf('%s+MB455RB_002.tif',SourceFile(1:11));   % Roger
parentfile{3} = sprintf('%s+MB470LB_003.tif',SourceFile(1:11));   % Roger
parentfile{4} = sprintf('%s+MB505Cy_004.tif',SourceFile(1:11));   % Roger
parentfile{5} = sprintf('%s+MB535Gr_005.tif',SourceFile(1:11));   % Roger
parentfile{6} = sprintf('%s+MB570Am_006.tif',SourceFile(1:11));   % Roger
parentfile{7} = sprintf('%s+MB625Rd_007.tif',SourceFile(1:11));   % Roger
parentfile{8} = sprintf('%s+MB700IR_008.tif',SourceFile(1:11));   % Roger
parentfile{9} = sprintf('%s+MB735IR_009.tif',SourceFile(1:11));   % Roger
parentfile{10} = sprintf('%s+MB780IR_010.tif',SourceFile(1:11));   % Roger
parentfile{11} = sprintf('%s+MB870IR_011.tif',SourceFile(1:11));   % Roger
parentfile{12} = sprintf('%s+MB940IR_012.tif',SourceFile(1:11));   % Roger
parentfile{13} = sprintf('%s+WBRBR25_020.tif',SourceFile(1:11));   % Roger
parentfile{14} = sprintf('%s+WBRBG61_021.tif',SourceFile(1:11));   % Roger
parentfile{15} = sprintf('%s+WBRBB47_022.tif',SourceFile(1:11));   % Roger
%}
parentpath = default.parentpath;% Roger
%}
 
% Parse to choose parent files
parentpath = default.parentpath;
if nargin == 1;
    bands = parse_filename( SourceFile );
else
    bands = varargin{2};
end
bands = [3 5 7 9 12 13 15 16 17 26 ];

parentfile = cell(1,numel(bands));
for b = 1:numel(bands);
    cd(sprintf('%s%s/',parentpath, SourceFile(1:11)));
    D = dir(sprintf('*_%03.0f_F*',bands(b)));
    if numel(D)>1;
        D = D(2);
    end
    parentfile{b} = D.name;
end

% Print output
fprintf('\n***********************************************************\n');
fprintf('Selected source and parent files\n')
fprintf('Source:      %s%s\n', source_path,SourceFile);
fprintf('\n***********************************************************\n');
fprintf('\nChoosing parent files\n');
fprintf('             %s\n', parentpath);
fprintf('      \n' );
if iscell(parentfile)
    for i = 1:size(parentfile,2)
        fprintf('             %s\n', parentfile{i});
    end
else
    fprintf('             %s\n', parentfile);
end

%% Read Image
fprintf('\n***********************************************************\n');
fprintf('Reading source image\n');
fprintf('\n***********************************************************\n');

I = imread([source_path,SourceFile],'tif');
if size(I,3) > 3;
    I = I(:,:,1:3);
end
%% Set up Metadata 

fprintf('Reading Metadata\n');
fprintf('\n***********************************************************\n');

parentfile_capture = cell(size(parentfile));
for f = 1:size(parentfile,2);
    parentfile_capture{f} = sprintf('%sF.tif',parentfile{f}(1:24));
end

%Get list of 7-character wavelength codes and order by capture sequence
if iscell(parentfile)
    parentfilepath = sprintf('%s%s',parentpath,parentfile{1});
    nparent = size(parentfile,2);
    wavelength = cell(nparent,1);
    bandindex = zeros(nparent,1);
    for i=1:nparent;
        wavelength{i}=parentfile{i}(13:19);
        bandindex(i) = str2double(parentfile{i}(21:23));
    end
    [~,index] = sort(bandindex);
    wavlength = wavelength(index);
    parentfile = parentfile(index);
else
    parentfilepath = sprintf('%s%s',parentpath,parentfile);
    wavelength = parentfile(13:19);
end

foo = parentfile{1};
    parentfilepath = sprintf('%s%s/%s',parentpath,foo(1:11),parentfile{1});
%exiffilepath = sprintf('%sF.tif',parentfilepath(1:end-11));
%options.exifpath = sprintf('%s%s/',default.exifpath,SourceFile(1:11)); %28
%options.exifpath = default.exifpath; %21
exiffilepath = parentfilepath;
%cd(options.exifpath);
%D = dir(sprintf('%s*',SourceFile(1:11)));
%exiffilepath = sprintf('%s%s',options.exifpath,D(32,:).name);
% Get IPTC information
[~,exifout1] = system(sprintf('%s -s -ObjectName "%s"',exiftoolcall,exiffilepath));
[~,exifout2] = system(sprintf('%s -s -Source "%s"',exiftoolcall,exiffilepath));
[~,exifout3] = system(sprintf('%s -s -Keywords "%s"',exiftoolcall,exiffilepath));
[~,exifout4] = system(sprintf('%s -s -Orientation "%s"',exiftoolcall, exiffilepath));
[~,exifout5] = system(sprintf('%s -s -BitsPerSample "%s"',exiftoolcall, exiffilepath));
 
if numel(exifout4)==0;
    warning('Parent file has no Orientation information');
    warning('Double check GUI to confirm correct rotation: upper text UP');
    rotation_angle = 0;
end
k1 = strfind(exifout1, ':');
k2 = strfind(exifout2, ':');
k3 = strfind(exifout3, ':');
if numel(exifout4)>0;
k4 = strfind(exifout4, ':');
end
k5 = strfind(exifout5, ':');


info.ObjectName = strtrim(exifout1(k1+2:end));
info.Source = strtrim(exifout2(k2+2:end));
info.Keywords = strtrim(exifout3(k3+2:end));
if numel(exifout4)>0;
exifrotation = strtrim(exifout4(k4+2:end));
end
BitsPerSample = str2num(strtrim(exifout5(k5+2:end)));

% Get manuscript and folio information
IDX = strfind(info.ObjectName, ',');
manuscript = info.ObjectName(1:IDX-1);
folio = info.ObjectName(IDX+2:end);


% Get rotation information
if numel(exifout4)>0;
switch exifrotation;
    case 'Horizontal (normal)';
        rotation_angle = 0;
    case 'Rotate 180';
        rotation_angle = 180;
    case 'Rotate 90 CW'
        rotation_angle = 90;
    case 'Rotate 270 CW';
        rotation_angle = 270;
end
end

% Set up default <ap> and Exif metadata (EDIT for efficiency)
info.SourceFile = sprintf('%s%s',source_path,SourceFile);
info.ID_Parent_File = parentfile_capture;
info.DAT_Bits_Per_Sample = '8';
info.DAT_Samples_Per_Pixel = num2str(size(I,3));
info.DAT_File_Processing_Rotation = sprintf('%s',num2str(rotation_angle));
info.DAT_Joining_Different_Parts_of_Folio = 'false';
info.DAT_Joining_Same_Parts_of_Folio = 'true';
info.DAT_Processing_Comments = 'none'; 
if numel(exifout4)==0;
    info.DAT_Processing_Comments = 'Warning: Confirm correct orientation (Upper text UP)';
end
info.DAT_Processing_Program{1} = 'MATLAB R2012a (7.14.0.739)';
info.DAT_Processing_Program{2} = 'ENVI 5.0';
info.DAT_Processing_Program{3} = 'Adobe Photoshop CS4';
info.DAT_Software_Version = 'See DAT_Processing_Program';
info.DAT_Type_of_Contrast_Adjustment = 'none';
info.DAT_Type_of_Image_Processing = 'enhancement';
info.HostComputer = sprintf('Local directory: %s, Local filename: %s',source_path, SourceFile);

%% Create Output filename
% Output filename should have the following format:
% <SHOOT_LIST><SHOT_SEQ><PROCESSOR><PROCESSING_TYPE><MODIFIERS>.<EXT>
fprintf('Creating output filename\n');

% Get shoot list and shot sequence
shoot_list = info.Source(1:4);
shot_sequence = info.Source(6:11);

% Get initials of creator
if ~(strcmp(info.Creator,'Dave Kelbe') ||...
        strcmp(info.Creator,'Roger Easton')||...
        strcmp(info.Creator,'Bill Christens-Barry'))
    error('Creator not recognized');
end
switch info.Creator;
    case 'Dave Kelbe';
        processor = 'DJK';
    case 'Roger Easton';
        processor = 'RLE';
    case 'Bill Christens-Barry';
        processor = 'WCB';
end



% Check uniqueness of proposed filename and add modifier as necessary
% Set base filename without modifier
outfilename_base = sprintf('%s_%s_%s_%s',shoot_list,shot_sequence,...
    processor, processing_type);

% Check existence of base filename in output directory
D = dir(sprintf('%s%s*.tif',default.outdir,outfilename_base));
existing_files = cell(size(D,1),1);
if size(D)>0
    for i = 1:size(D,1);
        existing_files{i} = D(i).name;
    end
    IX = strfind(existing_files{1},processing_type);
    existing_files = char(existing_files);
    used_modifiers = existing_files(:,IX+numel(processing_type)+1:IX+numel(processing_type)+2);
    if strcmp(used_modifiers,'RG')
        used_modifiers = '';
    end
    nmod = size(used_modifiers,1);
    used_indices = zeros(nmod,1);
    for i = 1:nmod;
        used_indices(i) = str2num(used_modifiers(i,:));
    end
    if strcmp(used_modifiers,'RG')
        used_indices = 0;
    end
else
    used_indices = 0;
end

% Set modifier 
modifier_unique = sprintf('%02d',max(used_indices)+1);


fprintf('             %s_%s\n', outfilename_base,modifier_unique);
fprintf('\n***********************************************************\n');

% Append filename with final modifier (RGB or single band?)
% if strcmp(info.DAT_Samples_Per_Pixel, '3');
%     modifier_RGB = sprintf('_RGB');
% elseif strcmp(info.DAT_Samples_Per_Pixel, '1');
%     modifier_RGB = '';
% else
%     error('Image has the wrong number of bands');
% end

outfilename = sprintf('%s_%s_%s_%s_%s',shoot_list,shot_sequence,...
    processor, processing_type,modifier_unique);
%outfilename = sprintf('%s_%s_%s_%s',shoot_list,shot_sequence,...
%    processor, processing_type);

%% Prepare image and reference image for visual inspection

fprintf('Rotating Image\n');
fprintf('\n***********************************************************\n');

% Rotate source image 
%I = imrotate(I,-rotation_angle);

fprintf('Resizing image for viewing\n');
fprintf('\n***********************************************************\n');

%{
%Resize source image for viewing
I2 = I(1:10:end,1:10:end,:);

% Load reference image for checking rotation
if numel(exifout4)>0;
Ref = imread(parentfilepath);
Ref = Ref(1:50:end,1:50:end,:);
Ref = imadjust(Ref);
Ref = imrotate(Ref,-rotation_angle);
else
    Ref = 0;
end

%% Call GUI

fprintf('Opening GUI for saving\n');
fprintf('\n***********************************************************\n');
  
[info.ObjectName,...
    info.Source,...
    info.DAT_File_Processing,...
    info.ID_Parent_File,...
    info.DAT_Processing_Comments,...
    info.DAT_Type_of_Contrast_Adjustment,...
    info.Creator,...
    info.DAT_Processing_Program,...
    outdir,...
    outfilename,...
    saveTF,...
    nrot] = gui_metadata(I2,info,outfilename,default.outdir,Ref);
%% Rotate again if necessary
if nrot > 0
    
    I = imrotate(I,-90*nrot);
end
%} 
%% Change to 8-bit if necessary
if BitsPerSample ~= 8;
    I = double(I);
    I = uint8(I./max(I(:))*256);
end 
outdir = default.outdir;
if ~exist(outdir,'dir');
    mkdir(outdir);
end
%% Save images
fprintf('Saving images and embedding metadata\n');
%TIF
info.SourceFile = sprintf('%s%s.tif',outdir,outfilename);
%LocalOutFile = sprintf('%s%s%s.tif',source_path,outfilename,modifier_RGB);
imwrite(I,info.SourceFile,'tif','Compression', 'none')
fprintf('             %s\n', info.SourceFile);
jsonfilename = sprintf('%s.json', info.SourceFile);
[~] = savejson('', info, jsonfilename);
[~,~] = system(sprintf('%s -j=%s -overwrite_original %s',exiftoolcall, jsonfilename,info.SourceFile));
[~,~] = system(sprintf('%s %s',rmcall, jsonfilename));
%[~,~] = system(sprintf('%s %s %s',cpcall, info.SourceFile, LocalOutFile));


%JPEG
info.SourceFile = sprintf('%s%s.jpg',outdir,outfilename);
%LocalOutFile = sprintf('%s%s%s.jpg',source_path,outfilename,modifier_RGB);
imwrite(I,info.SourceFile,'jpeg','Quality', 100)
fprintf('             %s\n', info.SourceFile);
jsonfilename = sprintf('%s.json', info.SourceFile);
[~] = savejson('', info, jsonfilename);
[~,~] = system(sprintf('%s -j=%s -overwrite_original %s',exiftoolcall, jsonfilename,info.SourceFile));
[~,~] = system(sprintf('%s %s',rmcall, jsonfilename));
%[~,~] = system(sprintf('%s %s %s',cpcall, info.SourceFile, LocalOutFile));


% Save individual bands if necessary
%{
if saveTF;
    DAT_File_Processing = info.DAT_File_Processing;
    info.DAT_Samples_Per_Pixel = 1;
    for i = 1:3;
        info.DAT_File_Processing = sprintf('Grayscale Band %u of %s%s: %s',i,outfilename, modifier_RGB,DAT_File_Processing);
        info.SourceFile = sprintf('%s%s_%s.tif',outdir,outfilename,num2str(i));
        imwrite(I(:,:,i),info.SourceFile,'tif','Compression', 'none')
        fprintf('             %s\n', info.SourceFile);
        jsonfilename = sprintf('%s.json', info.SourceFile);
        [~] = savejson('', info, jsonfilename);
        [~,~] = system(sprintf('%s -j=%s -overwrite_original %s',exiftoolcall, jsonfilename,info.SourceFile));
        [~,~] = system(sprintf('%s %s',rmcall, jsonfilename));
        
        info.SourceFile = sprintf('%s%s_%s.jpg',outdir,outfilename,num2str(i));
        imwrite(I(:,:,i),info.SourceFile,'jpeg','Quality', 100)
        fprintf('             %s\n', info.SourceFile);
        jsonfilename = sprintf('%s.json', info.SourceFile);
        [~] = savejson('', info, jsonfilename);
        [~,~] = system(sprintf('%s -j=%s -overwrite_original %s',exiftoolcall, jsonfilename,info.SourceFile));
        [~,~] = system(sprintf('%s %s',rmcall, jsonfilename));
    end
end
%}

% %% Save Excel File
% fprintf('\n***********************************************************\n');
% fprintf('Saving excel file\n');
% 
% date_str = date;
% if saveTF;
%     n_images = 4;
% else
%     n_images = 1;
% end 
% A = cell(n_images,2);
% for i = 1:n_images;
%     A{i,1} = date_str;
%     A{i,3} = manuscript;
%     A{i,4} = folio;
%     A{i,5} = SourceFile; 
%     A{i,6} = source_path;
%     A{i,7} = outdir; 
%     A{i,8} = shoot_list;
%     A{i,9} = shot_sequence; 
%     A{i,11} = strrep(info.Creator,',','.');
%     A{i,13} = info.DAT_File_Processing_Rotation;
%     A{i,14} = info.DAT_Type_of_Contrast_Adjustment;
%     A{i,15} = strrep(info.DAT_Processing_Comments,',','.');
% end
% A{1,2} = sprintf('%s%s.tif',outfilename,modifier_RGB);%Filename
% A{1,10} = strrep(info.DAT_File_Processing,',','.');%Description
% A{1,12} = num2str(size(I,3));%Bands
% if saveTF;
%     for i=1:3
%         A{i+1,2} = sprintf('%s_%s.tif',outfilename,num2str(i));%Filename
%         A{i+1,10} = strrep(sprintf('Grayscale Band %u of %s%s: %s',i,outfilename, modifier_RGB,info.DAT_File_Processing),',','.');%Description
%         A{i+1,12} = num2str(1);%Bands
%     end
% end
% 
% if ~exist(sprintf('%sProcessing_Log.csv',default.excelpath));
%     fid = fopen(sprintf('%sProcessing_Log.csv',default.excelpath), 'wt');
%     fprintf(fid, 'Date,Submitted File,Manuscript,Folio,Local File,Local Directory,Submitted Directory,Shoot List, Shot Sequence, Description,Creator, Samples, Rotation, Contrast, Comments');
%     fclose(fid);
% end
% fid = fopen(sprintf('%sProcessing_Log.csv',default.excelpath), 'rt');
% Sheet = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 0, 'Delimiter', ',','CollectOutput',1);
% Sheet = [Sheet{1}; A];
% fclose(fid);
% fid = fopen(sprintf('%sProcessing_Log.csv',default.excelpath), 'wt');
% for i=1:size(Sheet,1)
%     fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', Sheet{i,:});
% end
% fclose(fid);
% fprintf('             Excel file saved to %s\n',default.excelpath);

%cd(current_path)
fprintf('\n***********************************************************\n');
fprintf('Completed Successfully');
fprintf('\n***********************************************************\n');


end

