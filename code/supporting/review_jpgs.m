% Review processed 
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

[m_filepath] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose folders to review');
%m_path_upper = m_path_upper';
n_m = numel(m_filepath); 
%for m = 1:n_m;
%    m_path_upper{m} = sprintf('%s%s',m_path_upper{m},info_slash);
%end
fprintf('\n***********************************************************\n');
fprintf('Folios to review: \n');

% Find upper level directory
ix_slash = strfind(m_filepath{1},info_slash);
path_source = m_filepath{1}(1:ix_slash(end-1));


% Update previous directory 
fid = fopen(filepath_source_previous, 'w+');
fprintf(fid, '%s', path_source); 
fclose(fid);

% Choose reference image
%path_ref_previous = '/Volumes/Timothy/KTK-jpg/';
%[ref_filepath] = uipickfiles('filterspec',path_ref_previous,'Prompt','Please choose reference image');
%ref_filepath = ref_filepath{1};
%Iref = imread(ref_filepath);
%[nx_ref, ny_ref,~] = size(Iref);

figure('position', [592 123 1315 982]);
for m = 1:n_m;
   I = imread(m_filepath{m});
   if m==1;
     [nx_I, ny_I] = size(I);
      xmin_I = 1;
      ymin_I = 1;
      xmax_I = nx_I;
      ymax_I = ny_I;
   end
   Icrop = I(xmin_I:xmax_I, ymin_I:ymax_I);
   [n_xcrop, n_ycrop] = size(Icrop);
   handles.himg = imagesc(Icrop);
   %axis image
   colormap gray
   ascii_code = 0;
   while ascii_code ~=99; % continue
       ascii_code = getkey();
       if ascii_code ~= 99; 
           if ascii_code == 99;
               key = 'c';
           end
           if ascii_code == 116;
               key = 't';
           end
           if ascii_code == 114;
               key = 'r';
           end
           if ascii_code == 107;
               key = 'k'; 
           end
           if ascii_code == 105;
               key = 'i';
           end
           if ascii_code == 8;
               key = 'd';
           end     
           if ascii_code == 100;
               key = 'd';
           end 
           if ascii_code == 120;
               key = 'x';
           end
           switch key;
               case 't';
                   xmin_ref = floor(xmin_I * nx_ref ./ nx_I);
                   xmax_ref = floor(xmax_I * nx_ref ./ nx_I);
                   ymin_ref = floor(ymin_I * ny_ref ./ ny_I);
                   ymax_ref = floor(ymax_I * ny_ref ./ ny_I);
                   refcrop = Iref(xmin_ref:xmax_ref, ymin_ref:ymax_ref);
                   handles.himg = imagesc(refcrop);
                   colormap('gray');
               case 'r';               
                   handles.himg = imagesc(Icrop);
                   colormap('gray');       
               case 'k';
                   % Crop
                   fprintf('\n Crop\n');
                   h_roi = imrect;
                   mask = createMask(h_roi,handles.himg);
                   [xis, yis] = find(mask);
                   xmin_I = min(xis);
                   ymin_I = min(yis);
                   xmax_I = max(xis);
                   ymax_I = max(yis);
                   Icrop = I(xmin_I:xmax_I, ymin_I:ymax_I);
                   [n_xcrop, n_ycrop] = size(Icrop);
                   handles.himg = imagesc(Icrop);
                   colormap('gray');
               case 'i';
                   % invert
                   fprintf('\n Invert\n');
                   Icrop = imcomplement(Icrop);
                   handles.himg = imagesc(Icrop);
                   colormap(gray);
                   I = imcomplement(I);
                   imwrite(I,m_filepath{m},'jpeg','Quality', 50);
               case 'd';
                   % delete
                   fprintf('\n Delete\n');
                   command = sprintf('rm %s', m_filepath{m});
                   fprintf('\n%s\n', command);
                   [a,b] = system(command);
                   key = 'c';
                   break
               case 'c';
                   % Continue
                   fprintf('\n Continue\n');
               case 'x';
                   error('Halted by user')
                   break;
           end
       end
   end
  foo = 1;
   
end


