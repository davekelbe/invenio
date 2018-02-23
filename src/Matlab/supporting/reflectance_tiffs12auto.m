function [ w_wavelength, m_wavelength_file, m_wavelength_filepath, m_rotation_angle ] = reflectance_tiffs11( aux )
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
%% Preliminary setup
fprintf('\n***********************************************************\n');
fprintf('Contrast Enhancement: \n');

%options_same_stretch = true;
%filepath_stretchValMax = '/Volumes/Scythica/stretchValMax.mat';
%load(filepath_stretchValMax);
%stretchRef = w_stretchval;



m_path_upper = aux.m_path_upper;
m_folio = aux.m_folio;
m_mss = aux.m_mss;
m_name = aux.m_name;
is_band_subset = aux.is_band_subset;
bands = aux.bands;
info_rmcall = aux.info_rmcall;
info_slash = aux.info_slash;
info_user = aux.info_user;
n_m = aux.n_m;
options_delimiter = aux.options_delimiter;
options_delimiter_wavelength = aux.options_delimiter_wavelength;
options_folder_structure = aux.options_folder_structure;
options_movetonewfolder = aux.options_movetonewfolder;
path_source = aux.path_source;
path_target = aux.path_target;
subpath_tiff_dir = aux.path_tiff_dir;
subpath_jpg_dir = aux.path_jpg_dir;
subpath_matlab_dir = aux.path_matlab_dir;
subpath_envi_dir = aux.path_envi_dir;
exiftoolcall = aux.exiftoolcall;



clear aux
%% Get wavelength band names
cd(m_path_upper{1});
D = dir('*_F.tif');
w_wavelength = remove_hiddenfiles(D);
n_w = numel(w_wavelength);
for w = 1:n_w
    ix_delimiter = strfind(w_wavelength{w},options_delimiter_wavelength);
    w_wavelength{w} = w_wavelength{w}(ix_delimiter:end);
end

% Determine illumination type
w_istx = cellfun(@(x) ~isempty(strfind(x,'TX')), w_wavelength);
w_isref = cellfun(@(x) ~isempty(strfind(x,'MB')), w_wavelength);
w_isuv1 = cellfun(@(x) ~isempty(strfind(x,'WBUV')), w_wavelength);
w_isuv2 = cellfun(@(x) ~isempty(strfind(x,'+W')), w_wavelength);
w_isuv = w_isuv1 | w_isuv2;
w_isrb = cellfun(@(x) ~isempty(strfind(x,'WBRB')), w_wavelength);
w_isrs = cellfun(@(x) ~isempty(strfind(x,'RS')), w_wavelength);
w_isre = cellfun(@(x) ~isempty(strfind(x,'RE')), w_wavelength);

% Initialize ref_val and shutter_speed arrays
filepath_ref_val = sprintf('%s%s_%s_refval.mat',subpath_matlab_dir{1},m_mss{1}, m_folio{1});
filepath_shutter_speed = sprintf('%s%s_%s_shutter_speed.mat',subpath_matlab_dir{1},m_mss{1}, m_folio{1});
filepath_aperture = sprintf('%s%s_%s_aperture.mat',subpath_matlab_dir{1},m_mss{1}, m_folio{1});


clear D ix_delimiter w
% Output
% w_wavelength              - cell array of wavelength identifiers
% w_is*                     - logical identifying illumination type
% n_w                       - number of wavelength bands
% filepath_*                - paths to reference data
%% Select mask of parchment
m_path_upper_mask = cell(n_m,1);
m_wavelength_filepath = cell(n_m,1);
m_wavelength_file = cell(n_m,1);
for m = 1:n_m
    
    % Determine wavelength filepath
    m_wavelength_filepath{m} = cell(n_w,1);
    m_wavelength_file{m} = cell(n_w,1);
    switch options_folder_structure
        case 'mss_folio'
            for w = 1:n_w;
                m_wavelength_filepath{m}{w} = sprintf('%s%s%s', ...
                    m_path_upper{m},m_name{m},w_wavelength{w});
                m_wavelength_file{m}{w} = sprintf('%s%s', ...
                    m_name{m},w_wavelength{w});
            end
    end
    
    %specmask = detect_spectralon(filepath_uv)
    % Manually determine mask of parchment
    m_path_upper_mask{m} = sprintf('%s%s_parchment_mask.tif',subpath_tiff_dir{m}, m_name{m});
    if ~exist(m_path_upper_mask{m}, 'file')
        I = imread(m_wavelength_filepath{m}{1});
        I = imadjust(I,stretchlim(I), []);
        h = figure('name','Please choose parchment');
        imagesc(I);
        hFH = imfreehand();
        % Create a binary image ("mask") from the ROI object.
        mask = hFH.createMask();
        delete(h);
        imwrite(mask,m_path_upper_mask{m}, 'tif')
    end
end

clear m w I h hfH mask
% Output
% m_wavelength_filepath     - cell array of wavelength filepaths
% m_path_upper_mask           - mask
% n_w              - number of wavelength bands
% filepath_*                - paths to reference data
fprintf('Working on :\n');
%% Normalize each image and save
filepath_rotation_angle = cell(n_m,1);
m_rotation_angle = zeros(n_m,1);
for m = 1:n_m
    fprintf('                 \t\t%s\n', m_name{m});
    
    w_shutter_speed = zeros(n_w,1);
    w_aperture = zeros(n_w,1);
    w_stretchval = zeros(n_w,1);
    
    
    % Go through each wavelength image
    for w = 1:n_w;
       % continue 
        % Continue to next image if raking
        if w_isrs(w) || w_isre(w); continue; end
        
        % If the stretched file does not exist, begin
        filepath_tif = sprintf('%s%s_stretch.tif',...
            subpath_tiff_dir{m}, m_wavelength_file{m}{w}(1:end-4));
        filepath_jpg = sprintf('%s%s_stretch.jpg',...
            subpath_jpg_dir{m}, m_wavelength_file{m}{w}(1:end-4));
        filepath_summary_txt = sprintf('%s%s_%s_summary.txt',subpath_matlab_dir{m},...
            m_mss{m},m_folio{m});
        if  exist(filepath_tif,'file') && exist(filepath_jpg,'file')% &&...
               % exist(filepath_summary_txt, 'file')
            continue
        end
        
        I = double(imread(m_wavelength_filepath{m}{w}));
        %I_maxval = max(I(:));
        
        % Median filter 
        if w_isuv(w) || w_isrb(w);
            I = medfilt2(I,[3,3]);         
        end
        
        % Get shutter speed
        filepath_nospace = strrep( m_wavelength_filepath{m}{w}, ' ', '\ ');
        command = sprintf('%s -ShutterSpeed %s', exiftoolcall, filepath_nospace);
        [~, exifout] = system(command);
        k = strfind(exifout, ': ');
        if ~exist(filepath_shutter_speed, 'file');
            w_shutter_speed(w) = str2double(strtrim(exifout(k+2:end)));
        end;
        w_shutter_speed(w) = str2double(strtrim(exifout(k+2:end)));
        
        % Get aperture
        command = sprintf('%s -Aperture %s', exiftoolcall, filepath_nospace);
        [~, exifout] = system(command);
        k = strfind(exifout, ': ');
        if ~exist(filepath_shutter_speed, 'file');
            w_aperture(w) = str2double(strtrim(exifout(k+2:end)));
        end;
        w_aperture(w) = str2double(strtrim(exifout(k+2:end)));
        clear command exifout k
        
        % Determine stretch percentage based on illumination type
        if w_istx(w);
            TOL = [0 .99];
            sat = 0.75;
        elseif w_isref(w);
            TOL = [0 .99];
            sat = 0.75;
        elseif w_isuv(w);
            TOL = [0 .99];
            sat = 0.75;
        end
        
        % Determine stretch value
        mask = imread(m_path_upper_mask{m}, 'tif');
        parch1d = I(mask);
        max_parch = max(parch1d);
        Jparch1d = parch1d./max_parch;
        upper_stretch1 = stretchlim(Jparch1d,TOL);
        upper_stretch2(2) = upper_stretch1(2)*(1/sat);
        w_stretchval(w) =  upper_stretch2(2) * max_parch;
        clear mask parch1d max_parch Jparch1d upper_stretch1 upper_stretch2
        
        
        if exist('stretchValMax', 'var');
            J = I./stretchValMax(w);
        else 
          %  break;
            J = I./w_stretchval(w);
        end
        
        J(J>1) = 1;
        
        % Get rotation
%         if w==1;
%             command = sprintf('%s -Orientation %s', exiftoolcall, filepath_nospace);
%             [~, exifout] = system(command);
%              k = strfind(exifout, ': ');
%             if numel(exifout)>0;
%                 exifrotation = strtrim(exifout(k+2:end));
%             end
%             
%             if numel(exifout)>0;
%                 switch exifrotation;
%                     case 'Horizontal (normal)';
%                         m_rotation_angle(m) = 0;
%                     case 'Rotate 180';
%                         m_rotation_angle(m) = 180;
%                     case 'Rotate 90 CW'
%                         m_rotation_angle(m) = 90;
%                     case 'Rotate 270 CW';
%                         m_rotation_angle(m) = 270;
%                 end
%             end
%             filepath_rotation = sprintf('%srotation.mat', subpath_matlab_dir{m});
%             rotation_angle = m_rotation_angle(m);
%             save(filepath_rotation, 'rotation_angle'); 
%         end
%         J = imrotate(J,-m_rotation_angle(m));
%         clear command exifout k exifrotation
        
        % Write tif image
        %if  ~exist(filepath_tif,'file');
            imwrite(uint16(65536*J),filepath_tif,'tif');
        %end
        
        % If the jpg file does not exist, begin
        %if  ~exist(filepath_jpg,'file');
            % Write jpg
            Jjpg = imresize(J,.4);
            %iminfo = imfinfo(m_wavelength_filepath{m}{w});
            %bitdepth = iminfo.BitsPerSample(1);
            %maxval = 2^bitdepth;
            Jjpg = uint8(256*Jjpg);
            imwrite(Jjpg,filepath_jpg,'jpeg','Quality', 50);
        %end
        
        clear Jjpg J I
    end
    
    % Save auxiliary info
    filepath_shutter_speed= sprintf('%s%s_%s_shutter_speed.mat',subpath_matlab_dir{m},...
        m_mss{m},m_folio{m});
    filepath_aperture= sprintf('%s%s_%s_aperture.mat',subpath_matlab_dir{m},...
        m_mss{m},m_folio{m});
    filepath_stretchval= sprintf('%s%s_%s_stretchval.mat',subpath_matlab_dir{m},...
        m_mss{m},m_folio{m});
    if ~exist(filepath_shutter_speed, 'file');
        save(filepath_shutter_speed,'w_shutter_speed');
    end
    if ~exist(filepath_aperture, 'file');
        save(filepath_aperture,'w_aperture');
    end
    if ~exist(filepath_stretchval, 'file');
        save(filepath_stretchval,'w_stretchval');
    end
    
    cell_summary = cell(n_w, 4);
    cell_summary(:,1) = w_wavelength;
    cell_summary(:,2) = cellstr(num2str(w_aperture));
    cell_summary(:,3) = cellstr(num2str(w_shutter_speed));
    cell_summary(:,4) = cellstr(num2str(w_stretchval));
    
    %DJK if ~exist(filepath_summary_txt, 'file');
    %DJK     dlmcell(filepath_summary_txt,cell_summary);
    %DJK end
end

clear m w I
%fprintf('Completed %s successfully\n',m_path_upper{m}(slashindex(end)+1:end));

end




