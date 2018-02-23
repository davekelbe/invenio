function [ ] = resize_auxiliary_files( aux )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Preliminary setup
fprintf('\n***********************************************************\n');
fprintf('Resizing Auxiliary files: \n');

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
w_wavelength = aux.w_wavelength;
w_wavelength = aux.w_wavelength;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength_filepath = aux.m_wavelength_filepath;
%m_rotation_angle = aux.m_rotation_angle;

clear aux
%%
h = msgbox({'Add auxiliary files to image stack','  e.g., KTK_pseudo, MP_roi, etc.', 'Ensure consistent rotation'});
set(h, 'position', [600 500 250 100])
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 16 );
waitfor(h);
filepath_images = [];

for m = 1:n_m;
    fprintf('                 \t\t%s\n', m_name{m});

    cd(subpath_jpg_dir{m});
    
    %[filepath_images] = uipickfiles('Prompt','Please choose images files', 'FilterSpec', '*.jpg');
    
    counter = 1;
    
    D = dir('*KTK*');
    D = remove_hiddenfiles(D);
    for d = 1:numel(D);
        filepath_images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
        counter = counter + 1;
    end
    
    D = dir('*MP*');
    D = remove_hiddenfiles(D);
    for d = 1:numel(D);
        filepath_images{counter} = sprintf('%s%s',subpath_jpg_dir{m}, D{d});
        counter = counter + 1;
    end
    
    if isempty(filepath_images);
        n_f = 0;
    else
        n_f = size(filepath_images,2);
    end
    
    filepath_ref = sprintf('%s%s_stretch.jpg', subpath_jpg_dir{m},m_wavelength_file{m}{1}(1:end-4));
    inforef = imfinfo(filepath_ref);
    dimref = [inforef.Width inforef.Height];
    [~,ixsortref] = sort(dimref);
    dimrefsort = dimref(ixsortref);
    bitsref = inforef.BitDepth./inforef.NumberOfSamples;
    
    for b = 1:n_f;
        
        if ~isempty(strfind(filepath_images{b}, 'resize'));
            continue
        end
        
        info = imfinfo(filepath_images{b});
        dim = [info.Width info.Height];
        [~,ixsort] = sort(dim);
        dimsort = dim(ixsort);
        bits = info.BitDepth./info.NumberOfSamples;
        scale = dimrefsort(1)./dimsort(1);
        
        if scale ==1 && (bits==bitsref);
            filepath_target = sprintf('%s_resize.jpg', filepath_images{b}(1:end-4));
            command = sprintf('%s %s %s',  'mv', filepath_images{b}, filepath_target);
            system(command);
            continue
        end
        
        I = imread(filepath_images{b});
        Iscl = imresize(I,scale);
        
        if bits ~=8;
            Iscl = uint8(255*(double(Iscl)./(2^bits)));
        end
        
        filepath_target = sprintf('%s_resize.jpg', filepath_images{b}(1:end-4));
        imwrite(Iscl, filepath_target, 'jpg', 'quality', 50);
        
        command = sprintf('%s %s',  info_rmcall, filepath_images{b});
        system(command);
    end
    
end

end

