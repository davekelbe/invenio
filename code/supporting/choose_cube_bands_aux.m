function [  ] = choose_cube_bands_aux( aux )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Add to cube text file

m_path_upper = aux.m_path_upper;
m_folio = aux.m_folio;
m_mss = aux.m_mss;
m_name = aux.m_name;
%m_wavelength_file_new = aux.m_wavelength_file_new;
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
subpath_tiff_mask_dir = aux.path_tiff_mask_dir;
subpath_jpg_mask_dir = aux.path_jpg_mask_dir;
subpath_matlab_dir = aux.path_matlab_dir;
subpath_envi_dir = aux.path_envi_dir;
%w_wavelength = aux.w_wavelength;
%w_wavelength = aux.w_wavelength;
%m_wavelength_file = aux.m_wavelength_file;
%m_wavelength_filepath = aux.m_wavelength_filepath;
%rotation_angle = aux.m_rotation_angle;
info_colormap = aux.info_colormap;
m_wavelength_filepath = aux.m_wavelength_filepath;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength = aux.m_wavelength;
m_wavelength_file_new = aux.m_wavelength_file_new;

%% Setup

for m = 1:n_m;
    
    docontinue = true;
    w_filepath_band = {};
    while docontinue;
        prompt_str = sprintf('Please choose bands for %s',m_name{m});
        [w_filepath_band] = uipickfiles('filterspec',subpath_tiff_dir{m},'Prompt',prompt_str, 'Append', w_filepath_band);
        
        w_filepath_band = w_filepath_band';
        n_w = numel(w_filepath_band);
        
        if n_w < 2;
            docontinue = false;
            continue;
        end 
        
        % Determine manuscript name
        bandnames = cell(n_w,1);
        for w = 1:n_w;
            filepath_full = w_filepath_band{w};
            slash_ix = strfind(filepath_full,info_slash);
            name = filepath_full(slash_ix(end-2)+1:slash_ix(end-1)-1);
            ix_name = strfind(filepath_full,name);
            bandname = filepath_full(ix_name(end)+numel(name):end-4);
            bandnames{w} = bandname;
        end
        
        % cube file
        filepath_cube = sprintf('%s%s_cube.txt',subpath_matlab_dir{m}, m_name{m});
        
        if exist(filepath_cube, 'file');
            % cubeno = 0;
            fid = fopen(filepath_cube, 'r');
            text = textscan(fid,'%s');
            text = text{1};
            iscube = ~cellfun(@isempty,strfind(text,'Cube'));
            ix = find(iscube);
            cubeno = str2double(text{ix(end)}(5:6))+1;
            %    for l = 1:numel(iscube);
            %        if ~iscube(l);
            %            continue
            %        end
            %        cubeno_curr = str2double(text{l}(5:6))
            %        cubeno = max(cubeno,cubeno_curr);
            %    end
            fclose(fid);
        else
            cubeno = 1;
        end
        
        
        fid = fopen(filepath_cube, 'a+');
        fprintf(fid, '\r\n');
        fprintf(fid, 'Cube%02.0f\r\n', cubeno);
        for w = 1:n_w
            fprintf(fid, '%s\r\n', bandnames{w});
        end
        %    fprintf(fid, '\n');
        
    end
end


end

