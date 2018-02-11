function [ ] = choose_band_subset( aux )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Preliminary setup
fprintf('\n***********************************************************\n');
fprintf('Load band subset: \n');

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
%w_wavelength = aux.w_wavelength;
%w_wavelength = aux.w_wavelength;
m_wavelength_file = aux.m_wavelength_file;
m_wavelength_filepath = aux.m_wavelength_filepath;
%m_rotation_angle = aux.m_rotation_angle;

clear aux
%%
for m = 1:n_m;
    fprintf('                 \t\t%s\n', m_name{m});
    
    filepath_cube = sprintf('%s%s_cube.txt', subpath_matlab_dir{m},m_name{m});
    fid = fopen(filepath_cube, 'w+');
    
    
    fprintf(fid,'Header.ncubes = 8\n');
    fprintf(fid,'\n');
    
      fprintf(fid,'Cube01\n');
    fprintf(fid,'MB365UV_001_F\n');
    fprintf(fid,'MB455RB_002_F\n');
    fprintf(fid,'MB470LB_003_F\n');
    fprintf(fid,'MB505Cy_004_F\n');
    fprintf(fid,'MB535Gr_005_F\n');
    fprintf(fid,'MB570Am_006_F\n');
    fprintf(fid,'MB625Rd_007_F\n');
    fprintf(fid,'MB700IR_008_F\n');
    fprintf(fid,'MB735IR_009_F\n');
    fprintf(fid,'MB780IR_010_F\n');
    fprintf(fid,'MB870IR_011_F\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVB47_021_F\n');
    fprintf(fid,'WBUVG58_020_F\n');
    fprintf(fid,'WBUVO22_019_F\n');
    fprintf(fid,'WBUVUVP_022_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'Cube02\n');
    fprintf(fid,'MB365UV_001_F\n');
    fprintf(fid,'MB455RB_002_F\n');
    fprintf(fid,'MB470LB_003_F\n');
    fprintf(fid,'MB505Cy_004_F\n');
    fprintf(fid,'MB535Gr_005_F\n');
    fprintf(fid,'MB570Am_006_F\n');
    fprintf(fid,'MB625Rd_007_F\n');
    fprintf(fid,'MB700IR_008_F\n');
    fprintf(fid,'MB735IR_009_F\n');
    fprintf(fid,'MB780IR_010_F\n');
    fprintf(fid,'MB870IR_011_F\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVB47_021_F\n');
    fprintf(fid,'WBUVG58_020_F\n');
    fprintf(fid,'WBUVO22_019_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'Cube03\n');
    fprintf(fid,'MB455RB_002_F\n');
    fprintf(fid,'MB470LB_003_F\n');
    fprintf(fid,'MB570Am_006_F\n');
    fprintf(fid,'MB870IR_011_F\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVB47_021_F\n');
    fprintf(fid,'WBUVG58_020_F\n');
    fprintf(fid,'WBUVO22_019_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'Cube04\n');
    fprintf(fid,'MB455RB_002_F\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVB47_021_F\n');
    fprintf(fid,'WBUVG58_020_F\n');
    fprintf(fid,'WBUVO22_019_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'\n');

    fprintf(fid,'Cube05\n');
    fprintf(fid,'MB455RB_002_F\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVB47_021_F\n');
    fprintf(fid,'WBUVG58_020_F\n');
    fprintf(fid,'WBUVO22_019_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'WBUVUVP_022_F\n');
    fprintf(fid,'TX940IR_031_F\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'Cube06\n');
    fprintf(fid,'TX450RB_026_F\n');
    fprintf(fid,'TX535Gr_027_F\n');
    fprintf(fid,'TX570Am_028_F\n');
    fprintf(fid,'TX635Rd_029_F\n');
    fprintf(fid,'TX700IR_030_F\n');
    fprintf(fid,'TX780IR_031_F\n');
    fprintf(fid,'TX870IR_032_F\n');
    fprintf(fid,'TX940IR_033_F\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'Cube07\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBUVUVP_022_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'TX450RB_026_F\n');
    fprintf(fid,'TX535Gr_027_F\n');
    fprintf(fid,'TX570Am_028_F\n');
    fprintf(fid,'TX635Rd_029_F\n');
    fprintf(fid,'TX700IR_030_F\n');
    fprintf(fid,'TX780IR_031_F\n');
    fprintf(fid,'TX870IR_032_F\n');
    fprintf(fid,'TX940IR_033_F\n');
    fprintf(fid,'\n');

    fprintf(fid,'Cube08\n');
    fprintf(fid,'MB940IR_012_F\n');
    fprintf(fid,'WBRBG58_018_F\n');
    fprintf(fid,'WBRBO22_017_F\n');
    fprintf(fid,'WBUVUVP_022_F\n');
    fprintf(fid,'WBUVUVb_023_F\n');
    fprintf(fid,'TX450RB_026_F\n');
    fprintf(fid,'TX535Gr_027_F\n');
    fprintf(fid,'TX570Am_028_F\n');
    fprintf(fid,'TX635Rd_029_F\n');
    fprintf(fid,'TX700IR_030_F\n');
    fprintf(fid,'TX780IR_031_F\n');
    fprintf(fid,'TX870IR_032_F\n');
    fprintf(fid,'TX940IR_033_F\n');
    fprintf(fid,'\n');
end

end

