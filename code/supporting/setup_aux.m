function [  aux] = setup_aux( aux )
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
%Greek
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

m_path_upper = aux.m_path_upper;
m_name = aux.m_name;
n_m = aux.n_m;
options_delimiter_wavelength = aux.options_delimiter_wavelength;

m_wavelength_filepath =cell(n_m,1);
m_wavelength_file =cell(n_m,1);
m_wavelength =cell(n_m,1);
m_wavelength_file_new =cell(n_m,1);
for m = 1:n_m
    path_nospace = strrep(m_path_upper{m}, '\ ', ' ');
    path_nospace = path_nospace(1:end-1);
%     command = sprintf('cd '%s', m_path_upper{m});
%     [a,b] = system(command);
    cd(path_nospace);
    D = dir('*_F.tif');
    D = remove_hiddenfiles(D);
    n_w = numel(D);
    m_wavelength_filepath{m} = cell(n_w,1);
    m_wavelength_file{m} = cell(n_w,1);
    m_wavelength{m} = cell(n_w,1);
    m_wavelength_file_new{m} = cell(n_w,1);
    for w = 1:n_w
        ix_delimiter = strfind(D{w},options_delimiter_wavelength);
        m_wavelength_filepath{m}{w} = sprintf('%s%s',m_path_upper{m},D{w});
        m_wavelength_file{m}{w} = D{w};
        m_wavelength{m}{w} = D{w}(ix_delimiter:end);
        band_number = m_wavelength_file{m}{w}(end-8:end-6);
        wavelength = m_wavelength_file{m}{w}(end-16:end-10);
        m_wavelength_file_new{m}{w} = sprintf('%s+%s_%s.tif', m_name{m}, band_number, wavelength);
    end
end
aux.m_wavelength_filepath = m_wavelength_filepath;
aux.m_wavelength_file = m_wavelength_file;
aux.m_wavelength = m_wavelength;
aux.m_wavelength_file_new = m_wavelength_file_new;

end
    


%end




