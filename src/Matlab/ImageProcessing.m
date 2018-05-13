function [  ] = ImageProcessing(  )
%IMAGEPROCESSING Wrapper for image processing steps
%
%   There is no input to this function. Typing ImageProcessing in the
%   command line brings up a series of user interfaces which allow the user
%   to select file (directories) for processing. It is recommended that
%   the user change the source code directly to adjust default paths
%
%
% Image Processing   Tool
% Dave Kelbe <dave.kelbe@gmail.com>
% Rochester Institute of Technology
% Created for Early Manuscripts Electronic Library
% Sinai Pailimpsests Project
%
% V0.0 - Initial Version - February 6 2015
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
%% Setup  
aux.info_user = 'generic';

[aux] = processing_setup(aux);

%% Create reflectance tiffs 

[aux] = reflectance_tiffs11(aux);

%% Create truecolor images 

reflectance_tiffs9_rgb(aux);

%% Resize auxiliary files 

%resize_auxiliary_files(aux)

%% Select ROI for statistics 

gui_batch_select_ROI_wrapper2(aux)
 
%% Select band subsets 

choose_cube_bands_aux(aux);
%foo = 1;

end

