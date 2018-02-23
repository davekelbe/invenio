function [ ] = gui_batch_select_cube_bands_wrapper(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here 
default_dir = '/Volumes/'; 
cd(default_dir);
file.folders = uipickfiles('Prompt','Please choose folders');
n_f = size(file.folders,2);
info_slash = '/';

for f = 1:n_f;
    ix_slash = strfind(file.folders{f}, info_slash);
    name = file.folders{f}(ix_slash(end)+1:end);
    
    filepath_jpg = sprintf('%s/%s+jpg/',file.folders{f},name);
    cd(filepath_jpg);
    
    counter = 1;
    D = dir();
    for d = 1:numel(D);
        if D(d).name(1) == '.';
            continue
        end
       % if ~isempty(strfind(D(d).name,'KTK'));
       %     continue
       % end
        if ~isempty(strfind(D(d).name,'DJK'));
            command = sprintf('exiftool -colorcomponents %s',D(d).name);
            [~,exifout] = system(command);
            k = strfind(exifout, ':');
            n_bands = str2num(strtrim(exifout(k+2:end)));
            if n_bands>1;
                continue
            end
        end
        file.images{counter} = sprintf('%s%s',filepath_jpg, D(d).name);
        counter = counter + 1;
    end
    
    n_f = size(file.images,2);
    
    counter = 1;
    I = cell(n_f,1);
    for b = 1:n_f;
        B = imread(file.images{b});
        B = imresize(B,0.5);
        % [~,~,n_b] = size(B);
        I{counter} = B;
        counter = counter + 1;
    end
    
    bandnames = cell(n_f,1);
    for b = 1:n_f;
        filepath_full = file.images{b};
        slash_ix = strfind(filepath_full,'/');
        name = filepath_full(slash_ix(end-2)+1:slash_ix(end-1)-1);
        ix_name = strfind(filepath_full,name);
        bandname = filepath_full(ix_name(end)+numel(name):end-4);
        bandnames{b} = bandname;
    end
    
    
    %dir_upper = filepath_full(1:slash_ix(end-2));
    %target_dir = sprintf('%scube/',dir_upper);
    %
    target_dir = filepath_full(1:slash_ix(end));    
    cube = filepath_full(slash_ix(end-2)+1:slash_ix(end-1)-1);
    if ~exist(target_dir, 'dir');
        mkdir(target_dir);
    end
    gui_select_bands(I,target_dir,cube,bandnames)
end

end

