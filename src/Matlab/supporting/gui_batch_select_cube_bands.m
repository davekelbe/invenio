function [ ] = gui_batch_select_ROI_wrapper(  ) 
%UNTITLED Summary of this function goes here 
%   Detailed explanation goes here     
default_dir = '/Volumes/Philippi/Processed3/';
cd(default_dir);      
file.folders = uipickfiles('Prompt','Please choose folders');
 n_f = size(file.folders,2);
   
for f = 1:n_f;
    
   filepath_jpg = sprintf('%s/jpg/',file.folders{f});
cd(filepath_jpg);

%[file.images] = uipickfiles('Prompt','Please choose images files', 'FilterSpec', '*.jpg');

counter = 1;
D = dir('*KTK*');
for d = 1:numel(D);
    if D(d).name(1) == '1';
        continue
    end
    file.images{counter} = sprintf('%s%s',filepath_jpg, D(d).name);
    counter = counter + 1;
end

D = dir('*DJK_true*');
if numel(D)==2;
   D = D{2};
end
file.images{counter} = sprintf('%s%s',filepath_jpg, D.name);
counter = counter + 1;

D = dir('*WBUVO22_019_F*');
if numel(D)==2;
   D = D{2}; 
end
file.images{counter} = sprintf('%s%s',filepath_jpg, D.name);
counter = counter + 1;

D = dir('*WBUVG58_020_F*');
if numel(D)==2;
   D = D{2}; 
end
file.images{counter} = sprintf('%s%s',filepath_jpg, D.name);
counter = counter + 1;
 
D = dir('*WBUVB47_021_F*');
if numel(D)==2;
   D = D{2}; 
end
file.images{counter} = sprintf('%s%s',filepath_jpg, D.name);
counter = counter + 1;
 
D = dir('*TX940IR_031_F*'); 
if numel(D)==2; 
   D = D{2};
end
file.images{counter} = sprintf('%s%s',filepath_jpg, D.name);

n_f = size(file.images,2); 
  
counter = 1;
for b = 1:n_f; 
    B = imread(file.images{b});
    B = imresize(B,0.4); 
    [~,~,n_b] = size(B);
    I{counter} = B; 
    counter = counter + 1;
end      

bandnames = cell(n_b);
for b = 1:n_f;
    filepath_full = file.images{b};
    slash_ix = strfind(filepath_full,'/');
    filename = filepath_full(slash_ix(end)+1:end-4);
    bandnames{b} = filename;
end

target_dir = filepath_full(1:slash_ix(end));
gui_select_roi2(I,target_dir,bandnames) 
end

end

