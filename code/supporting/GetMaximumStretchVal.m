% Get maximum values 
info_slash = '/';

% Get folders 
path_source_previous = '/Volumes/Scythica/Processed/';
[m_path_upper] = uipickfiles('filterspec',path_source_previous,'Prompt','Please choose folders to process');
m_path_upper = m_path_upper';
n_m = numel(m_path_upper); 
for m = 1:n_m;
    m_path_upper{m} = sprintf('%s%s',m_path_upper{m},info_slash);
end
fprintf('\n***********************************************************\n');
fprintf('Folios to process: \n');



% Get stretchVal

for m = 1:n_m;
    if m == 1;
    nbands = numel(w_stretchval);
    M_stretchval = zeros(n_m,nbands);
    end
    ix_slash = strfind(m_path_upper{m}, info_slash);
    name = m_path_upper{m}(ix_slash(end-1)+1:ix_slash(end)-1);
    filepath_stretchVal = sprintf('%s%s+matlab%s%s_stretchval.mat', m_path_upper{m},name,info_slash, name);
    load(filepath_stretchVal);
    M_stretchval(m,:) = w_stretchval;
   % m_stretch_val = sprintf('%s%s%s'
end

stretchValMax = max(M_stretchval,[],1);

%foo = [M_stretchval; stretchValMax];
filepath_stretchValMax = sprintf('%sstetchValMax.mat',path_source_previous);
save(filepath_stretchValMax, 'stretchValMax');

foo = 1
