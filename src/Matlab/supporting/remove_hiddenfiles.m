function [ C ] = remove_hiddenfiles( D )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_D = numel(D);
is_valid = true(n_D,1);

for d = 1:n_D;
    if strcmp(D(d).name(1), '.');
        is_valid(d) = false;
    end
end

n_C = sum(is_valid);
C = cell(n_C,1);
counter = 1;

for d = 1:n_D;
    if ~is_valid(d); continue; end
    C{counter} = D(d).name;
    counter = counter + 1;
end


%enviwrite_bandnames(J,output_filename, bandnames);


end

