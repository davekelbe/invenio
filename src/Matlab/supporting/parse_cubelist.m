function [ cube ] = parse_cubelist( filepath_cubelist )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filepath_cubelist);

data = textscan(fid, '%s','headerlines', 0, 'collectoutput', true);
n_D = numel(data{1});

C = cell(n_D,1);
counter = 1;

for d = 1:n_D;
    C{counter} = char(data{1}(d));
    counter = counter + 1;
end
n_list = numel(C);

cube_ix = strfind(C,'Cube');
is_cube = ~cellfun(@isempty,cube_ix);
n_cubes = sum(is_cube);
cube_start = find(is_cube);
cube_end = cube_start-1;
cube_end = [cube_end(2:end); n_list];

cube = cell(n_cubes,1);
for c = 1:n_cubes;
    is_wavelength = false(n_list,1);
    is_wavelength(cube_start(c)+1:cube_end(c)) = true;
    cube{c} = C(is_wavelength);
end


end

