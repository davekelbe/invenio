function [ mypreferences ] = parse_EMEL_preferences( filepath_preferences )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(filepath_preferences);
    mypreferences = [];
    return
end
fid = fopen(filepath_preferences);
if fid <0;
    mypreferences = [];
    return
end

preferences = textscan(fid, '%s', 'Delimiter', '\n');
preferences = preferences{1};

is_empty = cellfun(@isempty,preferences);
preferences = preferences(~is_empty);

search = 'Creator';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
path_up = preferences(is);
mypreferences.Creator = path_up{1}(numel(search)+2:end);

search = 'processor';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
path_source = preferences(is);
mypreferences.initials = path_source{1}(numel(search)+2:end);

search = 'cubelist';
ix = strfind(preferences,search);
is = ~cellfun(@isempty, ix);
path_source = preferences(is);
temp = path_source{1}(numel(search)+2:end);
if strcmp(temp, 'true');
    mypreferences.cubelist = true;
else
    mypreferences.cubelist = false;
end

end

