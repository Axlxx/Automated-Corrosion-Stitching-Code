% GetAllWildcardNamesFromFolder('Datasets/Realistic Scans/ ', '22July', 1)
function names = GetAllWildcardNamesFromFolder(folder, wildcard, extension, how_many)
%GetAllWildcardNamesFromFolder Get all the files from given folder that
%contain wildcard in the name
%   Detailed explanation goes here
    files = dir(fullfile(folder, ['*', wildcard, '*', extension]));
    if how_many > 0 && how_many < numel(files)
        files = files(1:1:how_many);
    end
    names = {files.name};
    names = fullfile(folder, names);
end 