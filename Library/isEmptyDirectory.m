% An empty directory contains no files nor subdirectories. 
% With Unix or Windows systems, every directory contains an entry for ?.? and almost every directory contains ?..? (except for a root directory); 
% an empty directory contains no other entries.


function flag_empty_dir = isEmptyDirectory(path)
if isdir(path)
    dir_content_list = dir(path); % list content of dir
    flag_empty_dir = length(dir_content_list) <= 2;  
        % returns 'true'/1 if the number of files (length(f)) is *greather than two*
        % returns 'false'/0 if the number of files (length(f)) is *less than or equal to two*
else
    error('Error: %s is not a directory', path);
end;
end