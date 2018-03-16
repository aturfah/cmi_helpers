% User selects a directory from their computer to load the files from.
% This should be the parent directory of the patient data.
function display_str = select_dir_fn()
directory = uigetdir;
if directory == 0
    display_str = 0;
    return
end
directory = strrep(directory, '\', '/');
%display_str = sprintf('Directory: %s', directory);
display_str = directory;
end