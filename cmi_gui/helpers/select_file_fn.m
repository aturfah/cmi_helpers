function filename = select_file_fn()
[filename, pathname] = uigetfile({'*.csv';'*.xlsx';'*.xls'}, 'Select a file');
if filename == 0
    return
end
filename = strcat(pathname, filename);
end