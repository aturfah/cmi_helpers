function C = read_csv(file_location)
% Read data from csv
C = table2cell(readtable(file_location));
% Remove empty rows
to_del = [];
% Find the empty rows that we need to delete and delete them
for i = 1:size(C, 1)
    row = C(i,:);
    if strcmp(row{1},'')
        to_del = [i, to_del];
    end
end
for i = 1:size(to_del,2)
    C(to_del(i), :) = [];
end
end