function list = get_list_csv(filename)
%reads list of files from filename.csv and returns cell array with list
%names

file_id = fopen(filename , 'r');
rows = textscan(file_id , '%s' , 'Delimiter' , '\n');
rows = rows{1};
fclose(file_id);

for i=1:length(rows)
    r = textscan(rows{i} , '%s' , 'Delimiter' , ',');
    rows{i} = r{1}';
end
rowlens = cellfun(@length , rows);
cols = max(rowlens);
list = cell(length(rows) , cols);

for i=1:length(rows)
    [list{i , 1:rowlens(i)}] = rows{i}{:};
end
