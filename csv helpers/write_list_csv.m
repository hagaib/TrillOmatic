function write_list_csv(filename , cell_list , filemode)
%writes cell array with filenames into csv file
%filename - name of file to write
% cell_list - the data to write. cell array in string format 
% filemode- either 'w' or 'a'. Default is 'w'


%example: write_list_csv(filenamecsv , data1 , 'w');

if(nargin<3)
    filemode = 'w';
end


fileid = fopen(filename , filemode);
[col , row] = size(cell_list);
for i = 1:col
    for j=1:row
        fprintf(fileid , '%s' , cell_list{i,j});
        if(j<row) 
            fprintf(fileid , ','); 
        else
%             j
        end
    end
    fprintf(fileid , '\n');
end
fclose(fileid);