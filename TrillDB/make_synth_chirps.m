names = get_list('filenames.csv');
new_names = cell(size(names));

for i=1:length(names)
    [xx,fs] = audioread(names{i});
    [t , A , f0] = get_chirp_parameters(xx , fs , 0.01 , 0.001, 0.3, 1500);
    params.t = t;
    params.A = A;
    params.f0 = f0;
    yy = generate_chirp(fs , f0 , t , A, [30 45]);
    filetype_index = strfind(names{i} , '.');
    filetype_index = filetype_index(end);
    name = [names{i}(1:filetype_index-1) , '-synth' , names{i}(filetype_index:end)];
    audiowrite(name , yy , fs);
    new_names{i} = name;
    name = [name(1:filetype_index+5) , '.mat'];
    save(name,'params');
end
write_list_csv('filesDB.csv' , new_names);