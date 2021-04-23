%%
%trying to figure out if line fitting of pitch in syllable might help to
%eliminate false positives

filenames =  get_list_csv('filenames.csv');

filename = filenames{8};

[x , fs] = audioread(filename);
t = (0:(length(x)-1))/fs;

if(sum(size(x)==2))
    if(size(x,1)==2)
        x = transpose(x);
    end
    x = x(:,1);
end

yin = yin_wrapper(x , fs);
