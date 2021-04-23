filenames =  get_list_csv('filenames.csv');

filename = filenames{14};

[x , fs] = audioread(filename);
t = (0:(length(x)-1))/fs;

if(sum(size(x)==2))
    if(size(x,1)==2)
        x = transpose(x);
    end
    x = x(:,1);
end

figure(4)
spectrogram(x,2048,1800,[],fs,'yaxis')


yin = yin_wrapper(x , fs);

%%
% Get basic segmentation and data
[detect , segs , env] = longtrill_syllable_detection41(x , fs , yin , filename);