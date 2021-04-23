filenames =  get_list_csv('filenames.csv');

filename = filenames{1};


[x , fs] = audioread(filename);
x = x(:,1);
N= length(x);
t = 1:N;
t = (t-1)/fs;

[bh,ah] = butter(7,1500/fs*2 , 'high');
xhigh = filtfilt(bh,ah,x);

[b,a] = butter(7, [1800 , 3500]/fs*2);
xpass = filtfilt(b,a,x);

%%
% plotting the stuff
figure(1)
p1 = subplot(2,1,1);
plot(t,x);
title(filename)
% title([filename , ' , high pass filtered']);

% 
% p3 = subplot(3,1,2);
% plot(t,xfilt)

p2 = subplot(2,1,2);
spectrogram(x , 440 , 350 , 1024 , fs , 'yaxis')
colorbar off
title('Spectrogram')



linkaxes([p1,p2] , 'x');
