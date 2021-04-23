filenames =  get_list_csv('filenames.csv');

filename = filenames{1};
% filename =  '1605-00.44.57.983-lotr-4.wav';


[x , fs] = audioread(filename);
t = (0:(length(x)-1))/fs;

if(sum(size(x)==2))
    if(size(x,1)==2)
        x = transpose(x);
    end
    x = x(:,1);
end

% yin = yin_wrapper(x , fs);

f0=yin.f0;
dips = yin.dips;
time = yin.time;
perprobs = 0.2;
mindur = 0.025;

[periodictime , nonpertime , persegs , nonpersegs] = ...
    pertime(dips , time , mindur );
% 
%%
% plotting the stuff
figure(1)
p1 = subplot(2,1,1);
plot(t,x);
title('Periodic Time Example');

p2 = subplot(2,1,2);
spectrogram(x , 440 , 350 , 1024 , fs , 'yaxis')
h = plot_yin(f0/1000 , dips , time , p2);
colorbar off
hold on
plot(yin.time , periodictime , 'k' , 'Linewidth' , 2)
hold off
set(h(end) , 'Visible' , 'off')
title('L = 0.025')
linkaxes([p1,p2] , 'x');
