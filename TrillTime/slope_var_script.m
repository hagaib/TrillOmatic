filenames =  get_list_csv('filenames.csv');

filename = filenames{14};
% filename =  '1605-00.44.57.983-lotr-4.wav';


[x , fs] = audioread(filename);
t = (0:(length(x)-1))/fs;

if(sum(size(x)==2))
    if(size(x,1)==2)
        x = transpose(x);
    end
    x = x(:,1);
end

yin = yin_wrapper(x , fs);

dyin = yin.f0(2:end) - yin.f0(1:end-1);
dyin = [0 , dyin];

% window length for dyin averaging
wint = 0.02;

fig = figure(1);
yinslopvar = yinslopevarfilt(yin , wint);
threshold = 100;
tt = yin.time(sqrt(yinslopvar) < threshold);
[tt(1) , tt(end)] + 0.5*[-wint , wint]
% 
%%
% plotting the stuff
figure(1)
p1 = subplot(3,1,1);
% plot(t,x);
spectrogram(x , 440 , 350 , 1024 , fs , 'yaxis')
colorbar off
plot_yin(yin.f0/1000 , yin.dips , yin.time , p1)
title(filename);

p2 = subplot(3,1,2);
plot(yin.time , dyin);%plot(yin.time , yin.f0);
title('f_0(t) Derivative - Forward Difference Approximation');

p3 = subplot(3,1,3);
plot(yin.time , sqrt(yinslopvar))
hold on
line([0 , t(end)] , [threshold , threshold] , 'Color', 'g');
hold off
legend('Variance Filter of f_0''(t)' , 'Threshold')
title('Variance Filter')
linkaxes([p1,p2,p3] , 'x');
