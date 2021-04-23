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

%%
% Get basic segmentation and data
harmonics = zeros(length(x) , 2);
[ detect , segs , env ] = longtrill_syllable_detection2(x , fs , yin , harmonics , [] , []);
params = extract_trill_parameters(x , fs , yin , segs , env);

%%
% Create data for polynomial fitting
syl_t = mean(params.syllable_timestamps ,1);
syl_f = params.syllable_f0med;
poly = polyfit(syl_t , syl_f , 1);
extrapol = polyval(poly , yin.time);

%%
%plotting stuff
%%
L=floor(0.01*fs);
figure(17)
p1 = subplot(2,1,1);
plot(t , x , 'b' , t ,detect*max(x) , 'c')
p2 = subplot(2,1,2);
spectrogram(x , hamming(L) , floor(L/8) , 2^10 , fs ,'yaxis')
colorbar off
plot_yin(yin.f0/1000 , yin.dips , yin.time , p2);
hold on
plot(yin.time , extrapol/1000 , 'k' , 'LineWidth' , 3)
ylim([0 10])
hold off
linkaxes([p1 , p2] , 'x')
