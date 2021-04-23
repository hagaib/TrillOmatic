%% Actually , this script used to develop trill time algorithm. 
%% But now its used to find a 3rd method for trill syllable detection

filenames =  get_list_csv('filenames.csv');

filename = filenames{12};
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

% wint = 0.02;
% threshold_std = 100;
% trilltime = trilltime_estimate(x , fs , yin , 0.02 , 100);



%% yinslopes estimation
wint = 0.02;
threshold_std = 200;
yinslopvar = yinslopevarfilt(yin , wint);

L = floor(yin.steprate*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
yinslopvar = filter(g , 1 , yinslopvar);
yinslopvar = yinslopvar(:);
yinslopvar = [yinslopvar((L+1)/2:end) ; zeros((L-1)/2 , 1)];


%% TEO estimation
[b,a] = butter(7 , [1800 3500]/fs*2 , 'bandpass');
xx0 = filtfilt(b, a , x);

[~ , threshold_teo , teo , teo_orig] = trilltime_TEO(xx0 , fs);


%% Hamming Energy
hammeng = hamming_energy(xx0 , fs , 0.02 , 1);


%% Analysis & Cleanup
[dips , locs]  = findpeaks(-sqrt(yinslopvar), ...
                    'MinPeakHeight' , -threshold_std);
locs = locs(mean(yin.dips([locs-1 , locs , locs+1]) , 2) < 0.4);
diptime = yin.time(locs)';
dipxxsamp = floor(diptime*fs);
diptime = diptime(teo(dipxxsamp)>threshold_teo);
%%


%% TEO Area Calculation
teo_area = teo;
teo_area(teo_area < 20*threshold_teo) = 0;

this = 1;
area = 0;
for i=1:length(teo_area)
    if(teo_area(i)==0)
        teo_area(this:i-1) = area / (i-this);
        this = i+1;
        area = 0;
    end
    area=area+teo_area(i);
end


%% Envelope Calculation
absfft = abs(fft(teo));
[p , peakindex] = findpeaks(absfft(1:floor(length(absfft)/2)));
[~ , maxt] = max(p);
maxindex = peakindex(maxt);
PRI = 1/((maxindex-1)/length(absfft)*fs);

env = envelope(teo , floor(fs*PRI*0.8),'peak' );
%%


%% Envelope Rejects
localav = 0.012;
L = floor(localav*fs);
envrejects = zeros(size(diptime));
for i=1:length(diptime)
    dipsamp = find(t>diptime(i) , 1);
    maxte = max(teo(dipsamp+(-L:L)));
    meante = mean(teo(dipsamp+(-L:L)));
    if( maxte < env(dipsamp)*0.7 )
        envrejects(i) = diptime(i);
        diptime(i) = -1;
    end
end
diptime = diptime(diptime>0);
envrejects = envrejects(envrejects>0);


% %% todo : estimate mean gap by better methods (FFT on envelope or energy?)
% %% Variance Rejects
% gaptime = diptime(2:end) - diptime(1:end-1);
% [mg,ig] = max(gaptime);
% while mg > mean(gaptime) + 4*std(gaptime)
%     center = find(diptime > mean(diptime) , 1);
%     if(ig>center)
%         ig = ig+1;
%     end
%     diptime(ig) = [];
%     gaptime = diptime(2:end) - diptime(1:end-1);
%     [mg,ig] = max(gaptime);
% end
% %%



%%
% Total Trill Time Estimation
trilltime = [-1 , 1]*wint*2 + [diptime(1) , diptime(end)]
%%

% 
%%
% plotting the stuff
figure(1)
p1 = subplot(4,1,1);
plot(t,xx0);
title(filename);

p2 = subplot(4,1,2);
plot(t , teo , 'b' , t , teo_area , 'c' , t , env , 'm' );
hold on
line([0 , t(end)] , threshold_teo*[1 , 1] , 'Color', 'g');
plot(diptime , teo(floor(diptime*fs)) , '*')
plot(envrejects , teo(floor(envrejects*fs)) , 'o')
hold off
legend({'smooth teo'});

p3 = subplot(4,1,3);
plot(t , hammeng)

p4 = subplot(4,1,4);
plot(yin.time , sqrt(yinslopvar)  , 'b')
hold on
line([0 , yin.time(end)] , threshold_std*[1 ,1] , 'Color', 'g');
hold off
linkaxes([p1,p2,p3,p4] , 'x');
%%