function [trilltime , diptime] = trilltime_estimate( x , fs , yin  , wint , threshold_std)
%TRILLTIME finds "exact" time of trill from a segment in which trill is
%known to be present
%uses yinslopes filter and TEO operator for initial estimation
%after initial estimation , elimination is performed according to gap variance

%output: 
%trilltime: estimated trill time interval
%diptime: dips of yinslopevarfilt suspected to indicate syllables

% example : trilltime = trilltime_estimate(x , fs , yin , 0.02 , 100);



%% yinslopes estimation
% wint = 0.02;
% threshold_std = 100;

yinslopvar = yinslopevarfilt(yin , wint);
 
%% TEO estimation
% [b,a] = butter(7 , [1800 3500]/fs*2 , 'bandpass');
% xx0 = filter(b, a , x);

[~ , threshold_teo , teo] = trilltime_TEO(x , fs);
%%


%%
%% Analysis
[dips , locs]  = findpeaks(-sqrt(yinslopvar), ...
                    'MinPeakHeight' , -threshold_std ,...
                    'MinPeakProminence' , sqrt(threshold_std));
diptime = yin.time(locs)';
yindips = yin.dips(locs)';


%% YIN Rejects
diptime = diptime(yindips < 0.3);


%% TEO Rejects
dipxxsamp = floor(diptime*fs);
diptime = diptime(teo(dipxxsamp)>threshold_teo);



%% PRI Estimate
absfft = abs(fft(teo));
[p , t] = findpeaks(absfft(1:floor(length(absfft)/2)));
[~ , maxt] = max(p);
maxt = t(maxt);
PRI = 1/((maxt-1)/length(absfft)*fs)


%% Variance Rejects
gaptime = diptime(2:end) - diptime(1:end-1);
[mg,ig] = max(gaptime);
while mg > PRI + 3*std(gaptime)
    center = find(diptime > mean(diptime) , 1);
    if(ig>center)
        ig = ig+1;
    end
    diptime(ig) = [];
    gaptime = diptime(2:end) - diptime(1:end-1);
    [mg,ig] = max(gaptime);
end

trilltime = [-1 , 1]*wint*2 + [diptime(1) , diptime(end)];
margin = 0.05;
trilltime(1) = max(trilltime(1) , margin);
trilltime(2) = min(trilltime(2) , length(x)/fs-margin);

