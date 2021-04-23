function [detect , segs , env , teo] = longtrill_syllable_detection41(x , fs , yin , filename , f0band)
% Algorithm 3 implementation

if(nargin < 4)
    filename = [];
end

if(nargin<5)
    f0band = [1800 3500]; % halcyon smyrnensis
end

t = ((0:length(x)-1)/fs)';


%% TEO estimation
n = cheb2ord(f0band/fs*2 , [f0band(1)*0.95 , f0band(2)*1.05]/fs*2 , 0.5 , 40);
[z,p,k] = cheby2(n , 40, f0band/fs*2 , 'bandpass');
sos = zp2sos(z,p,k);
xx0 = sosfilt(sos , x);

[~ , threshold_teo , teo] = trilltime_TEO(xx0 , fs);

%% Envelope Calculation
absfft = abs(fft(teo , length(teo)*2));%fft interpolation so real peak is not confused with harmonic multiplicities
[p , peakindex] = findpeaks(absfft(1:floor(length(absfft)/2)));
[~ , maxt] = max(p);
maxindex = peakindex(maxt);
PRI = 1/((maxindex-1)/length(absfft)*fs);

env = envelope(teo , floor(fs*PRI*0.8),'peak' );
%%

%% yinslopes estimation
wint = PRI*0.2; % alittle less than half! %wint = 0.02;
threshold_std = 1000;
yinslopvar = yinslopevarfilt(yin , wint);

L = floor(yin.steprate*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
yinslopvar = filter(g , 1 , yinslopvar);
yinslopvar = yinslopvar(:);
yinslopvar = [yinslopvar((L+1)/2:end) ; zeros((L-1)/2 , 1)];

%% Calculate VDF0 Threshold
margin_t = 0.1;
yinsteps = floor( margin_t / (yin.time(2) - yin.time(1) ) );
[dips1 , ~] = findpeaks(-yinslopvar(1:yinsteps));
[dips2 , ~] = findpeaks(-yinslopvar(end-yinsteps:end));
dips = [dips1;dips2];
vdf0_thresh = -prctile(dips , 80);
% Alternative approach: 1)calculate the mean / median at the margin 
% 2) claculate the max at the whole track
% 3) set threshold as the mean of both values


%% Analysis & Cleanup
[dips , locs]  = findpeaks(-sqrt(yinslopvar), ...
                    'MinPeakHeight' , -sqrt(vdf0_thresh));
locs = locs(mean(yin.dips([locs-1 , locs , locs+1]) , 2) < 0.9);
diptime = yin.time(locs)';
%%


%to delete:
olddiptime = diptime;
%% Envelope Rejects
localav = 0.012;
L = floor(localav*fs);
envrejects = zeros(size(diptime));
newenvrejects = zeros(size(diptime));
for i=1:length(diptime)
    dipsample = find(t>diptime(i) , 1);
    [maxte , imax] = max(teo(dipsample+(-L:L)));
%     meante = mean(teo(dipsamp+(-L:L)));
    if( maxte < env(dipsample)*0.7 )
        newenvrejects(i) = (dipsample-L+imax)/fs;
        envrejects(i) = diptime(i);
        diptime(i) = -1;
    else
        diptime(i) = (dipsample-L+imax)/fs;
    end
end
olddiptime = olddiptime(diptime>0);
diptime = diptime(diptime>0);
envrejects = envrejects(envrejects>0);
newenvrejects = newenvrejects(newenvrejects>0);

dipxxsamp = floor(diptime*fs);
diptime = diptime(teo(dipxxsamp)>threshold_teo);



% %% Detecting pulse interval
% detect = zeros(size(x));
% windur = 0.05;
% low = prcfilt(teo , floor(windur*fs) , 30 , 'omitnan')+threshold_teo;
% for i=1:length(diptime)
%     istart = floor(diptime(i)*fs);
%     iend = floor(diptime(i)*fs);
%     while(istart>0 && teo(istart)>low(istart)) , istart = istart-1; end
%     while(iend<length(x) && teo(iend)>low(iend)) , iend = iend+1; end
%         
%     detect(istart:iend) = 1;
% end
% 
% segs = logical2segments(detect , fs);

%% Segmentation
energy = xx0.^2;
periods = diptime(2:end) - diptime(1:end-1);
winlen = median(periods);
s_winlen = floor(winlen/2*fs); s_winlen = s_winlen + 1-mod(s_winlen , 2);
s_winlenShort = floor(s_winlen/4) + 1 - mod(floor(s_winlen/4) , 2);
hammingEnergy = filter(hamming(s_winlen)/sum(hamming(s_winlen)) , 1 , energy);
hammingEnergyShort = filter(hamming(s_winlenShort)/sum(hamming(s_winlenShort)) , 1 , energy);
hammingEnergySteadystate = hammingEnergy(ceil(s_winlen/2):end);
hammingEnergy = [hammingEnergy(ceil(s_winlen/2):end) ; zeros(floor(s_winlen/2) , 1)];
hammingEnergyShort = [hammingEnergyShort(ceil(s_winlenShort/2):end) ; zeros(floor(s_winlenShort/2) , 1)];
[detect , segs , plotparams] = ...
    longtrill_segmentation(x , fs , hammingEnergy, diptime , 'adaptive' , hammingEnergyShort);



%% Length Rejects
lenthresh = 0.001;
seglen = segs(2,:) - segs(1,:);
segs = segs(: , seglen>lenthresh);
detect = segments2logical(segs , t , fs);

%%
% plotting the stuff for ISCEE2018
if(~isempty(filename))
f = figure(1);
f.Position = [300 , 800 , 560 , 180];
% p1 = subplot(3,1,1);
plot(t,xx0 , 'b' , t, detect , 'm');
% title(filename);

f = figure(2);
f.Position = [300 , 800 , 560 , 180];
plot(t , teo , 'b' , t , env , 'm' , t , 0.7*env , 'c');
% plot(t , teo , 'b');
hold on
% line([0 , t(end)] , threshold_teo*[1 , 1] , 'Color', 'g');
% plot(olddiptime , teo(floor(olddiptime*fs)) , 'dg')
% plot(diptime , teo(floor(diptime*fs)) , 'db')
% plot(envrejects , teo(floor(envrejects*fs)) , 'd', 'Color' , 'g')
% plot(newenvrejects , teo(floor(newenvrejects*fs)) , 'd', 'Color' , 'b')en
plot(olddiptime , teo(floor(olddiptime*fs)) , '*r')
plot(diptime , teo(floor(diptime*fs)) , '*b')
plot(envrejects , teo(floor(envrejects*fs)) , 'o', 'Color' , 'r')
plot(newenvrejects , teo(floor(newenvrejects*fs)) , 'o', 'Color' , 'b')
hold off
legend({'TEO Energy' , 'Energy Envelope' , '0.7 Energy Envelope','Original' , 'Accepted' ,'Original' , 'Rejected'});
% legend({'TEO Energy' , 'Original Dips' , 'Improved Peaks'});
title('Energy Function')

sqrinslopvar = sqrt(yinslopvar);

f = figure(3);
f.Position = [300 , 800 , 560 , 180];
plot(yin.time , sqrinslopvar  , 'b')
hold on

line([0 , yin.time(end)] , sqrt(vdf0_thresh)*[1 ,1] , 'Color', 'g');
plot(yin.time(locs) , sqrinslopvar(locs) , 'dg')
% plot(olddiptime , sqrinslopvar(floor(olddiptime*yin.steprate)) , '*r')
% plot(envrejects , sqrinslopvar(floor(envrejects*yin.steprate)) , 'o' ,'Color' ,'r')
hold off
legend('DV(t)' , 'DV Threshold' , 'detected points')
title('Fundamental Derivative Variance');
% linkaxes([p1,p2,p3] , 'x');
end
% plotting the stuff
if(~isempty(filename))
figure(1)
p1 = subplot(3,1,1);
plot(t,xx0 , 'b' , t, detect , 'm');
title(filename);

p2 = subplot(3,1,2);
plot(t , teo , 'b' , t , env , 'm' , t , 0.7*env , 'c');
% plot(t , teo , 'b');
hold on
line([0 , t(end)] , threshold_teo*[1 , 1] , 'Color', 'g');
% plot(olddiptime , teo(floor(olddiptime*fs)) , 'dg')
% plot(diptime , teo(floor(diptime*fs)) , 'db')
% plot(envrejects , teo(floor(envrejects*fs)) , 'd', 'Color' , 'g')
% plot(newenvrejects , teo(floor(newenvrejects*fs)) , 'd', 'Color' , 'b')
plot(olddiptime , teo(floor(olddiptime*fs)) , '*r')
plot(diptime , teo(floor(diptime*fs)) , '*b')
plot(envrejects , teo(floor(envrejects*fs)) , 'o', 'Color' , 'r')
plot(newenvrejects , teo(floor(newenvrejects*fs)) , 'o', 'Color' , 'b')
hold off
legend({'TEO Energy' , 'Threshold TEO' , 'Energy Envelope' , '0.7 Energy Envelope','Original' , 'Accepted' ,'Original' , 'Rejected'});
% legend({'TEO Energy' , 'Original Dips' , 'Improved Peaks'});
title('Energy Function')

sqrinslopvar = sqrt(yinslopvar);

p3 = subplot(3,1,3);
plot(yin.time , sqrinslopvar  , 'b')
hold on

line([0 , yin.time(end)] ,sqrt(vdf0_thresh)*[1 ,1] , 'Color', 'g');
plot(yin.time(locs) , sqrinslopvar(locs) , 'dg')
% plot(olddiptime , sqrinslopvar(floor(olddiptime*yin.steprate)) , '*r')
% plot(envrejects , sqrinslopvar(floor(envrejects*yin.steprate)) , 'o' ,'Color' ,'r')
hold off
legend('DV(t)' , 'DV Threshold' , 'detected points')
title('Fundamental Derivative Variance');
linkaxes([p1,p2,p3] , 'x');
end
