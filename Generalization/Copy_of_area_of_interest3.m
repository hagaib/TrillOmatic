function result_band = area_of_interest3( audiofilename , isplot)
%GET_AOA GEt Area of Interest of imput track 'audiofilename'
% freqband = is a 2x1 array with low and high frequencies in band
% This version has envelope cross correlation implemented


%% Constants and parameters
maxFreqsCount = 5;
minfreq = 0.5; %kHz
minBandwidth = 0.15; % kHz
additiveMargin = 0.5; % kHz
gaussFitCount = 2;

%% Prepare Spectrogram
[x, fs] = audioread(audiofilename);
% soundsc(x,fs);
if(size(x,2)>1)
    x = x(:,1);
end
time = (0:length(x)-1) / fs;
winsize = 512;
zeropad = 4;
overlap = 0.75;
[S,F,T] = spectrogram(x,hanning(winsize),floor(winsize*overlap),...
                    winsize*zeropad, fs, 'yaxis');


%% Smooth image with 2d median filtering
origspect = 20*log10(abs(S));
medfiltspect = medfilt2(origspect , [5 5]);

%% Median Clipping
medclipspect = medfiltspect;
medcols = median(medclipspect);
medrows = median(medclipspect , 2);
counter = 0;
for i=1:size(medclipspect , 1)
    for j=1:size(medclipspect , 2)
        if(medclipspect(i,j) < 1*max(medrows(i) , medcols(j)))
            medclipspect(i,j) = -inf;
            counter = counter + 1;
        end
    end
end

if(isplot)
fignum = 1;figure(fignum); fignum=fignum+1;
p=subplot(3,1,1);
imagesc(T,F,origspect);
title('Spectrogram') , set(p,'YDir','Normal');
ylim([0 5*10^3]);
p=subplot(3,1,2);
imagesc(T,F,medfiltspect);
title('Median Filtered Spectrogram') , set(p,'YDir','Normal');
ylim([0 5*10^3]);
p=subplot(3,1,3);
imagesc(T,F,medclipspect);
title('Median Clipped Spectrogram') , set(p,'YDir','Normal');
ylim([0 5*10^3]);
end


%% Creating Frequency Histogram
firstbin = find(minfreq*10^3 < F, 1);
fhist = zeros(size(F));
energyframe = nan(size(T));
for i=1:length(T)
    fftframe = medclipspect(firstbin:end,i);
    energyframe(i) = sum(fftframe.^2);
    for j=1:maxFreqsCount
        [val,ind] = max(fftframe);
%         fprintf('max frequency is %f\n', f(firstbin+ind-1));
        fftframe(ind) = 0;
        fhist(firstbin+ind-1) = fhist(firstbin+ind-1) + 1;
    end
end

if(isplot)
figure(fignum); fignum=fignum+1;
plot(F ,fhist) , hold on
end
fhist = medfilt1(fhist , 5);
if(isplot)
plot(F ,fhist) , hold off
title('Frequency Histogram');
end


%% Gaussian Model Fitting
minPeakDist = 1000;
mu = zeros(gaussFitCount*100 , 1);
[pks, locs] = findpeaks(fhist, 'SortStr' , 'descend');
minPeakHeight = pks(1)*0.05;%10;YZ
locs = locs(pks > minPeakHeight);
minPeakDistBins = find(F>minPeakDist , 1);
i=1;
while true
    inRange = abs(locs - locs(1)) < minPeakDistBins;
    outRange = ~inRange;
%     nextindex = find(outRange,1);
    mu(i) = floor(mean(locs(inRange)));
    if(sum(outRange)==0 || i>2) , i=i-1; break; end
    locs = locs(outRange);
    i=i+1;
end
gaussFitCount = max(i , gaussFitCount);
mu = mu(1:gaussFitCount);
while(sum(mu==0))
    %mu(find(mu==0,1)) = min(mu(mu~=0)) - minPeakDistBins;
    mu(find(mu==0,1)) = floor(0.5*min(mu(mu~=0)));%YZ
end
mu = F(mu);

if(isplot)
fun = @(x) gaussian_mix_residual(x , F, fhist , fignum); fignum=fignum+1;
else
fun = @(x) gaussian_mix_residual(x , F, fhist);
end
x0 =   [mu(1) 1000, mu(2) 1000]';
% options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt');
% [gaussparams,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);
options = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective');
[gaussparams,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[500 , 200 , 500 , 200],[],options);


% Out of Bounds fix (inspired by white noise)
global coeffs;
gaussianHeightThresh = 1;
sigmaThresh = 4*10^3;
sigma = gaussparams(2:2:end);
inBounds = coeffs > gaussianHeightThresh & sigma < sigmaThresh;
gaussFitCount = sum(inBounds);
inBounds = repmat(inBounds' , 2 , 1);
gaussparams = gaussparams(inBounds(:));



%% Find Frequency Bands of Interest
freq_bands = zeros(gaussFitCount, 2);
freq_bands_large = zeros(gaussFitCount, 2);
for i=1:gaussFitCount
    mu = gaussparams(2*i-1);
    sigma = abs(gaussparams(2*i));
    freq_bands(i,:) = mu + 2*sigma*[-1,1];
    freq_bands_large(i,:) = mu + 3*sigma*[-1,1];
end
freq_bands = freq_bands(freq_bands(:,2)>minfreq*10^3 , :);
freq_bands_large = freq_bands_large(freq_bands_large(:,2)>minfreq*10^3 , :);
freq_bands(:,1) = max(freq_bands(:,1) , minfreq*10^3);
freq_bands_large(:,1) = max(freq_bands_large(:,1) , minfreq*10^3);

band_count = size(freq_bands , 1);
if(band_count > 1)
bandA = freq_bands_large(1,:);
bandB = freq_bands_large(2,:);
inter = interval_intersect(bandA , bandB);
if(~isempty(inter) && ~isequal(inter, bandA) && ~isequal(inter, bandB))
    bandA = freq_bands(1,:);
    bandB = freq_bands(2,:);
    freq_bands = [freq_bands ; min(bandA(1) , bandB(1))  max(bandA(2) , bandB(2))];
    band_count = band_count+1;
end
end
if(isplot) , freq_bands , end

% 
% %% Minimum Bandwidth Thresholding
% bwapprove_bool = (freq_bands(:,2)-freq_bands(:,1)) > minBandwidth*10^3;
% freq_bands = freq_bands(bwapprove_bool,:);
% maxAreasOfInterest = length(find(bwapprove_bool));
% 
%% Filter Frequency Bands of Interest
xfilt = zeros(length(x) , band_count);
for i=1:band_count
    wpass = freq_bands(i,:) /fs *2;
    wstop = [wpass(1)*0.98 , wpass(2)*1.02];
    Rp = 0.05; Rs = 30;
    n = cheb2ord(wpass , wstop , Rp , Rs);
    [z,p,g] = cheby2(n,Rs , wstop);
    sos = zp2sos(z,p,g);
%     fvt = fvtool(sos,'Fs',fs);
    xfilt(:,i) = sosfilt(sos, x);
end

if(isplot)
figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xfilt(:,i));
    title(sprintf('%f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end
end

%% TEO Calculation and Hilbert envelope

xteo = zeros(length(x) , band_count);
for i=1:band_count
    teo = teager(xfilt(:,i));
    xteo(:,i) = [teo(1) ; teo ; teo(end)];
end

if(isplot)
figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xteo(:,i));
    title(sprintf('TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end
end

L = floor(fs*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
xteo = filter(g,1,xteo);

if(isplot)
figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xteo(:,i));
    title(sprintf('Lowpass TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end
end

henv = abs(hilbert(xfilt));

%% Energy Elimination
e_threshold = -13; %db
bw_energy = sum(xteo)';% ./ ( freq_bands(:,2) - freq_bands(:,1) );
[~ , m1] = max(bw_energy); [~ , m2] = min(bw_energy);
e_diff = 10*log10(bw_energy(m2) / bw_energy(m1));
if(isplot) , fprintf('Energy differece is %f db\n' , e_diff); end
if(floor(e_diff) <= e_threshold)
    band_count = band_count - 1;
    xteo = xteo(:,1:end~=m2);
    henv = henv(:,1:end~=m2);
    freq_bands = freq_bands(1:end~=m2 , :);
end

%% Peak Rate Curve Fitting
if(isplot) , figure(fignum); fignum=fignum+1;end
l2err = nan(band_count,1);
l1err = nan(band_count,1);
linferr = nan(band_count,1);
% for i=1:band_count
%     %find peaks to calculate pulse rate
%     [p , t] = findpeaks(xteo(:,i) , 'SortStr','descend');
%     t = t/fs;
%     t = sort(t(1:min(10 , floor(length(t)*0.7))));
%     tx = (1:length(t))';
%     dt = t(2:end)-t(1:end-1);
%     [pol, S] = polyfit(tx,t,2);
%     [pval,delta] = polyval(pol,tx,S);
%     l2err(i) = sum((t-pval).^2) / length(t);
%     l1err(i) = sum(abs(t-pval)) / length(t);
%     linferr(i) = max(abs(t-pval)) / length(t);
%     
%     if(isplot)
%     subplot(band_count,1,i)
%     plot(t,'b'), hold on; plot(pval,'r'); hold off;
%     title(sprintf('Quad Poly Fit for band: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
%     end
% end
% [~,minband] = min(l2err);


%% periodicity estimation using envelope cross correlation
periodmean = nan(band_count , 1);
periodstd = nan(band_count , 1);

for i=1:band_count
    periods = env_xcorr(x , fs, time, xteo(:,i));
    periodmean(i) = mean(periods);
    periodstd(i) = std(periods);
end

[~ , minband] = min(periodstd);


minband = freq_bands(minband,:);
result_band = minband + additiveMargin * 10^3 * [-1 1];
if(isplot)
l2err
% l1err
% linferr
fprintf('Band Width of Interest: %f,%f\n',result_band(1), result_band(2));
end