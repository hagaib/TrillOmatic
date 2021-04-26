function [result_band , tpeaks] = area_of_interest4( audiofilename , isplot , x , fs)
%GET_AOA Get Area of Interest of imput track 'audiofilename'
% freqband = is a 2x1 array with low and high frequencies in band
% This version has envelope cross correlation implemented


%% Constants and parameters
maxFreqsCount = 10;
minfreq = 0.5; %kHz
minBandwidth = 0.15; % kHz
additiveMargin = 0.3; % kHz
gaussFitCount = 3;
modelValidatePercent = 0.9;

%% Prepare Spectrogram
if(nargin < 4)
    [x, fs] = audioread(audiofilename);
end
% soundsc(x,fs);
if(size(x,2)>1)
    x = x(:,1);
end
time = (0:length(x)-1) / fs;
winsize = 256;
zeropad = 4;
overlap = 0.75;
[S,F,T] = spectrogram(x,hanning(winsize),floor(winsize*overlap),...
                    winsize*zeropad, fs, 'yaxis');


%% Smooth image with 2d median filtering
origspect = 20*log10(abs(S));
% medfiltspect = medfilt2(origspect , [5 5]);
medfiltspect = origspect;
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
noisefloor = medfilt1(fhist , floor(winsize*zeropad/2 * 0.5));
if(isplot)
plot(F ,fhist) , plot(F , noisefloor) , hold off
title('Frequency Histogram');
end
% fhist = medfilt1(fhist , 3);

fhist = max(fhist - noisefloor , 0);

% figure(fignum); fignum=fignum+1;
% plot((0:length(x)-1)*fs/length(x) , abs(fft(x)));
Fvalidate = zeros(size(F));
mu = [];
[pks, locs] = findpeaks(fhist, 'SortStr' , 'descend');
minPeakDist = 500;%1000;
minPeakDistBins = find(F>minPeakDist , 1);
while(true)
    %% Gaussian Model Fitting
    
%     mu = zeros(gaussFitCount*100 , 1);
    F_remain = F(~Fvalidate);
    b_locs_iter = ~Fvalidate(locs);
    
    locs_iter = locs(~Fvalidate(locs));
    pks_iter = pks(~Fvalidate(locs));
    
    if(length(pks_iter)<1) , break; end
    minPeakHeight = pks_iter(1)*0.05;%10;YZ
    locs_iter = locs_iter(pks_iter > minPeakHeight);
    
    i=length(mu);
    while true
        inRange = abs(locs_iter - locs_iter(1)) < minPeakDistBins;
        outRange = ~inRange;
    %     nextindex = find(outRange,1);
        muloc = floor(mean(locs_iter(inRange)));
        mu = [mu ; F(muloc)];
%         if(sum(outRange)==0 || i>2) , i=i-1; break; end
        if(sum(outRange)==0 || length(mu)>=gaussFitCount) , break; end
        locs_iter = locs_iter(outRange);
        i=i+1;
    end
    gaussFitCount = max(length(mu) , gaussFitCount);
    while(length(mu) < gaussFitCount)
        %mu(find(mu==0,1)) = min(mu(mu~=0)) - minPeakDistBins;
        if(min(mu) >= minPeakDist)
            muval = floor(0.5*min(mu(mu>0)));%YZ
        else
            muval = floor(2*max(mu(mu>0)));
        end
        mu = [mu ; muval];
    end
    while(length(mu) > gaussFitCount)
        mu = mu(mu ~= min(mu)); 
    end

    if(isplot)
    fun = @(x) gaussian_mix_residual(x , F, fhist , fignum);
    else
    fun = @(x) gaussian_mix_residual(x , F, fhist);
    end
    % x0 =   [mu(1) 1000, mu(2) 1000]';
    x0 = 1000 * ones(1 , gaussFitCount*2);
    x0(1:2:2*gaussFitCount) = mu;
    constraints_low = 100 * ones(1 , gaussFitCount*2);
    constraints_low(1:2:2*gaussFitCount) = 500; % [500 , 200 , 500 , 200]
    constraints_high = 400*ones(1 , gaussFitCount*2);
    constraints_high(1:2:2*gaussFitCount) = inf;
    % options = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective');
    % [gaussparams,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,constraints_low,constraints_high,options);
    [gaussparams,resnorm,fval_1,~,extra_arguments] = LevenbergMarquardt(fun, x0, constraints_low, constraints_high);

    % Out of Bounds fix (inspired by white noise)
    global coeffs;
    gaussianHeightThresh = 1;
    sigmaThresh = 4*10^3;
    sigma = gaussparams(2:2:end);
    inBounds = coeffs > gaussianHeightThresh & sigma' < sigmaThresh;
    gaussFitCount = sum(inBounds);
    inBounds = repmat(inBounds' , 2 , 1);
    gaussparams = gaussparams(inBounds(:));


    %% Find Frequency Bands of Interest
    freq_bands = zeros(gaussFitCount, 2);
    freq_bands_large = zeros(gaussFitCount, 2);
    for i=1:gaussFitCount
        mumu = gaussparams(2*i-1);
        sigma = abs(gaussparams(2*i));
        freq_bands(i,:) = mumu + 2.5*sigma*[-1,1]; %freq_bands(i,:) = mumu + 2*sigma*[-1,1];
        freq_bands_large(i,:) = mumu + 3*sigma*[-1,1];
    end
    freq_bands = freq_bands(freq_bands(:,2)>minfreq*10^3 , :);
    freq_bands_large = freq_bands_large(freq_bands_large(:,2)>minfreq*10^3 , :);
    freq_bands(:,1) = max(freq_bands(:,1) , minfreq*10^3);
    freq_bands_large(:,1) = max(freq_bands_large(:,1) , minfreq*10^3);

    %% Model validation to data
    for i=1:size(freq_bands , 1)
        Fvalidate = Fvalidate | ( F > freq_bands(i , 1) & F < freq_bands(i , 2) );
    end
    
    %     if(isplot) , fprintf('%f data fitted\n', norm(residual , 1) / norm(fhist , 1)); end
    if(isplot) , fprintf('%f data fitted\n', norm(fhist(Fvalidate) ,1) / norm(fhist , 1)); end


    % Loop termination condition
%     if(norm(residual , 1) / norm(fhist , 1) < 1 - modelValidatePercent  || gaussFitCount == 4) ,  break;
    if(norm(fhist(Fvalidate) , 1) / norm(fhist , 1) > modelValidatePercent  || gaussFitCount == 4) ,  break;
    else gaussFitCount = gaussFitCount + 1; end
end


if(isplot) , freq_bands , fignum=fignum+1; end

% 
% %% Minimum Bandwidth Thresholding
% bwapprove_bool = (freq_bands(:,2)-freq_bands(:,1)) > minBandwidth*10^3;
% freq_bands = freq_bands(bwapprove_bool,:);
% maxAreasOfInterest = length(find(bwapprove_bool));
% 
%% Filter Frequency Bands of Interest
band_count = size(freq_bands, 1);
xfilt = zeros(length(x) , band_count);
for i=1:band_count
    wpass = freq_bands(i,:) /fs *2;
    wstop = [wpass(1)*0.98 , wpass(2)*1.02];
    ftype = 'bandpass';
    if(wpass(2) >=1 || wstop(2) >=1 ) %high pass
        wpass = wpass(1); wstop = wstop(1);
        ftype = 'high';
    end
    Rp = 0.05; Rs = 40;
    n = cheb2ord(wpass , wstop , Rp , Rs);
    [z,p,g] = cheby2(n,Rs , wstop , ftype);
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

%% ENERGY CALCULATION
%% TEO Calculation

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
xteo = filter(g,1,xteo(end:-1:1 , :) );
xteo = xteo(end:-1:1 , :);
sigmabins = 1/((L-1)/5*sqrt(2))/(2*pi)*fs; %(new sigma is 1/(sqrt(2)*sigma)
economic_samplerate = floor(sigmabins*5*2^4);
subsample_factor = floor(fs / economic_samplerate);

if(isplot)
figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xteo(:,i));
    title(sprintf('Lowpass TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end
end

% %% OR: Short Time Energy Calculation
% L = floor(fs*0.02);
% if(mod(L,2)==0) , L=L+1; end
% g = gausswin(L);
% stenergy = filter(g , 1 , xfilt.^2);
% xteo = stenergy;



%% Energy Elimination
e_discard_threshold = -13; %db
e_overlap_threshold = -0.05; % db
overlap_threshold = 0.8; %percent
bw_energy = max(sum(xteo-repmat(median(xteo) , size(xteo,1),1)) , 0)';% ./ ( freq_bands(:,2) - freq_bands(:,1) );
[~ , m1] = max(bw_energy); [~ , m2] = min(bw_energy);
e_diff = 10*log10(bw_energy / bw_energy(m1));
if(isplot) , fprintf('Energy differece is %f db\n' , e_diff); end
highest_e_bands = e_diff > e_overlap_threshold;
highest_freqbands = freq_bands(highest_e_bands , :);
overlap_matrix = zeros(sum(highest_e_bands) , sum(highest_e_bands));
intervals_toremove = false(sum(highest_e_bands) , 1);
for i = 1:size(overlap_matrix , 1)
    for j = i+1:size(overlap_matrix , 2)
        interval1 = highest_freqbands(i,:);
        interval2 = highest_freqbands(j,:);
        inter_inter = interval_intersect(interval1 , interval2);
        overlap_matrix(i,j) = inter_inter(2) - inter_inter(1);
        overlap_matrix(j,i) = overlap_matrix(i,j);
        overlap_matrix(i,j) = overlap_matrix(i,j) / (interval1(2) - interval1(1));
        overlap_matrix(j,i) = overlap_matrix(j,i) / (interval2(2) - interval2(1));
        
        if(overlap_matrix(i,j) > overlap_threshold)
            fprintf('Interval [%0.1f , %0.1f] is contained in interval [%0.1f , %0.1f]!\n' , ...
                interval1(1) , interval1(2) , interval2(1) , interval2(2));
            intervals_toremove(j) = true;
        elseif(overlap_matrix(j,i) > overlap_threshold)
            fprintf('Interval [%0.1f , %0.1f] is contained in interval [%0.1f , %0.1f]!\n' , ...
                interval2(1) , interval2(2) , interval1(1) , interval1(2));
            intervals_toremove(i) = true;
        end
    end
end

highest_indices = find(highest_e_bands);
intervals_toremove = highest_indices(intervals_toremove);
intervals_remain = setdiff(1:size(freq_bands , 1) , intervals_toremove);
freq_bands = freq_bands(intervals_remain , :);


high_e_bands = floor(e_diff(intervals_remain)) > e_discard_threshold;
xteo = xteo(: , high_e_bands);
freq_bands = freq_bands(high_e_bands , :);
band_count = size(freq_bands , 1);

xteo_subsampled = xteo(1:subsample_factor:end,:);


% if(floor(e_diff) <= e_threshold)
%     band_count = band_count - 1;
%     xteo = xteo(:,1:end~=m2);
%     freq_bands = freq_bands(1:end~=m2 , :);
% end

%% Peak Rate Curve Fitting
if(isplot) , figure(fignum); fignum=fignum+1; end
l2err = nan(band_count,1);
l1err = nan(band_count,1);
linferr = nan(band_count,1);
periodmean = nan(band_count , 1);
peaklocsstd = nan(band_count , 1);
pdratio = nan(band_count , 1);
meanvariation = nan(band_count , 1);
peakloc_bands = cell(band_count , 1); peakloc_bands(:,1) = {nan};
for i=1:band_count
    %find peaks to calculate pulse rate
    [p , t] = findpeaks(xteo(:,i) , 'SortStr','descend');
    t = t/fs;
    t = sort(t(1:min(10 , floor(length(t)*0.7))));
    tx = (1:length(t))';
    dt = t(2:end)-t(1:end-1);
    [pol, S] = polyfit(tx,t,2);
    [pval,delta] = polyval(pol,tx,S);
    l2err(i) = sum((t-pval).^2) / length(t);
    l1err(i) = sum(abs(t-pval)) / length(t);
    linferr(i) = max(abs(t-pval)) / length(t);
    
%     if(isplot)
%     subplot(band_count,1,i)
%     plot(t,'b'), hold on; plot(pval,'r'); hold off;
%     title(sprintf('Quad Poly Fit for band: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
%     end


%     [periods , peaklocs] = env_xcorr2(x , fs, time, xteo(:,i) , 10^-1, isplot);
    %use economic version: 
    [periods , peaklocs] = env_xcorr3(x , economic_samplerate, time(1:subsample_factor:end), xteo_subsampled(:,i) , 10^-1, isplot);
    peakloc_bands{i} = peaklocs;
    peaklocs = floor(peaklocs*economic_samplerate);
    p = 40; phigh = 5;
    b_periods_main = periods > prctile(periods , p) & ...
                   periods < prctile(periods , 100-phigh);
    b_peaklocs_main = [b_periods_main; 0] | [0 ; b_periods_main]; % use only peaks with near neighbours
    periods_main = periods(b_periods_main);
    peaklocs_main = peaklocs(b_peaklocs_main);
    periodmean(i) = mean(periods_main);
    
    if(length(periods_main) < 2 && band_count > 1)
        periodmean(i) = inf;
        meanvariation(i) = 0;
        peaklocsstd(i) = inf;
        pdratio(i) = 0;
        continue;
    end
    
%     p=25;
%     peaklocs_main = peaklocs(peaklocs > prctile(peaklocs , p) & ...
%                             peaklocs< prctile(peaklocs , 100-p));
    ind_periods_main = find(b_periods_main);
    dipheights = nan(length(ind_periods_main) , 1);
    ind_dipheights = nan(length(ind_periods_main) , 1);
    heightDiffs = nan(length(dipheights)*2 , 1);

    for k=1:length(ind_periods_main)
            [dipheights(k) , ind_dipheights(k)] = min(xteo_subsampled(peaklocs(ind_periods_main(k)):peaklocs(ind_periods_main(k)+1) , i) );
            ind_dipheights(k) = ind_dipheights(k) + peaklocs(ind_periods_main(k));
            heightDiffs(2*k-1:2*k) = xteo_subsampled(peaklocs([ind_periods_main(k) , ind_periods_main(k)+1]) , i) - dipheights(k);
    end

    
    peakheights = xteo_subsampled(peaklocs_main , i);
%     peakheights = mean([peakheights(1:end-1) , peakheights(2:end)] ,2);
    
    pdratio(i) = mean(peakheights) / mean(dipheights);
%     heightDiffs = nan(size(peakheights));
%     heightDiffs(1:end-1) = peakheights(1:end-1) - dipheights;
%     heightDiffs(end) = peakheights(end) - dipheights(end);
    meanvariation(i) = mean(heightDiffs);
    peaklocsstd(i) = std(periods);
end
energy_thresh = 1;
periodmean_thresh = 10^-2;
periodstd_thresh = 15;
beststd_thresh = 20^10^-3;
meanvar_thresh = 0.1; %percents
[~ , perstd_sorted] = sort(peaklocsstd);
[~ , meanvar_sorted] = sort(meanvariation , 'descend');
minstd_band = perstd_sorted(1);
[~ ,pdratio_sorted] = sort(pdratio , 'descend');
[~ ,energy_sorted] = sort(e_diff , 'descend');
energy_good = abs(e_diff) < energy_thresh;
best_band = meanvar_sorted(1);
% while(periodstd(best_band) < beststd_thresh) 
% end
goodbandsmean = abs(periodmean - periodmean(best_band)) < periodmean_thresh | ...
            abs(periodmean - 2*periodmean(best_band)) < periodmean_thresh | ...
            abs(periodmean - periodmean(best_band)/2) < periodmean_thresh;
goodbandsstd = peaklocsstd/peaklocsstd(best_band) < periodstd_thresh;
goodbandsmeanvar = meanvariation/meanvariation(best_band) > meanvar_thresh;
goodbands = goodbandsmean & goodbandsstd & goodbandsmeanvar;
% [~,minband] = min(l2err);

goodbands = freq_bands(goodbands,:);
if(isempty(goodbands)) , result_band = [-inf , inf]; tpeaks = {nan}; else
    result_band = [min(goodbands(:,1)) max(goodbands(:,2))];
    result_band = result_band + additiveMargin * 10^3 * [-1 1];
    result_band(1) = max(result_band(1) , minfreq*10^3);
    tpeaks = peakloc_bands{best_band};
end
if(isplot)
l2err
% l1err
% linferr
fprintf('Band Width of Interest: %f,%f\n',result_band(1), result_band(2));
end