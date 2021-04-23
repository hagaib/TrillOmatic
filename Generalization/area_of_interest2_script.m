clear
% close all
%% Constants and parameters
maxFreqsCount = 5;
minfreq = 0.5; %kHz
minBandwidth = 0.15; % kHz
additiveMargin = 0.2; % kHz
gaussFitCount = 2;
filename = ...
...'XC414678 - Whimbrel - Numenius phaeopus.wav';
...'XC193392 - Tawny Owl - Strix aluco.wav';
...'XC279231 - Varied Triller - Lalage leucomela yorki-2.wav';
...'XC35251 - Trilling Tailorbird - Orthotomus chloronotus chloronotus.wav';...
...'XC279350 - Chestnut-breasted Cuckoo - Cacomantis castaneiventris weiskei.wav';
'XC279350 - Chestnut-breasted Cuckoo - Cacomantis castaneiventris weiskei-2.wav';
...'XC337848 - Little Bronze Cuckoo - Chrysococcyx minutillus peninsularis.wav';

%% Prepare Spectrogram
[x, fs] = audioread(filename);
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
fignum = 1;figure(fignum); fignum=fignum+1;
p=subplot(2,1,1);
imagesc(T,F,origspect);
title('Spectrogram') , set(p,'YDir','Normal');
ylim([0 5*10^3]);
p=subplot(2,1,2);
imagesc(T,F,medfiltspect);
title('Median Filtered Spectrogram') , set(p,'YDir','Normal');
ylim([0 5*10^3]);



%% Creating Frequency Histogram
firstbin = find(minfreq*10^3 < F, 1);
fhist = zeros(size(F));
energyframe = nan(size(T));
for i=1:length(T)
    fftframe = medfiltspect(firstbin:end,i);
    energyframe(i) = sum(fftframe.^2);
    for j=1:maxFreqsCount
        [val,ind] = max(fftframe);
%         fprintf('max frequency is %f\n', f(firstbin+ind-1));
        fftframe(ind) = 0;
        fhist(firstbin+ind-1) = fhist(firstbin+ind-1) + 1;
    end
end


figure(fignum); fignum=fignum+1;
plot(F ,fhist) , hold on
fhist = medfilt1(fhist , 5);
plot(F ,fhist) , hold off
title('Frequency Histogram');



%% Gaussian Model Fitting
minPeakDist = 1000;
minPeakHeight = 10;
mu = zeros(gaussFitCount*100 , 1);
[pks, locs] = findpeaks(fhist, 'SortStr' , 'descend');
orig_locs = locs;
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
    mu(find(mu==0,1)) = min(mu(mu~=0)) - minPeakDistBins;
end
mu = F(mu);


fun = @(x) gaussian_mix_residual(x , F, fhist , fignum); fignum=fignum+1;
x0 =   [mu(1) 1000, mu(2) 1000]';
options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt');
[gaussparams,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);


%% Find Frequency Bands of Interest

freq_bands = zeros(gaussFitCount, 2);
freq_bands_large = zeros(gaussFitCount, 2);
for i=1:gaussFitCount
    mu = gaussparams(2*i-1);
    sigma = abs(gaussparams(2*i));
    freq_bands(i,:) = mu + 2*sigma*[-1,1];
    freq_bands_large(i,:) = mu + 3*sigma*[-1,1];
end

band_count = gaussFitCount;
bandA = freq_bands_large(1,:);
bandB = freq_bands_large(2,:);
inter = interval_intersect(bandA , bandB);
if(~isempty(inter) && ~isequal(inter, bandA) && ~isequal(inter, bandB))
    bandA = freq_bands(1,:);
    bandB = freq_bands(2,:);
    freq_bands = [freq_bands ; min(bandA(1) , bandB(1))  max(bandA(2) , bandB(2))];
    band_count = band_count+1;
end

freq_bands

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

figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xfilt(:,i));
    title(sprintf('%f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end


%% TEO Calculation

xteo = zeros(length(x) , band_count);
for i=1:band_count
    teo = teager(xfilt(:,i));
    xteo(:,i) = [teo(1) ; teo ; teo(end)];
end

figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xteo(:,i));
    title(sprintf('TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end

L = floor(fs*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
xteo = filter(g,1,xteo);

figure(fignum); fignum=fignum+1;
for i=1:band_count
    subplot(band_count,1,i)
    plot(time , xteo(:,i));
    title(sprintf('Lowpass TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end

%% Peak Rate Curve Fitting
figure(fignum); fignum=fignum+1;
l2err = nan(band_count,1);
l1err = nan(band_count,1);
linferr = nan(band_count,1);
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
    
    subplot(band_count,1,i)
    plot(t,'b'), hold on; plot(pval,'r'); hold off;
    title(sprintf('Quad Poly Fit for band: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
    
end

l2err
% l1err
% linferr
[~,minband] = min(l2err);
minband = freq_bands(minband,:);
result_band = minband + additiveMargin * 10^3 * [-1 1];
fprintf('Band Width of Interest: %f,%f\n',result_band(1), result_band(2));
    
%     
%     
% %     meddt = median(dt);
%     
% %     absfft = abs(fft(smoothenergy));
%    
% %     [p , t] = findpeaks(absfft(1:floor(length(absfft)/2)));
% %     [~ , maxt] = max(p);
% %     maxt = t(maxt);
%     PRI = 1/((maxt-1)/length(absfft)*fs)