clear
close all
%% Constants and parameters
maxFreqsCount = 5;
maxAreasOfInterest = 4;
minfreq = 0.5; %kHz
minBandwidth = 0.15; % kHz
noisefloorFrames = 10;
filename = ...
'XC193392 - Tawny Owl - Strix aluco.wav';
...'XC279231 - Varied Triller - Lalage leucomela yorki-2.wav';
...'XC35251 - Trilling Tailorbird - Orthotomus chloronotus chloronotus.wav';...
...'XC279350 - Chestnut-breasted Cuckoo - Cacomantis castaneiventris weiskei.wav';
...'XC279350 - Chestnut-breasted Cuckoo - Cacomantis castaneiventris weiskei-2.wav';
...'XC337848 - Little Bronze Cuckoo - Chrysococcyx minutillus peninsularis.wav';

%% Prepare Spectrogram
[x, fs] = audioread(filename);
% soundsc(x,fs);
if(size(x,2)>1)
    x = x(:,1);
end
time = (0:length(x)-1) / fs;
winsize = 512;
zeropad = 2;
overlap = 0.5;
[S,F,T] = spectrogram(x,hanning(winsize),floor(winsize*overlap),...
                    winsize*zeropad, fs, 'yaxis');
fignum = 1; figure(fignum); fignum=fignum+1;
spectrogram(x,hanning(winsize),floor(winsize*overlap),...
                    winsize*zeropad, fs, 'yaxis');
title(filename);

%% Creating Frequency Histogram
firstbin = find(minfreq*10^3 < F, 1);
fhist = zeros(size(F));
energyframe = nan(size(T));
for i=1:length(T)
    fftframe = abs(S(firstbin:end,i));
    energyframe(i) = sum(fftframe.^2);
    for j=1:maxFreqsCount
        [val,ind] = max(fftframe);
%         fprintf('max frequency is %f\n', f(firstbin+ind-1));
        fftframe(ind) = 0;
        fhist(firstbin+ind-1) = fhist(firstbin+ind-1) + 1;
    end
end

% fhist = medfilt1(fhist , 3);
figure(fignum); fignum=fignum+1;
plot(F ,fhist)
title('Frequency Histogram');

%% Noise Floor Histogram Estimate
nfhist = zeros(size(F));
energyframe(energyframe==0) = inf;
for i=1:noisefloorFrames
    [val,ind] = min(energyframe);
    energyframe(ind) = inf;
    fftframe = abs(S(firstbin:end,ind));
    for j=1:maxFreqsCount
        [val,ind] = max(fftframe);
%         fprintf('max frequency is %f\n', f(firstbin+ind-1));
        fftframe(ind) = 0;
        nfhist(firstbin+ind-1) = nfhist(firstbin+ind-1) + 1;
    end
end

figure(fignum); fignum=fignum+1;
hold on
plot(F ,nfhist , '.')
nfhist = medfilt1(nfhist , 3);
plot(F ,nfhist , '-')
hold off
title('Noise Floor Histogram');

%% Find Noise Bounds
inlobe = false;
lobearea = 0;
noiselobes = [];
for i=1:length(F)
    if(nfhist(i) > 0)
        lobearea = lobearea + nfhist(i);
        if(~inlobe) 
            inlobe = true;
            noiselobes = [noiselobes ; F(i) 0 0];
        end
    elseif(inlobe)
        noiselobes(end , 2:3) = [F(i-1) lobearea];
        inlobe = false; lobearea = true;
    end
end
% noiselobes is an nx3 matrix with columns: 1)Low F 2)High F 3) Area
noiselobes = noiselobes(noiselobes(:,3) > maxFreqsCount * noisefloorFrames * 0.5 , :);
[~ , ind] = max(noiselobes(:,3));
noise_bounds = noiselobes(ind , 1:2);

%% Find Peaks in Histogram
[pks, locs] = findpeaks(fhist, ...'MinPeakProminence' , max(fhist)*0.1, ...
                    'SortStr' , 'descend');
pks = pks(1:maxAreasOfInterest); locs = locs(1:maxAreasOfInterest);
[ peaks , ind_peaks , lobe_ampl ,lobe_bounds] = get_lobes(fhist);
sd = setdiff(locs,ind_peaks);
for el = sd'
    locs = locs(locs~=el);
end

locs_in_indpeaks = arrayfun(@(x)find(ind_peaks==x,1),locs);
freq_bands = lobe_bounds(locs_in_indpeaks,:);
freq_bands = F(freq_bands);


% %% Histogram Peaks and Noise Comparison
% lobe_noise_intersect = nan(maxAreasOfInterest,2);
% for i=1:maxAreasOfInterest
%     inter = interval_intersect(noise_bounds , freq_bands(i,:));
%     if(isempty(inter)) , inter = [nan , nan]; end
%     lobe_noise_intersect(i,:) = inter;
% end
% 
% %% Frequency Bands Acquisition
% for i=1:maxAreasOfInterest
%     for j=1:length(locs)
%         if(~isnan(lobe_noise_intersect(1)))
%             
%         end
%     end
% end




%% Minimum Bandwidth Thresholding
bwapprove_bool = (freq_bands(:,2)-freq_bands(:,1)) > minBandwidth*10^3;
freq_bands = freq_bands(bwapprove_bool,:);
maxAreasOfInterest = length(find(bwapprove_bool));

%% Filter Frequency Bands of Interest
xfilt = zeros(length(x) , maxAreasOfInterest);
for i=1:maxAreasOfInterest
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
for i=1:maxAreasOfInterest
    subplot(maxAreasOfInterest,1,i)
    plot(time , xfilt(:,i));
    title(sprintf('%f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end


%% TEO Calculation

xteo = zeros(length(x) , maxAreasOfInterest);
for i=1:maxAreasOfInterest
    teo = teager(xfilt(:,i));
    xteo(:,i) = [teo(1) ; teo ; teo(end)];
end

figure(fignum); fignum=fignum+1;
for i=1:maxAreasOfInterest
    subplot(maxAreasOfInterest,1,i)
    plot(time , xteo(:,i));
    title(sprintf('TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end

L = floor(fs*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
xteo = filter(g,1,xteo);

figure(fignum); fignum=fignum+1;
for i=1:maxAreasOfInterest
    subplot(maxAreasOfInterest,1,i)
    plot(time , xteo(:,i));
    title(sprintf('Lowpass TEO: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
end

%% Peak Rate Curve Fitting
figure(fignum); fignum=fignum+1;
l2err = nan(maxAreasOfInterest,1);
l1err = nan(maxAreasOfInterest,1);
for i=1:maxAreasOfInterest
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
    
    subplot(maxAreasOfInterest,1,i)
    plot(t,'b'), hold on; plot(pval,'r'); hold off;
    title(sprintf('Quad Poly Fit for band: %f , %f' , freq_bands(i,1) , freq_bands(i,2)));
    
end

l2err
l1err
[~,minband] = min(l1err);
minband = freq_bands(minband,:);
fprintf('Band Width of Interest: %f,%f\n',minband(1), minband(2));
    
    
    
%     meddt = median(dt);
    
%     absfft = abs(fft(smoothenergy));
   
%     [p , t] = findpeaks(absfft(1:floor(length(absfft)/2)));
%     [~ , maxt] = max(p);
%     maxt = t(maxt);
%     PRI = 1/((maxt-1)/length(absfft)*fs)