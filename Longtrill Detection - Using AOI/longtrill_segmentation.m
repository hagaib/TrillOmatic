function [detect , segs , plotparams] = longtrill_segmentation(x , fs , energy , t_energy_peaks , segmentation_type , short_energy)
%LONGTRILL_SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here


switch (segmentation_type)
    case 'median'
        [detect , low] = median_segmenation(x , fs , energy , t_energy_peaks);
    case 'margin'
        [detect , low] = marginal_energy_segmenation(x , fs , energy , t_energy_peaks);
    case 'adaptive'
        if(nargin<6) , error('When adaptive option is chosed, a shorter window energy should be supplied as well');end
        [detect , low] = syllable_adaptive_energy_segmenation(x , fs , energy , t_energy_peaks , short_energy);
        
end

outlier_dur_t = 0.001; %ms
outlier_dur_s = ceil(outlier_dur_t*fs);
outlier_dur_s = outlier_dur_s + ( 1 - mod(outlier_dur_s , 2) );

detect = medfilt1(detect , outlier_dur_s);

segs = logical2segments(detect , fs);
plotparams.low = low;




function [detect , low] = median_segmenation(x , fs , energy , t_energy_peaks)

%% Detecting pulse interval
detect = zeros(size(x));
% windur = 0.05;
% low = prcfilt(teotime, floor(windur*fs) , 5 , 'omitnan')+threshold_teo;
low = 1.5*median(energy);

for i=1:length(t_energy_peaks)
    istart = floor(t_energy_peaks(i)*fs);
    iend = floor(t_energy_peaks(i)*fs);
    while(istart>0 && energy(istart)>low) , istart = istart-1; end
    while(iend<length(x) && energy(iend)>low) , iend = iend+1; end
        
    detect(istart:iend) = 1;
end
low = low*ones(size(energy));

function [detect , low] = marginal_energy_segmenation(x , fs , energy , t_energy_peaks)

detect = zeros(size(x));

transient_time = 0.1;
transient_samples = floor(transient_time*length(x));
transient_samples = round(transient_samples, -1);
transient = 1:round(transient_samples, -1);
transient = reshape(transient , [floor(transient_samples/5) , 5]);

threshold = mean(prctile(energy(transient) , 90));
low = 1.5*threshold;

for i=1:length(t_energy_peaks)
    istart = floor(t_energy_peaks(i)*fs);
    iend = floor(t_energy_peaks(i)*fs);
    while(istart>0 && energy(istart)>low) , istart = istart-1; end
    while(iend<length(x) && energy(iend)>low) , iend = iend+1; end
        
    detect(istart:iend) = 1;
end
low = low*ones(size(energy));



function [detect , low] = syllable_adaptive_energy_segmenation(x , fs , energy , t_energy_peaks , energy_short)

p = 0.7;%0.6;
psoft = 0.3;

detect = zeros(size(x));
low = zeros(size(detect)); low2 = zeros(size(detect));
s_energy_peaks = floor(t_energy_peaks*fs);
s_energy_dips = zeros(length(s_energy_peaks)+1,1);
energy_dips = zeros(size(s_energy_dips));
energy_dips(1) = median(energy_short(1:s_energy_peaks(1))); energy_dips(1) = 20*log10(energy_dips(1));
energy_dips(end) = median(energy_short(s_energy_peaks(end):end));
s_energy_dips(1) = 1; s_energy_dips(end) = length(x);
energy_peaks = energy_short(s_energy_peaks);

for k=1:length(s_energy_peaks)
%     [highenergy(k) , hightime(k)] = max(energy(lowtime(k):lowtime(k+1)));
%     hightime(k) = hightime(k) + lowtime(k)-1;
%         disp(['high: ' , num2str(hightime(k)/fs)])
    if(k<length(s_energy_peaks))
        [energy_dips(k+1) , s_energy_dips(k+1)] = min(energy_short(s_energy_peaks(k):s_energy_peaks(k+1)));
        s_energy_dips(k+1) = s_energy_dips(k+1) + s_energy_peaks(k)-1;
    end
    energy_peaks(k) = 20*log10(energy_peaks(k));
    energy_dips(k+1) = 20*log10(energy_dips(k+1));

    istart = s_energy_peaks(k);
    iend = istart;

    rightthresh = 10^((p*energy_peaks(k) + (1-p)*energy_dips(k+1))/20);
    leftthresh = 10^((p*energy_peaks(k) + (1-p)*energy_dips(k))/20);
    softrightthresh = 10^((psoft*energy_peaks(k) + (1-psoft)*energy_dips(k+1))/20);
    softleftthresh = 10^((psoft*energy_peaks(k) + (1-psoft)*energy_dips(k))/20);

%     rightthresh = p*energy_peaks(k) + (1-p)*energy_dips(k+1);
%     leftthresh = p*energy_peaks(k) + (1-p)*energy_dips(k);
%     softrightthresh = psoft*energy_peaks(k) + (1-psoft)*energy_dips(k+1);
%     softleftthresh = psoft*energy_peaks(k) + (1-psoft)*energy_dips(k);

    goflag = true;
    while(goflag)
        goflag = false;
        if(istart> s_energy_dips(k))
            if(energy(istart) > leftthresh|| energy_short(istart)>leftthresh)
                goflag = true;

            elseif(energy_short(istart) < leftthresh)
                if(energy_short(istart) > softleftthresh)
                    goflag = true;
%                     elseif(de(istart)<0 || de(istart) > slopethresh) ...
%                             || (sde(istart)<0 || sde(istart) > slopethresh*20) 
%                         goflag = true;
                end
            end
        end            
        if(goflag)
            istart = istart-1;
        end
    end

    goflag = true;        
    while(goflag)
        goflag = false;
        if(iend < s_energy_dips(k+1))

            if(energy(iend) > rightthresh || energy_short(iend)>rightthresh)
                goflag = true;
            elseif(energy_short(iend) < rightthresh)

                if(energy_short(iend) > softrightthresh)
                    goflag = true;
%                     elseif( de(iend)>0 || de(iend) < -slopethresh )
%                         %%||(sde(iend)>0 || sde(iend) < -slopethresh*20)
%                         goflag = true;
                end
            end
        end

        if(goflag)
            iend = iend+1;
        end
    end

    detect(istart:iend) = 1;
    low(istart:s_energy_peaks(k)) = softleftthresh;
    low(s_energy_peaks(k):iend) = softrightthresh;
    low2(istart:s_energy_peaks(k)) = leftthresh;
    low2(s_energy_peaks(k):iend) = rightthresh;
end

t = (0:(length(energy)-1))/fs;
% figure(100000)
% plot(t,x/10/2/2); hold on;
% plot(t,energy);
% plot(t,energy_short);
% plot(t,low2);
% plot(t,low);
% plot(t(s_energy_peaks), 10.^(energy_peaks/20) , '*');
% plot(t(s_energy_dips), 10.^(energy_dips/20) , '*');
% hold off