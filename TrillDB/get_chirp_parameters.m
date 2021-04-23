function [t , A , f0] = get_chirp_parameters(xx , fs , windur , stepdur , threshold_energy , threshold_zcr)
% Extracts information for chirp generation


xx = xx(:,1);

bpfilt = designfilt('bandpassfir','FilterOrder',48, ...
         'CutoffFrequency1',threshold_zcr,'CutoffFrequency2',fs*0.45, ...
         'SampleRate',fs);
xx = filtfilt(bpfilt , xx);


AA = hilbert(xx);
lpfilt = fir1(48,0.05);
AA = filtfilt(lpfilt ,1 , abs(AA));




winsize = floor(windur * fs);
if(mod(winsize,2)) , winsize = winsize+1; end

t = 0:stepdur:(length(xx)/fs);
A = zeros(size(t));
f0 = zeros(size(t));
w = hamming(winsize+1);

av_E = sum(xx.^2)/length(xx);

N = 2048;

istart = find(t>=(winsize/2+1)/fs ,1, 'first');
iend = find(t<=t(end)-windur/2,1,'last');
for i= istart:iend

    n = 1 + floor(t(i)*fs);
    win = xx(n-winsize/2:n+winsize/2);
    win_E = sum(win.^2)/length(win);
    zcrwin = zcr(win , 1/fs , 1);
    if(win_E > av_E * threshold_energy && zcrwin > threshold_zcr)
        
        P = abs(fft(win.*w , N));
        P(2:(N/2)-1) = 2*P(2:(N/2)-1);
        P = P(1:N/2);
        
        [ peaks , ind_peaks , lobe_ampl ] = get_lobes(P);
        
        [spectE , bin] = max(P(3:winsize/2));
        bin = bin+1;
        f0(i) = bin/N*fs;
        A(i) = AA(n);
        
%         f0_lobe = lobe_ampl(ind_peaks == bin+1);
%         
%         
%         
%         j= 2*(bin+1);
%         range = [floor(-bin/2) , floor(bin/2)];
%         while (j<length(P))
%             Prange = range + j;
% %             if(range(end) > length(P)) ,  break ,  end
%             peak_range = ind_peaks>Prange(1) & ind_peaks<Prange(2);
%             [m ,ind] = max(peaks(peak_range)); 
%             bins = ind_peaks(peak_range);
%            
%            
%             
%             j = j + bin+1;
%         end
        
    end
end

not_voiced = ~(A | f0);
new_f0 = f0;
for i=1:length(not_voiced)
    if( not_voiced(i) )
        if(i==637)
        end
        if(i>1 &&f0(i-1)) , new_f0(i) = f0(i-1); 
        elseif(i<length(f0) && f0(i+1)) , new_f0(i) = f0(i+1); end
    end
end
f0 = new_f0;