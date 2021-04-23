function [ fadein , fadeout ] = fade_amplitude(amp , dur_p,fs , maxin , maxout)
% fade_amplitude calculates fade in and out aplitudes using lpc
% coefficients in both directions (forwards and backwards)

try
    samp = moving_average(amp,floor(fs * 0.0025));
catch ME
    if(strcmp(ME.identifier , 'moving_average:not_enough_samples'))
        disp('not enough samples for smoothing');
    end
    samp = amp;
end
rev_amp = wrev(amp);
L = length(samp);
[M,I] = max (samp);


N = min(floor(dur_p*fs) , L-I);
if (N<floor(0.2*dur_p*fs)) , N = floor(L/2); end % minimum samples for lpc calculation
a = lpc(amp(L-N+1:end),N);
i=0;
x=1;

while(x>0 && length(amp(L+1:end))<maxout) 
    x = -a(2:end)*wrev(amp(L-N+1+i:L+i));
    amp = [amp ; x];
    i=i+1;
end



N = min(floor(dur_p*fs) , I);
if (N<floor(0.2*dur_p*fs)) , N = floor(L/2); end
a_rev = lpc(rev_amp,N);
i=0;
x=1;

while(x>0 && length(rev_amp(L+1:end)) < maxin)
    x = -a_rev(2:end)*wrev(rev_amp(L-N+1+i:L+i));
    rev_amp = [rev_amp; x];
    i=i+1;
end



fadein = wrev(rev_amp(L+1:end));
fadeout = amp(L+1:end);
% figure(3);
% plot([fadein ;  amp(1:L) ; fadeout]);
end
