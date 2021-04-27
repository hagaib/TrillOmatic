function xxfilt = filter_signal(xx , fs , flow , fhigh)
%filtering the signal to desired spectrum band
if flow == 0
    if fhigh == fs/2; return; end
    % lowpass
    n = cheb2ord(fhigh/fs*2 , fhigh*1.02/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, fhigh/fs*2 , 'low');
elseif fhigh >= fs/2
    % highpass
    n = cheb2ord(flow/fs*2 , flow*0.98/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, flow/fs*2 , 'high');
else
    n = cheb2ord([flow, fhigh]/fs*2 , [flow*0.98 , fhigh*1.02]/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, [flow , fhigh]/fs*2 , 'bandpass');
end
sos = zp2sos(z,p,k);
xxfilt = sosfilt(sos , xx);
