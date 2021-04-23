function y = yinslopevarfilt(yin , windur)
%%filters signal based on yin structure.
%filtration is performed by calculating slopes of frequency vector 
%(output of yin f0 estimation) and calculating sample variance in window
%specified by windur.
%output is centered around window (not causal)
%transients are replaced by 'Inf'
%%ex: yinslopevarfilt(yin , 0.02);


winlen = floor(windur*yin.steprate);
if(~mod(winlen , 2)) , winlen=winlen+1; end
% hwin = hamming(winlen)';

dyin = yin.f0(2:end) - yin.f0(1:end-1);
% dyin(isnan(dyin)) = max(abs(dyin));
dyin = [0 , dyin];


y = inf(length(dyin) , 1);
for i=1:length(dyin)-winlen+1
    win = dyin(i:i+winlen-1);
    y(i) = var(win , 'omitnan');
end
% dyinloc = filter(win , 1 , dyin);
delay = (winlen-1)/2;
y = [inf(delay , 1) ; y(1:end-delay)];


end

