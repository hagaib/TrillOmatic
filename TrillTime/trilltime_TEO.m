function [trilltime , threshold , teo , orig_teo] = trilltime_TEO(xx , fs)
%calculate trilltime according to TEO thresholding
%threshold is found according to beginning of signal

transient_time = 0.01;
transient_samples = floor(transient_time*length(xx));
transient_samples = round(transient_samples, -1);
transient = 1:round(transient_samples, -1);
transient = reshape(transient , [floor(transient_samples/5) , 5]);

L = floor(fs*0.01);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);


orig_teo = teager(xx);
orig_teo = [orig_teo(1); orig_teo ; orig_teo(end)];

teo = filter(g , 1 , orig_teo);
teo = teo(:);
teo = [teo((L+1)/2:end) ; zeros((L-1)/2 , 1)];


%threshold is calculated according to original teo (non-smooth)
threshold = mean(prctile(orig_teo(transient) , 90));

istrill = find(teo > threshold);

trilltime = [istrill(1) , istrill(end)]/fs;
