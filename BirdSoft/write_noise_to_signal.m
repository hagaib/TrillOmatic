function write_noise_to_signal( filename , snr )
% yy is assumed to be a clean signal
% snr is a vector of desired SNR values in db
% there is no output.
% function will prepare the same number of files as length(snr)
% each filename will be same as original , with a __db postfix

snr = snr(:)';

[yy , fs] = audioread(filename);

noise = randn(size(yy)); %variance = 1 = mean energy when expectance=0

yye = mean(yy.^2);

noisecoeffs = 10.^(-snr./20)*sqrt(yye); 
% so that 10*log10(mean(yy.^2)/mean(noise.^2)) = snr

noisemat = repmat(noise , 1 , length(snr)) * diag(noisecoeffs);

postfix = [];
if(isequal(filename(end-3) , '.'))
    postfix = filename(end-3:end);
    filename = filename(1:end-4);
end
for i=1:length(snr)
    out = yy + noisemat(:,i);
    newname = sprintf([filename , '_' , num2str(snr(i)) , 'db' , postfix]);
    audiowrite(newname , out , fs);
end
