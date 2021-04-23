function energy = hamming_energy(xx , fs ,  windur , step)
%calculates energy of hamming window
%example:
%hamming_energy(weighted_bands, fs ,  0.01 , 1);
energy = zeros(size(xx));
winsize = floor(windur*fs);
win = hamming(winsize);
xx = xx(:);
L = length(xx);
transient = ceil(winsize/2);
xx = [zeros(transient , 1) ; xx ; zeros(transient , 1)];
for i=1:L
    xxwin = xx(i:i+winsize-1);
    e = win .* xxwin;
    e = sum(e.^2);
    energy(i) = e;
end
energy = energy/winsize;