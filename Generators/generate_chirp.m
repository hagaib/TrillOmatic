function yy = generate_chirp(fs , f0s , times , amps , harmonics_amp_ratio , init_phase_setting)
%Chirp generation based on parameters
% chirp is linear!
% f0s = vector of fundamental frequencies of analysis points
% times = time instances of analysis points
% amps = amplitudes in analysis points
% harmonics = number of harmonics generated
% harmonics_amp_ratio = array: log-amlitude ratio: fundamental to harmonic (in db) 
% init_phase_setting = 0 for aligned , 1 for random (default: aligned) 

if(nargin < 6) , init_phase_setting = 0; end
if (nargin < 5), harmonics_amp_ratio = 0; 
else harmonics_amp_ratio = horzcat(0 , harmonics_amp_ratio);
end
harmonics_amp_ratio = 10.^(-harmonics_amp_ratio/20);
harmonics = length(harmonics_amp_ratio)-1;
if (init_phase_setting) , phase = rand(1 ,harmonics+1)*2*pi;
else
    phase = zeros(1, harmonics+1);
end

samplecount = 1;

yy = zeros(1,floor(times(end)*fs));
Ts = 1/fs;

for i = 1:(length(times)-1)
dur = times(i+1) - times(i);
tt = 0:Ts:(dur-Ts);
k = (f0s(i+1) - f0s(i))/ dur;
ff = f0s(i) + k/2 * tt;
aa = (amps(i+1)*tt + amps(i)*(dur-tt))/dur;

for  j=1:harmonics+1
    
    y = cos(j*2*pi*ff.*tt + phase(j));

    phase(j) = acos(y(end));
    if(y(end-1) < y(end)) , phase(j) = 2*pi - phase(j); end
    phase(j) = phase(j) + 2*pi/(fs/f0s(i+1)/j);

    yy(samplecount:samplecount+length(tt)-1) = ...
        yy(samplecount:samplecount+length(tt)-1) + harmonics_amp_ratio(j) .* aa .* y;
end
samplecount = samplecount + length(tt);

end

