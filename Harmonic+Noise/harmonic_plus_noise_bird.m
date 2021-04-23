function [yy , f0 , time , harms ,  yin] = harmonic_plus_noise_bird(xx , fs ,windur , stepdur , min_freq , energy_threshold , yin)

if(size(xx,1) == 1) , xx = xx'; end


winsize = floor(windur*fs);
if mod(winsize , 2) , winsize = winsize+1; end 
hann = hanning(winsize);
stepsize = floor(stepdur*fs);


yy = zeros(size(xx));
time = nan(1,1+floor(length(xx)/stepsize));
f0 = nan(size(time));



bb = baseline(xx , fs , windur , min_freq , 'fir1');

if(exist('yin','var'))
    f0_yin = yin.f0;
    dips_yin = yin.dips;
    time_yin = yin.time;
else
%     [ f0_yin , dips_yin , time_yin ] = yin3(xx-bb , fs , windur , 0.1 , min_freq , stepsize);
    [ f0_yin , dips_yin , time_yin ] = yinbird(xx-bb , fs , windur , 0.1 ,stepsize, energy_threshold , min_freq);
%     f0(dips>0.3) = nan;
    yin = [];
    yin.f0 = f0_yin;
    yin.dips = dips_yin;
    yin.time = time_yin;
end
% low_pass = fir1(34 , 0.5);
% xx = filter(low_pass , 1 , xx);

harms = harmonics_of_signal(xx , fs , yin , 5);

av_energy = (xx-bb)'*(xx-bb)/length(xx);
j=1;
for i=winsize/2:stepsize:(length(xx)-winsize/2)
    if(j>length(time))
    end
    
    samples = i-winsize/2+1:i+winsize/2;
    win = xx(samples) - bb(samples);
    win_energy = win'*win/winsize;
    timestamp = i/fs;
    yindex = find(time_yin<timestamp, 1 , 'last');
    
    f0_yin =  yin.f0(yindex);
    
    time(j) = timestamp;
    
%     if (~isempty(yindex) && win_energy > energy_threshold * av_energy && f0_yin>min_freq && dips_yin(yindex) < 0.5)    
    if (~isempty(yindex) && f0_yin>min_freq)    
        [y , f0fit , A] = harmonic_plus_noise_win(xx(samples)-bb(samples) , samples, fs, hann , f0_yin ,  harms(yindex , :));
        f0(j) = f0fit;
        yy(samples) = hann.*y + yy(samples);
    end  
    j=j+1;

end

function [yy , f0fit , A ,phase , frequencies] = harmonic_plus_noise_win(xx , tt, fs , winfunc , f0 , frequencies)

frequencies = frequencies(:);
frequencies = [f0 ; frequencies];

 
L = length(frequencies);
harmonics = (1:L)';
harmonics = harmonics(~isnan(frequencies));
frequencies = frequencies(~isnan(frequencies));
L = length(frequencies);
 
f0fit = harmonics \ frequencies; %simple least squares regression on data


B = repmat(tt'/fs,1, 2*L+1);
B = B * f0fit;
B = B * diag([-wrev(harmonics); 0 ; harmonics]);
B = exp(2*pi*1i*B);

R = (B)'*diag(winfunc.^2)*B;
b = (B)'* diag(winfunc.^2)*xx;
A = R \ b;
% figure(2)
% plot(real(B*A)) , hold on
% plot(xx)  , hold off
% legend('Least Sq','Original');
% sum(abs(A(A~=A((length(A)+1)/2))))
yy = real(B*A);


