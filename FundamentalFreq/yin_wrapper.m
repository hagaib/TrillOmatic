function yin = yin_wrapper(xx , fs , stepdur , min_freq)
% wrapper function for YIN fundamental frequency estimator.
% inputs:
% xx: input 1D signal
% fs: sample rate
% stepdur: duration of interval between analysis points in seconds.
% min_freq: minimal fundamental frequency in Herz.
%
%outputs:
% yin: structure with the following fields:
% yin.f0 = estimated fundamental frequency at analysis points
% yin.dips = dip height at respective analysis points
% yin.time = times of analysis points
% yin.steprate = 1/stepdur

windur = 0.0015; % around 2.7 pitch periods (1/2000)

if(nargin<4)
    min_freq = 1800; % kinfisher case
end

if(nargin<3)
    stepdur = 0.001;
end

stepsize = floor(stepdur*fs);
if(min_freq>-inf)
    b = fir1(256 , min_freq/fs*2 , 'high');
    delay = mean(grpdelay(b,1));
    xx_foryin = filter(b, 1 , xx);
    xx_foryin = [xx_foryin(delay+1:end) ; zeros(delay , 1)];
else
    xx_foryin = xx;
end

[ f0 , dips , time ] = yin3(xx_foryin, fs , windur , 0.1 , min_freq, stepsize);
% [ f0 , dips , time ] = yinbird(xx, fs , windur , 0.1 ,stepsize, 0.1 , min_freq);
yin.f0 = f0;
yin.dips = dips;
yin.time = time;
yin.steprate = 1/stepdur;