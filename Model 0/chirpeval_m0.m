 function [y, freq , bl , samp , times]=chirpeval_m0(xx , fs , dur , stepDur , noiseCutoff , minFreq , fadeType , baselineType , pitchEval)
% freqeval evalutes the frequency contour of a given whistle signal and
% reconstructs its melody using a sinusoidal signal.
% Input arguments:
% xx - the input  signal (should be a whistle)
% fs - the sampling frequency of the signal
% dur - the duration of each frame (seconds)
% stepDur - duration of the gap between two adjacent frames. (seconds)
% noiseCutoff - minimal signal power(in precentage, out of the total input)
% minFreq - minimal frequency for a non noise frame (in Herz). Any frame with
% frequenct below minFreq is considered noise.
% fadeType - 4 possible fade types:
%     'lpc-hamming' - lpc with hamming smoothing
%     'lpc-poly'  - lpc with polynomial smoothing
%     'exponential' - exponential fade
%     'hermite' - using hermite polynomial to fade in/out
% baseType - 3 possible baseline types:
%     'ideal low pass'
%     'moving average'
%     'amplitude average sampling'
% PitchEval - Choice of algorithm for pitch evaluation. Default is Zero
% Crossing Rate. Possible values:
%     'zcr'
%     'yin'
%     'yinbird'
% Output arguments:
% y - the output signal, a concatenation of sinusoids segments which imitate the whistle.
% freq - the evaluated frequency contour of the signal.
%
% Usage: [y, freq]=freqeval(xx,fs,dur);
% Example: [xx,fs]=audioread('KingFisher 1005 2016-02-15 12-47 1.wav'); [y, freq]=chirpeval_m0(xx,fs,0.05);


if nargin<2, fs=44100; end
if nargin<3, dur=0.005; end
if nargin<4, stepDur=0.0025; end
if nargin<5, noiseCutoff = 0.1; end
if nargin<6, minFreq = 1800; end
if nargin<7, fadeType = 'lpc-hamming'; end
if nargin<8, baselineType = 'amplitude average sampling'; end
if nargin<9, pitchEval = 'zcr'; end
[row col]=size(xx);
if col==2, xx=xx(:,1); end
Ts=1/fs;
winlen=floor(dur*fs);
amp_win_len = floor(fs*0.3*10^-3); % 3 ms amplitude winlen
if(rem(winlen,2)~=0)
    winlen = winlen +1;
    display('changing window size for smooth overlap'); 
end
tt=0:Ts:stepDur-Ts;
N=length(tt);
    
zz=zeros(size(xx));
freq = zeros(floor(length(xx)/N),1);
yin_freq = zeros(floor(length(xx)/N),1);
zcr_freq = zeros(floor(length(xx)/N),1);
indices = zeros(floor(length(xx)/N),1);
noise = zeros(floor(length(xx)/N),1);
phase = zeros(floor(length(xx)/N),1);
times = zeros(floor(length(xx)/N),1);

amp_max = zeros(length(xx),1);
amp_min = zeros(length(xx),1);
amp = zeros(length(xx),1);
y_offset = zeros(length(xx),1);

avrg_Energy=sum(xx.^2)/length(xx); % average energy of the signal.
phi=0;%pi/2;


dd = floor(amp_win_len/2);
fill = linspace(0,1,amp_win_len+1);
fill = fill(2:end-1);

prevk = 1;
for k=1+dd:amp_win_len:length(xx)-dd

    max_s = max(extrema(xx(k-dd:k+dd),1 , 1));
    min_s = min(extrema(xx(k-dd:k+dd),-1 , 1));

    if(isempty(max_s)),k , max_s = max(xx(k-dd:k+dd)); end
    if(isempty(min_s)),k , min_s = min(xx(k-dd:k+dd)); end
    amp_max(k) = max_s;
    amp_min(k) = min_s;    

    amp(k) = 0.5*(amp_max(k) - amp_min(k));
    y_offset(k) = 0.5*(amp_max(k) + amp_min(k));

    if (prevk > 1) 
    amp(prevk+1:k-1) = amp(prevk) + (amp(k) - amp(prevk)).*fill;
    y_offset(prevk+1:k-1) = y_offset(prevk) + (y_offset(k) - y_offset(prevk)).*fill;
    end
    
    

    prevk = k;
end

samp = moving_average(amp,4);
bl = moving_average(y_offset,4);

if (strcmp(baselineType , 'moving average'))
    bl = bassline(xx);
elseif (strcmp(baselineType , 'ideal low pass'))
    bl = baseline(xx,fs,0.05,150 , 'fourier');
%     bb_half = baseline(xx,fs,0.025,150 , 'fourier');
elseif (strcmp(baselineType , 'non'))
    bl = zeros(size(bl));
elseif (~strcmp(baselineType , 'amplitude average sampling'))
    error('baselineType parameter unknown value. Possible values:''moving average'' , ''ideal low pass'' , ''amplitude average sampling''');
end

xx_normalised = xx - bl;

fadeout = 0;
noise(1) = 1;
indices(1) = winlen/2;

% bl(1:winlen/2) = 0;

yinbird_flag = 0;
if (strcmp(pitchEval, 'yinbird'))
    [ f0_yinbird , dips_yinbird , time_yinbird ] = yinbird(xx_normalised , fs ,  dur , 0.1 , 16 , noiseCutoff , minFreq);
    yinbird_flag = 1;
elseif (strcmp(pitchEval, 'yin'))
    [ f0_yinbird , dips_yinbird , time_yinbird ] = yin3(xx_normalised , fs ,  dur , 0.1 , minFreq , 16);
    yinbird_flag = 1;
end

j=2;
for i=winlen/2+1:N:length(xx)-winlen+1;
    win=xx_normalised(i:i+winlen-1);
    
    indices(j) = i+N;
    times(j) = indices(j) / fs;
    
    if (yinbird_flag)
        yin_index = find(times(j)>=time_yinbird, 1 , 'last');
        if(isempty(yin_index)) yin_index = 1; end
        f0 = f0_yinbird(yin_index);
    else
        f0=zcr(win,Ts,1);
        freq(j) = f0;
    end
    
    freq(j) = f0;
    
    T0=1/f0;
    E_win=sum(win.^2)/winlen; % the average energy of each segment

    if E_win>avrg_Energy*noiseCutoff && f0 > minFreq
        %signal is not noise
        if (noise(j-1)==1 || fadeout==1)
            fadein = 1;
        else
            fadein = 0;
        end
        fadeout = 0;
        
    else
        noise(j) = 1;
        fadein = 0;
        if(noise(j-1) == 0 && fadeout == 0)
            fadeout = 1;
        else
            fadeout = 0;
        end
    end
    
    if ~noise(j) %|| fadeout
        
        if(freq(j-1)~=0 || fadeout)
            k = (freq(j) - freq(j-1))/ stepDur;
            ff = freq(j-1) + k/2 * tt;
            framesin=cos(2*pi*ff.*tt+phi);
           
        elseif(f0~=0)
            framesin=cos(2*pi*f0*tt+phi);
            
        else
            
            %signal is noise => should never actually happen unless theres
            %something wrong because was determined earlier
            framesin = zeros(size(tt));
            noise (j) = 1;
            amp(j) = 0;
        end
        
%         if (i>2000 && framesin*framesin' > 0)
%             plot(zz(i-300:i+300)), figure(1),  end % for checking phase continuity
        
        % phase correction between adjacent sinusoid segments
        if framesin(end)<framesin(end-1)
            phi=acos(framesin(end))+Ts/T0*2*pi;
            
        else
            phi=2*pi-acos(framesin(end))+Ts/T0*2*pi;
        end
        
        phase(j) = phi;
        
        %amplitude correction
        if j>1

            
            ampl = samp(i:indices(j)-1);
            framesin = framesin.*ampl';
            
        end
    
        
        zz(i:i+N-1)=framesin;
        
    else
        zz(i:i+N-1)=zeros(size(tt));
        samp(i:indices(j)-1) = 0;
%         bl(i:indices(j)-1) = 0;
        freq(j) = 0;
        zcr_freq(j) = 0;
    end

    j=j+1;

end


% figure(5)
% plot( freq , 'b')
% hold on
% plot( zcr_freq, 'g')
% hold off

%fade ins and outs:
new_amp = samp;
chirps = [];

ttin = (tt+Ts)/stepDur;   %for fade ins

%H(0)=H'(0)=H'(1)=0 , H(1)=1
hermite_env = (3*ttin.*ttin - 2*ttin.*ttin.*ttin)';

% fadein_envelope = hermite; 
exp_env = exp(6*(ttin-1))';
% fadeout_envelope = exp(-2000*tt);

%find time location of all different "utterances" or chirps in synthesized
%signal
for j=2:length(noise)
    if(noise(j)==0 && noise(j-1)==1) %start chirp
        j_start = j;
    elseif(noise(j)==1 && noise(j-1)==0) %end chirp
        chirps = [chirps ; [j_start j-1]];
    end
end

disp('Generating Chirps');
for ch = 1:size(chirps,1)
    j_start = chirps(ch,1);
    j = chirps(ch,2);
    chirp_start = indices(j_start-1);
    chirp_end = indices(j)-1;
    phase_out = phase(j);
    
    if size(chirps,1) == 1
        phase_in = 0;
        gap_left = chirp_start - 1;
        gap_right = length(xx) - chirp_end;
    elseif (ch==1)
        phase_in = 0;
        gap_left = chirp_start - 1;
        gap_right = indices(chirps(ch+1,1)) - chirp_end;
    elseif(ch==size(chirps,1))
        phase_in = phase(chirps(ch-1,2));
        gap_left = chirp_start - indices(chirps(ch-1,2));
        gap_right = length(xx) - chirp_end;
    else
        phase_in = phase(chirps(ch-1,2));
        gap_left = chirp_start - indices(chirps(ch-1,2));
        gap_right = indices(chirps(ch+1,1)) - chirp_end;
    end
    
    disp([num2str(ch) ,') ',num2str(chirp_start/fs),':',num2str(chirp_end/fs) ,' - windows:'...
                ,num2str(chirps(ch,1)),':',num2str(chirps(ch,2))]);
    switch(fadeType)
        case 'exponential'
            env_in = exp_env * samp(chirp_start);
            env_out = wrev(exp_env) * samp(chirp_end);
        case 'hermite'
            env_in = hermite_env * samp(chirp_start);
            env_out = wrev(hermite_env)* samp(chirp_end);
        otherwise
            %30 ms window for p
            [env_in , env_out] = fade_amplitude(samp(chirp_start : chirp_end) , 0.03 , fs , gap_left , gap_right);

            env_in = windowed_env(env_in , gap_left , 'in' , fadeType);
            env_out = windowed_env(env_out , gap_right , 'out' , fadeType);
    end
    
    if(gap_left < length(ttin))
        env_in = env_in((length(ttin)-gap_left+1):end);
    end
    if(gap_right< length(ttin))
        env_out = env_out(1:gap_right);
    end
        
    new_indices = chirp_end+1:chirp_end+length(env_out);
    [a , s] = ...
        create_fade_sinusoid(env_out , 'out', freq(j) , phase_out ,  fs);
    new_amp(new_indices) = new_amp(new_indices)+a;
    zz(new_indices) = zz(new_indices) + s;
    
    new_indices = chirp_start-length(env_in):chirp_start-1;
    [a , s] = ...
        create_fade_sinusoid(env_in , 'in', freq(j_start) , phase_in ,  fs);
    new_amp(new_indices) = new_amp(new_indices)+a;
    zz(new_indices) = zz(new_indices) + s;

end

% soundsc(zz,fs)
y=zz;

end

function env_out = windowed_env(env , gap , inorout , wintype)
    if (~ (strcmp(inorout,'in') || strcmp(inorout,'out'))) , error('wrong value for inorout parameter') ; end
    if (~ (strcmp(wintype,'lpc-hamming') || strcmp(wintype,'lpc-poly'))) , error('wrong value for wintype parameter') ; end
    if (strcmp(wintype , 'lpc-hamming') )
        p = 0.1;
        win = hamming(floor(gap*(0.5+p)));
        winlen = floor(length(win)/2);
    elseif(strcmp(wintype , 'lpc-poly'))
        alpha = 8;
        win = polynomial_window(gap , alpha , inorout);
        winlen = length(win);
    else
        error('unspecified window type')
    end
    
    if (length(env) > winlen) 
        if ( strcmp(inorout , 'in') ) 
            env = env(end-winlen+1:end);
        else
            env = env(1:winlen);
        end
    end
    
    
    % here we can assume length(env)<=length(win)
        
    if (strcmp(wintype , 'lpc-hamming'))
        if ( strcmp(inorout , 'in') )
            win = win(1:winlen);
        else
            win = win(end-winlen+1:end);
        end
    end
    
    if ( strcmp(inorout , 'in') )
        win = win(end-length(env)+1:end);
    else
        win = win(1:length(env));
    end
    
    env_out = win .* env;
end

function [amp_out , sig_out] = create_fade_sinusoid(env , inorout , freq , phase ,  fs)
    Ts = 1/fs;
    amp_out = env;
    durtt = length(env) / fs;

    if(strcmp(inorout , 'out'))
        tt = 0:Ts:durtt-Ts;
    elseif(strcmp(inorout , 'in'))
        tt = -Ts:-Ts:-durtt;
    end
    
    framesin=cos(2*pi*freq*tt+phase);
    if(strcmp(inorout , 'in'))
        framesin = wrev(framesin);
    end
    sig_out= framesin'.*env;
end

