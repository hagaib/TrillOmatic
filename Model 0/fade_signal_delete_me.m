function y = fade_signal( x , fs , yin)
%FADE_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

if(size(x,2)==2), x = x(:,1); end

analyticx = hilbert(x);
i_phase = angle(analyticx);
i_amp = abs(analyticx);

dx = [x(2:end) - x(1:end-1) ; 0];
indicator = medfilt1(1*(dx~=0),5);

segs_t = logical2segments(indicator , fs)';
segs_s= floor(segs_t*fs);

phase = i_phase(segs_s);
amp = i_amp(segs_s);
yin.f0(arrayfun(@(z)find(yin.time<=z , 1 , 'last') , segs_t))

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

