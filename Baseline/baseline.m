function bb = baseline(xx , fs , dur , freq_cutoff , filter_type)
%baseline returns baseline of signal xx, based on spectral analysis.
% dur is length of window on which analysis is preformed.
% freqcutoff is the maximal frequency in the baseline. (all frequeny
% components above freqcutoff are zeroed out).

    
    winlen = floor(dur * fs);
    if (mod(winlen,2)) , winlen = winlen+1; end
    halflen = winlen/2;
   
    hamm = hamming(winlen);
    
    switch(filter_type)
        case 'fourier'
            bb = bl_fourier(xx , fs , hamm , winlen , halflen , freq_cutoff);
            
        case 'fir1'
            bb = bl_fir1(xx , fs , hamm , winlen , halflen , freq_cutoff);
            
        case 'median'
            bb = bl_median(xx , hamm , winlen , halflen);
%             bb = moving_average(bb , 3);
    end

end

function bb = bl_fir1(xx , fs , hamm , winlen , halflen , freq_cutoff)
    
    
    w = freq_cutoff/fs * 2;
    blo = fir1(248 , w);
    
    bb = filter(blo , 1 , xx);
    
    transient = mean(grpdelay(blo));
    bb(1:end-transient) = bb(transient+1:end);
    bb(end-transient+1:end) = 0;
    
end

function bb = bl_median(xx , hamm , winlen , halflen)
    bb = zeros(length(xx),1);

    for i=1:halflen:length(xx)-winlen
        x = xx(i:i+winlen-1) .* hamm;
        b = medfilt1(x , 300);
        bb(i:i+winlen-1) = bb(i:i+winlen-1) + b.*hamm;
    end
end


function bb = bl_fourier(xx , fs , hamm , winlen , halflen , freq_cutoff)
    bb = zeros(length(xx),1);
    k=0;
    while ((k+1)/winlen*fs <= freq_cutoff) , k=k+1;  end

    for i=1:halflen:length(xx)-winlen
        x = xx(i:i+winlen-1);
        f = fft(x .* hamm);
        f(1+k+1:end-k)=0;
        b = ifft(f);
        bb(i:i+winlen-1) = bb(i:i+winlen-1) + b.*hamm;
    end
end
