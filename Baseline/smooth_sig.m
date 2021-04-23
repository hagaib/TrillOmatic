function [sigout, counter] = smooth_sig(sigin , win_width , average_width , min_zcr)
%smoothes out a signal by averaging over windows with high volatility
%   win_width     = sample window width 
%   average_width = calculate moving average with 2*average_width+1 samples
%   min_zcr       = averaging is preformed on every sample on window only if zcr of window larger than min_zcr
    N = length(sigin);
    sigout = sigin;
    counter = 0;
    for i=1:win_width:N-win_width+1
        win = sigin(i : i+ win_width -1);
        diff = win(2:win_width) - win(1:win_width-1);
        
        if(zcr(diff,0,0) > min_zcr)
            
            x = moving_average(sigin,average_width,i,i+win_width-1);
            sigout(i:i+win_width-1) = x(i:i+win_width-1);
            counter = counter+1;
        end
        
    end
    
%     figure(2)
%     plot(sigin,'b')
%     hold on
%     plot(sigout,'g')
%     hold off
end

