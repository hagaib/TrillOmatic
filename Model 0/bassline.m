function bl = bassline( y )
    %input: signal y
    %output: "smooth" mean of signal = bassline
%     y_average = moving_average(y,50);
    [bl, counter] = smooth_sig(y,100,60,5);
    while (counter > 0)
        [bl, counter] = smooth_sig(bl,100,60,5);
    end

    % figure(2)
    % plot(y_average,'b')
    % hold on
    % plot(y_bassline,'g')

    for i=4:-1:3
        counter=1;
        while (counter > 0)
            [bl, counter] = smooth_sig(bl,50,20,i);
            bl = moving_average(bl,4);
        end
    end
    
%     bl = moving_average(bl,10);
end

