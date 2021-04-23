function [env , detectlogic , detectsegs , PRI] = energy_envelope2(energy, fs , time , axes , bulktime  , xx)

    detectlogic = zeros(size(energy));
    detectsegs = [];
    
    
%     smallsyl = filter(b,1,wenv);
    L =floor(0.02*fs);
    b = ones(1,L)/L;
    if(mod(L,2)==0) , L=L+1; end
    win = gausswin(L);
    win = win/sum(win);
    smoothenergy = filter(win , 1 , energy);
    smoothenergy= [smoothenergy((L+1)/2:end) ; zeros((L-1)/2 , 1)];
    
    
    %find peaks to calculate pulse rate
    [p , t] = findpeaks(smoothenergy(floor(bulktime(1)*fs):floor(bulktime(2)*fs)) , 'SortStr','descend');
    t = bulktime(1) + t/fs;
    t = sort(t(1:min(10 , floor(length(t)*0.7))));
    dt = t(2:end)-t(1:end-1);
    
    
    
    meddt = median(dt)
    
    absfft = abs(fft(smoothenergy));
%     psdefilt = medfilt1(psde , 3);
    
%     figure (5)
%     plot((0:length(smoothenergy)-1)/length(smoothenergy)*fs , psde)
%     hold on
%     plot((0:length(smoothenergy)-1)/length(smoothenergy)*fs , psdefilt)
%     hold off
    [p , t] = findpeaks(absfft(1:floor(length(absfft)/2)));
    [~ , maxt] = max(p);
    maxt = t(maxt);
    PRI = 1/((maxt-1)/length(absfft)*fs)
    
    env = envelope(energy , floor(fs*PRI*0.8),'peak' );
%     smallsyl = envelope(wenv , floor(fs*0.001) , 'peak' );
    longsyl = envelope(energy , floor(fs*PRI*0.25) , 'peak' );
    henv = envelope(xx);
    
      
    [high , hightime] = extrema(longsyl , 1 , 0);
    [low , lowtime] = extrema(longsyl , -1 , 0);
    
    
    M=hightime;
    m=lowtime;
    
    lowtime = lowtime/fs;
    hightime = hightime/fs;
    
    hightime = hightime(hightime>bulktime(1) & hightime < bulktime(2));
    if isempty(find(lowtime>hightime(end) , 1))
        hightime = hightime(1:end-1);
    end
    if isempty(find(lowtime<hightime(1) , 1 , 'last'))
        hightime = hightime(2:end);
    end
    lowind = [find(lowtime<hightime(1) , 1 , 'last') , find(lowtime>hightime(end) , 1)];
    lowtime = lowtime(lowind(1):lowind(2));
    
    
    hightime = floor(hightime*fs);
    lowtime = floor(lowtime*fs);
    highenergy = energy(hightime);
    lowenergy = energy(lowtime);
    
    de = energy(2:end) - energy(1:end-1);
    sde = smoothenergy(2:end) - smoothenergy(1:end-1);
    slopethresh = 10^-7;
    
    p = 0.8;
    psoft = 0.6;
    
    lowenergy(1) = 20*log10(lowenergy(1));
    
    
    rightthresharr = zeros(size(hightime));
    leftthresharr = zeros(size(hightime));
    softrightthresharr = zeros(size(hightime));
    softleftthresharr = zeros(size(hightime));
    
    for k=1:length(hightime)
        [highenergy(k) , hightime(k)] = max(energy(lowtime(k):lowtime(k+1)));
        hightime(k) = hightime(k) + lowtime(k)-1;
%         disp(['high: ' , num2str(hightime(k)/fs)])
        if(k<length(hightime))
            [lowenergy(k+1) , lowtime(k+1)] = min(energy(hightime(k):hightime(k+1)));
            lowtime(k+1) = lowtime(k+1) + hightime(k)-1;
%             disp(['low: ' , num2str(lowtime(k+1)/fs)])
        end
        highenergy(k) = 20*log10(highenergy(k));
        lowenergy(k+1) = 20*log10(lowenergy(k+1));
        
        istart = hightime(k);
        iend = istart;
        
        rightthresh = 10^((p*highenergy(k) + (1-p)*lowenergy(k+1))/20);
        leftthresh = 10^((p*highenergy(k) + (1-p)*lowenergy(k))/20);
        softrightthresh = 10^((psoft*highenergy(k) + (1-psoft)*lowenergy(k+1))/20);
        softleftthresh = 10^((psoft*highenergy(k) + (1-psoft)*lowenergy(k))/20);
     
        rightthresharr(k) = rightthresh;
        leftthresharr(k) =leftthresh;
        softrightthresharr(k) =softrightthresh;
        softleftthresharr(k) =softleftthresh ;
        
        goflag = true;
        while(goflag)
            goflag = false;
            if(istart> lowtime(k))
                if(smoothenergy(istart) > leftthresh|| energy(istart)>leftthresh)
                    goflag = true;

                elseif(energy(istart) < leftthresh)
                    if(energy(istart) > softleftthresh)
                        goflag = true;
%                     elseif(de(istart)<0 || de(istart) > slopethresh) ...
%                             || (sde(istart)<0 || sde(istart) > slopethresh*20) 
%                         goflag = true;
                    end
                end
            end            
            if(goflag)
                istart = istart-1;
            end
        end
        
        goflag = true;        
        while(goflag)
            goflag = false;
            if(iend < lowtime(k+1))
    
                if(smoothenergy(iend) > rightthresh || energy(iend)>rightthresh)
                    goflag = true;
                elseif(energy(iend) < rightthresh)
                    
                    if(energy(iend) > softrightthresh)
                        goflag = true;
%                     elseif( de(iend)>0 || de(iend) < -slopethresh )
%                         %%||(sde(iend)>0 || sde(iend) < -slopethresh*20)
%                         goflag = true;
                    end
                end
            end
            
            if(goflag)
                iend = iend+1;
            end
        end
        
        detectlogic(istart:iend) = 1;
   
        detectsegs(:,k) = [istart/fs , iend/fs];
    end
    
    if(~isempty(axes))
        hold(axes , 'on');
        plot(axes, time , detectlogic , 'g')
%         plot(axes, time , energy*100 , time , env*100 , time , smoothenergy*100 , time , longsyl*100 , time , detectlogic*0.4 , time , henv);
%         legend(axes , 'xx' , 'energy' , 'env 0.05' , 'smooth e' , 'longsyl' , 'detect' , 'hilbert')
        legend(axes , 'Signal' , 'Indicator')
        hold(axes , 'off');
    end
    
%% Plot for the presentation
%%
%     
%     figure(100)
%     E_mult = 150;
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , energy*E_mult ,'k', time , smoothenergy*E_mult, 'g', time , longsyl*E_mult , 'r' , 'Linewidth' , 2)
%     hold off
%     legend('Band Pass Signal' , 'Energy' , 'Gaussian Filtered Energy' , 'Spline Smoothed Energy')
%     
%     figure(101)
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , longsyl*E_mult , 'r' ,'Linewidth' , 2);
%     plot(time , energy*E_mult , 'k' ,'Linewidth' , 2);
%     plot(time(m) , longsyl(m)*E_mult, 'dg')
%     plot(time(M) , longsyl(M)*E_mult, 'dg')
%     hold off
%     legend('Bandpass Signal' , 'Spline Smoothed Energy' , 'Energy', 'Local Extrema')
% 
%     
%     figure(102)
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , energy*E_mult ,'k' , time , longsyl*E_mult , 'r' ,'Linewidth' , 2);
%     plot(time(m) , longsyl(m)*E_mult, 'dg' ,'Linewidth' , 2)
%     plot(time(lowtime) , energy(lowtime)*E_mult, 'dm' ,'Linewidth' , 2)
%     plot(time(M) , longsyl(M)*E_mult, 'dg')
%     plot(time(hightime) , energy(hightime)*E_mult, 'dm')
%     %%lines
%     plot([time(lowtime) , time(lowtime)] , get(gca,'ylim') , 'm')
%     hold off
%     legend('Bandpass Signal' ,'Energy', 'Spline Smoothed Energy' , 'Spline Extrema' , 'Energy Extrema')
%     
%     
%     figure(103)
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , energy*E_mult ,'k' , time , smoothenergy*E_mult, 'g' ,'Linewidth' , 2);
%     %%lines
%     plot(get(gca , 'xlim')  ,[rightthresharr(7) , rightthresharr(7)]*E_mult , 'c')
%     plot(get(gca , 'xlim')  ,[softrightthresharr(7) , softrightthresharr(7)]*E_mult , 'r')
%     plot(time(lowtime) , energy(lowtime)*E_mult, 'dm' ,'Linewidth' , 2)
%     plot(time(hightime) , energy(hightime)*E_mult, 'dm')
%     plot(time(lowtime(2:end)) ,rightthresharr*E_mult , 'oc')
%     plot(time(lowtime(2:end)) ,softrightthresharr*E_mult , 'or')
%     hold off
%     legend('Bandpass Signal' ,'Energy', 'Gaussian Smooth Energy' , 'Right Threshold' , 'Soft Right Threshold')
%     
%     
%     
%     figure(104)
%     p1 = subplot(2,1,1);
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , energy*E_mult ,'k' , time , smoothenergy*E_mult, 'g' ,'Linewidth' , 2);
%     %%lines
%     plot(get(gca , 'xlim')  ,[rightthresharr(7) , rightthresharr(7)]*E_mult , 'c')
%     plot(get(gca , 'xlim')  ,[softrightthresharr(7) , softrightthresharr(7)]*E_mult , 'r')
%     plot(time(lowtime) , energy(lowtime)*E_mult, 'dm' ,'Linewidth' , 2)
%     plot(time(hightime) , energy(hightime)*E_mult, 'dm')
%     plot(time(lowtime(2:end)) ,rightthresharr*E_mult , 'oc')
%     plot(time(lowtime(2:end)) ,softrightthresharr*E_mult , 'or')
%     plot(time ,[zeros(size(detectlogic(1:hightime(7)-1))); detectlogic(hightime(7):ceil(0.7624*fs)); zeros(length(detectlogic)-ceil(0.7624*fs),1)]*energy(hightime(7))*E_mult , 'm')
%     hold off
%     legend('Bandpass Signal' ,'Energy', 'Gaussian Smooth Energy' , 'Right Threshold' , 'Soft Right Threshold')
%     title('Slide Right Step 1')
%     
%     p2 = subplot(2,1,2);
%     plot(time ,  xx ,'b')
%     hold on
%     plot(time , energy*E_mult ,'k' , time , smoothenergy*E_mult, 'g' ,'Linewidth' , 2);
%     %%lines
%     plot(get(gca , 'xlim')  ,[rightthresharr(7) , rightthresharr(7)]*E_mult , 'c')
%     plot(get(gca , 'xlim')  ,[softrightthresharr(7) , softrightthresharr(7)]*E_mult , 'r')
%     plot(time(lowtime) , energy(lowtime)*E_mult, 'dm' ,'Linewidth' , 2)
%     plot(time(hightime) , energy(hightime)*E_mult, 'dm')
%     plot(time(lowtime(2:end)) ,rightthresharr*E_mult , 'oc')
%     plot(time(lowtime(2:end)) ,softrightthresharr*E_mult , 'or')
%     plot(time ,[zeros(size(detectlogic(1:hightime(7)-1))); detectlogic(hightime(7):end)]*energy(hightime(7))*E_mult , 'm')
%     hold off
%     title('Slide Right Step 2')
%     linkaxes([p1,p2],'x');


%% Article Images
    E_mult = 150;
    figure(101)
    plot(time ,  xx ,'b')
    hold on
    plot(time , longsyl*E_mult , 'r' ,'Linewidth' , 2);
    plot(time , energy*E_mult , 'k' ,'Linewidth' , 2);
    plot(time(m) , longsyl(m)*E_mult, 'dr' , 'MarkerSize' , 16 , 'Linewidth' , 2)
    plot(time(M) , longsyl(M)*E_mult, 'dr' , 'MarkerSize' , 16 , 'Linewidth' , 2)
    hold off
    legend('Bandpass Signal' , 'Spline Smoothed Energy' , 'Energy', 'Local Extrema')
    
    figure(102)
    plot(time ,  xx ,'b')
    hold on
    plot(time , energy*E_mult ,'k' , time , longsyl*E_mult , 'r' ,'Linewidth' , 2);
    plot(time(m) , longsyl(m)*E_mult, 'dr'  , 'MarkerSize' , 16 , 'Linewidth' , 2)
    plot(time(lowtime) , energy(lowtime)*E_mult, 'dk'  , 'MarkerSize' , 16 , 'Linewidth' , 2)
    plot(time(M) , longsyl(M)*E_mult, 'dr' , 'MarkerSize' , 16 , 'Linewidth' , 2)
    plot(time(hightime) , energy(hightime)*E_mult, 'dk' , 'MarkerSize' , 16 , 'Linewidth' , 2)
    %%lines
    plot(time , detectlogic*0.08 , 'm');
%     plot([time(lowtime) , time(lowtime)] , get(gca,'ylim') , 'm')
    hold off
    legend('Bandpass Signal' ,'Energy', 'Spline Smoothed Energy' , 'Spline Extrema' , 'Energy Extrema')

    



    
    disp('what')

%% End of plot
%%



%     [exte , exti] = extrema(longsyl , 0 , 0);
%     
%     extt = exti/fs;
%     
%     extintrill = find(extt> bulktime(1) & extt < bulktime(2));
%     
%     
%     for k = extintrill(1):extintrill(end)
%         if(exttype(exte , k)==1)
%             syllable_domain = [exti(k-1), exti(k+1)];
% %             peaksi = find_energy_peaks(shortsyl , syllable_domain);
% %             syltime = select_peaks(peaksi , 2)/fs;
%             syltime = findbigpeak(shortsyl , syllable_domain)/fs; 
%             hold (axes, 'on')
% %             syltime = (peaksi(1,1):peaksi(end,end))/fs;
%             plot(axes , syltime , ones(size(syltime))*0.05 , 'm');
%             hold (axes, 'off')
%         end
%     end

function syltime = select_peaks(peaks , threshdb)
    peaks(:,3) = 20*log10(max(peaks(:,3))./peaks(:,3));
    selected = find(peaks(:,3) < threshdb);
    syltime = [peaks(selected(1),1) , peaks(selected(end) , 2)];
    
    
function peakdomain = findbigpeak(env , domain) 
    [maxv, maxi] = extrema(env(domain(1):domain(2)) , 1 , 0);
    maxi = maxi + domain(1) -1;
    [minv, mini] = extrema(env , -1 , 0);
    [maxv , ind] = max(maxv);
    maxi = maxi(ind);
    start = find(mini < maxi , 1 , 'last');
    stop = find(mini > maxi , 1 , 'first');
    peakdomain = mini([start , stop]);
    
    
function peaksi = find_energy_peaks(energy , domain)
    %peaksi is a list 
    % each row is: 
    % [leftmin_index , rightmin_index , max_value]
    [exte , exti] = extrema(energy(domain(1):domain(2)) , 0 , 0);
    exti = exti + domain(1) -1;
    n=1;

    peaksi = zeros(length(exti) , 3);
    for k=1:length(exte)

        if(exttype(exte , k)==1)
            area = sum(exte(k-1:k+1));
            peaksi(n,:) = [exti(k-1); exti(k+1) ; exte(k)];
            n=n+1;
        end

    end

    peaksi = peaksi(1:n-1 , :);


    
function  decision = exttype(extarr , index)
%         return values:
%         max == 1
%         min == -1
    if(index<=1 || index>=length(extarr))
        decision = 0;
    elseif(extarr(index)<extarr(index+1) && extarr(index)<extarr(index-1))
        decision = -1;
    elseif(extarr(index)>extarr(index+1) && extarr(index)>extarr(index-1))
        decision = 1;
    end
    
