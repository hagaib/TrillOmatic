function env = energy_envelope(energy, fs , time , axes , bulktime)


    env = envelope(energy , floor(fs*0.05),'peak' );
%     smallsyl = envelope(wenv , floor(fs*0.001) , 'peak' );
    longsyl = envelope(energy , floor(fs*0.02) , 'peak' );
    
    
    
%     smallsyl = filter(b,1,wenv);
    L =floor(0.02*fs);
    b = ones(1,L)/L;
    if(mod(L,2)==0) , L=L+1; end
    win = gausswin(L);
    win = win/sum(win);
    shortsyl = filter(win , 1 , energy);
    shortsyl= [shortsyl((L+1)/2:end) ; zeros((L-1)/2 , 1)];
    
    
    
    hold(axes , 'on');
    plot(axes, time , energy*100 , time , env*100 , time , shortsyl*100 , time , longsyl*100);
    legend(axes , 'xx' , 'energy' , 'env 0.05' , 'shortsyl' , 'longsyl')
    hold(axes , 'off');
    
    [exte , exti] = extrema(longsyl , 0 , 0);
    
    extt = exti/fs;
    
    extintrill = find(extt> bulktime(1) & extt < bulktime(2));
    
    
    for k = extintrill(1):extintrill(end)
        if(exttype(exte , k)==1)
            syllable_domain = [exti(k-1), exti(k+1)];
%             peaksi = find_energy_peaks(shortsyl , syllable_domain);
%             syltime = select_peaks(peaksi , 2)/fs;
            syltime = findbigpeak(shortsyl , syllable_domain)/fs; 
            hold (axes, 'on')
%             syltime = (peaksi(1,1):peaksi(end,end))/fs;
            plot(axes , syltime , ones(size(syltime))*0.05 , 'm');
            hold (axes, 'off')
        end
    end

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
    
