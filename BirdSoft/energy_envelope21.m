function [env , detect , detectsegs] = energy_envelope21(energy, fs , time , axes , bulktime , xx)

    detect = zeros(size(energy));
    detectsegs = [];

    env = envelope(energy , floor(fs*0.05),'peak' );
%     smallsyl = envelope(wenv , floor(fs*0.001) , 'peak' );
    longsyl = envelope(energy , floor(fs*0.02) , 'peak' );
    
    henv = envelope(xx);
    
    
%     smallsyl = filter(b,1,wenv);
    L =floor(0.02*fs);
    if(mod(L,2)==0) , L=L+1; end
    win = gausswin(L);
    win = win/sum(win);
    smoothenergy = filter(win , 1 , energy);
    smoothenergy= [smoothenergy((L+1)/2:end) ; zeros((L-1)/2 , 1)];
    
    
      
    [high , hightime] = extrema(longsyl , 1 , 0);
    [low , lowtime] = extrema(longsyl , -1 , 0);
    lowtime = lowtime/fs;
    hightime = hightime/fs;
    
    hightime = hightime(hightime>bulktime(1) & hightime < bulktime(2));
    lowind = [find(lowtime<hightime(1) , 1 , 'last') , find(lowtime>hightime(end) , 1)];
    lowtime = lowtime(lowind(1):lowind(2));
    
    
    hightime = floor(hightime*fs);
    lowtime = floor(lowtime*fs);
    highenergy = energy(hightime);
    lowenergy = energy(lowtime);
    
    de = energy(2:end) - energy(1:end-1);
    sde = smoothenergy(2:end) - smoothenergy(1:end-1);
    slopethresh = 10^-7;
    
    p = 0.5;
    psoft = 0.4;
    
    lowenergy(1) = 20*log10(lowenergy(1));
    
    for k=1:length(hightime)
        [highenergy(k) , hightime(k)] = max(energy(lowtime(k):lowtime(k+1)));
        hightime(k) = hightime(k) + lowtime(k)-1;
        disp(['high: ' , num2str(hightime(k)/fs)])
        if(k<length(hightime))
            [lowenergy(k+1) , lowtime(k+1)] = min(energy(hightime(k):hightime(k+1)));
            lowtime(k+1) = lowtime(k+1) + hightime(k)-1;
            disp(['low: ' , num2str(lowtime(k+1)/fs)])
        end
        highenergy(k) = 20*log10(highenergy(k));
        lowenergy(k+1) = 20*log10(lowenergy(k+1));
        
        istart = hightime(k);
        iend = istart;
        
        rightthresh = 10^((p*highenergy(k) + (1-p)*lowenergy(k+1))/20);
        leftthresh = 10^((p*highenergy(k) + (1-p)*lowenergy(k))/20);
        softrightthresh = 10^((psoft*highenergy(k) + (1-psoft)*lowenergy(k+1))/20);
        softleftthresh = 10^((psoft*highenergy(k) + (1-psoft)*lowenergy(k))/20);
        
        goflag = true;
        while(goflag)
            istart = istart-1;
            goflag = false;
            if(istart> lowtime(k))
                if(smoothenergy(istart) > leftthresh|| energy(istart)>leftthresh)
                    goflag = true;

                elseif(energy(istart) < leftthresh)
%                     if(energy(istart) > softleftthresh)
%                         goflag = true;
                    if (de(istart)<0 || de(istart) > 10*slopethresh)
                           % || (sde(istart)<0 || sde(istart) > slopethresh*20)
                        goflag = true;
                    end
                end
            end
        end
        goflag = true;
        while(goflag)
            iend = iend+1;
            goflag = false;
            if(iend < lowtime(k+1))
                if(smoothenergy(iend) > rightthresh || energy(iend)>rightthresh)
                    goflag = true;
                elseif(energy(iend) < rightthresh)
                    if(smoothenergy(iend) > softrightthresh || energy(iend) > softrightthresh)
                        goflag = true;
                    elseif (de(iend)>0 || de(iend) < -slopethresh) %%...
                            %%||(sde(iend)>0 || sde(iend) < -slopethresh*20)
                        goflag = true;
                    end
                end
            end
            
            
        end
        detect(istart:iend) = 1;
        detectsegs(:,k) = [istart/fs , iend/fs];
    end
    hold(axes , 'on');
    plot(axes, time , energy*100 , time , env*100 , time , smoothenergy*100 , time , longsyl*100 , time , detect*0.4 , time , henv);
    legend(axes , 'xx' , 'energy' , 'env 0.05' , 'smooth e' , 'longsyl' , 'detect' , 'hilbert')
    hold(axes , 'off');
    
