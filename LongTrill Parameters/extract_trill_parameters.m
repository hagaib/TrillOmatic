function params = extract_trill_parameters(xx , fs , yin , syllable_timestamps , envelope)
%%extracts useful information about long trills. (for biologic research)
% xx = original signal
% fs = sample rate
% f0 = fundamental frequency information
% times = trame time info
%
% returns structure of desired parameters


if(isempty(syllable_timestamps) || sum(sum(isnan(syllable_timestamps)))) ,...
        params = nan; return; end

params.syllable_timestamps = syllable_timestamps;
params.syllable_dur = syllable_timestamps(2,:)-syllable_timestamps(1,:);
params.trill_dur = syllable_timestamps(end) - syllable_timestamps(1);
params.syllable_count = length(params.syllable_dur);
params.syllable_rate = params.syllable_count / params.trill_dur;
params.syllable_dur_avg = mean(params.syllable_dur);
params.syllable_dur_std = std(params.syllable_dur);

params.gaps_dur = syllable_timestamps(1 , 2:end) - syllable_timestamps(2 , 1:end-1);
params.gaps_dur_avg = mean(params.gaps_dur);
params.gaps_dur_std = std(params.gaps_dur);

params.f0 = yin.f0;
params.t = yin.time;

params.syllable_f0max = zeros(size(params.syllable_dur));
params.syllable_f0min = zeros(size(params.syllable_dur));
params.syllable_f0med = zeros(size(params.syllable_dur));
params.syllable_f0mean = zeros(size(params.syllable_dur));
params.syllable_bw = zeros(size(params.syllable_dur));
params.attack_dur = zeros(size(params.syllable_dur));
params.release_dur = zeros(size(params.syllable_dur));

env_max = zeros(size(syllable_timestamps));

syllables = cell(params.syllable_count,1);

time_max = zeros(size(params.syllable_dur));
time_min = zeros(size(params.syllable_dur));

for i=1:params.syllable_count
%     frames = syllable_frames(1,i):syllable_frames(2,i);
    frames = yin.time > syllable_timestamps(1 , i) & yin.time < syllable_timestamps(2 , i);
    w.dips = yin.dips(frames);
    w.f0 = yin.f0(frames);
    w.time = yin.time(frames);
    syllables{i} = w;
    fvect = medfilt1(w.f0(w.dips<0.2),3);
    if(isempty(fvect)) , fvect = nan; end
    [params.syllable_f0max(i) , time_max(i)] = max(fvect);
    [params.syllable_f0min(i) , time_min(i)] = min(fvect);
    params.syllable_f0med(i) = median(w.f0 , 'omitnan');
    params.syllable_f0mean(i) = mean(w.f0 , 'omitnan');
    params.syllable_bw(i) = params.syllable_f0max(i)-params.syllable_f0min(i);
    
    if(~isempty(w.time))
        time_max(i) = w.time(time_max(i));
        time_min(i) = w.time(time_min(i));
    else
        time_max(i) = nan; time_min(i) = nan;
    end
    
    samples = floor(syllable_timestamps(:,i)*fs);
    samples = samples(1):samples(2);
    try
    [env_max(1,i) , env_max(2,i)] = max(envelope(samples));
    catch
%         i
    end
    params.attack_dur(i) = (env_max(2,i)-1)/fs;
    params.release_dur(i) = (length(samples)-env_max(2,i))/fs;
end

params.syllables = syllables;

[~, loudest] = max(env_max(1,:));
maxt = env_max(2,loudest);
maxt = maxt/fs + syllable_timestamps(1 , loudest);

params.trill_attack_dur = maxt - syllable_timestamps(1,1);
params.trill_release_dur = params.trill_dur - params.trill_attack_dur;

params.trill_bw = max(params.syllable_f0max) - min(params.syllable_f0min);
params.syllable_f0max_derivative = params.syllable_f0max(2:end) - params.syllable_f0max(1:end-1);
deltaT = time_max(2:end) - time_max(1:end-1);
params.syllable_f0max_derivative  = params.syllable_f0max_derivative ./ deltaT;

params.syllable_f0min_derivative = params.syllable_f0min(2:end) - params.syllable_f0min(1:end-1);
deltaT = time_min(2:end) - time_min(1:end-1);
params.syllable_f0min_derivative  = params.syllable_f0min_derivative ./ deltaT;