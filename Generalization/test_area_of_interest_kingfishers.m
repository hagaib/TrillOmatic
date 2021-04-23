clear
gt_filename = 'kingfisher_filenames.csv';
isplot = true;
gt_path = fullfile('..','Algorithm Evaluation','Synth Trills','HNM Trills','Synth HNM' ,'Clean Final');
audio_path = fullfile('..','Algorithm Evaluation','Synth Trills','HNM Trills','Synth HNM' ,'White Noise Final' , 'SNR10db');

gt_list = get_list_csv(gt_filename);

filecount = size(gt_list,1);

low_error = nan(filecount , 1);
high_error = nan(filecount , 1);
estimated_bw = nan(filecount , 2);

for i=8%1:filecount
    fprintf('Processing %d..\n',i);
    audiofilename = fullfile(audio_path , gt_list{i});
    gt_filename = fullfile(gt_path , gt_list{i});
    gt_filename = strrep(gt_filename , '.wav' , '.mat');
    load(gt_filename);
    trill_data = params;
    %% GT bandwidth from params struct
    find_start = 1;
    is_syl = false(size(params.f0));
    for j=1:size(params.segs,2)
        s = trill_data.segs(:,j);
        indices = find_start-1 + ...
            find(trill_data.time(find_start:end) > s(1) & trill_data.time(find_start:end) < s(2));
        is_syl(indices) = true;
        find_start = indices(end)+1;
    end
    gt_band = [min(trill_data.f0(is_syl)) ,  max(trill_data.f0(is_syl))]*10^-3;
    
    [x , fs] = audioread(audiofilename);
    gt_params = extract_trill_parameters(x , fs , trill_data.yin , trill_data.segs , []);
    
    
    result_band = area_of_interest4(audiofilename , isplot)*10^-3;
    low_error(i) = gt_band(1) - result_band(1); 
    high_error(i) = gt_band(2) - result_band(2); 
    estimated_bw(i,:) = result_band;
end

figure(101)
subplot(2,1,1)
h = histogram(low_error , 'BinWidth',0.05 , 'Normalization' , 'count');
title('Errors at Bottom Edge of Frequency Band. Positive Errors Mean Over Detection (Good).');
subplot(2,1,2)
histogram(-high_error , 'BinWidth',0.05 , 'Normalization' , 'count');
title('Errors at Top Edge of Frequency Band. Negative Errors Mean Under Detection (Bad)');