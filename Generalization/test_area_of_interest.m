clear
isplot = true;

% for debugging
% foldername = '';
% gt_filename = 'bandwidth_trills.csv';

% for testing
foldername = '../Evaluation/Generalization/Synth Trills/Trill DB - Union/Natural';
gt_filename = 'trill_bandwidth_gt.csv';

gt_path = fullfile(foldername , gt_filename);
gt_list = get_list_csv(gt_path);

filecount = size(gt_list,1);

low_error = nan(filecount , 1);
high_error = nan(filecount , 1);
estimated_bw = nan(filecount , 2);

for i=32:filecount
    fprintf('Processing %d..\n',i);
    audiofilename = fullfile(foldername , gt_list{i , 1} );
    gt_band = [str2num(gt_list{i , 2}) ,  str2num(gt_list{i , 3})]; 
    result_band = area_of_interest4(audiofilename , isplot)*10^-3;
    low_error(i) = gt_band(1) - result_band(1); 
    high_error(i) = gt_band(2) - result_band(2); 
    estimated_bw(i,:) = result_band;
end

figure(202)
subplot(2,1,1)
h = histogram(low_error , 'BinWidth',0.05 , 'Normalization' , 'count');
title('Errors at Bottom Edge of Frequency Band. Positive Errors Mean Over Detection (Good).');
subplot(2,1,2)
histogram(-high_error , 'BinWidth',0.05 , 'Normalization' , 'count');
title('Errors at Top Edge of Frequency Band. Negative Errors Mean Under Detection (Bad)');

disp('Mean bandwidth:');
mean(cellfun(@str2num,gt_list(:,3))-cellfun(@str2num,gt_list(:,2)))
disp('Mean center frequency:');
mean((cellfun(@str2num,gt_list(:,3))+cellfun(@str2num,gt_list(:,2)))/2)
% results.high_error = high_error;
% results.low_error = low_error;
% results.estimated_bw = estimated_bw;
% save('results5.mat' , 'results');