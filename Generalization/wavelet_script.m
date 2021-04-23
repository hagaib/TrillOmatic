foldername = 'Single Trills';
gt_filename = 'single_trills_bandwidth.csv';

gt_path = fullfile(foldername , gt_filename);
gt_list = get_list_csv(gt_path);


filenum = 1;
filecount = size(gt_list,1);
audiofilepath = fullfile(foldername, gt_list{filenum});
[x, fs] = audioread(audiofilepath);


dwtlevel = 8;
wavetype = 'coif5';

% cw1 = cwt(x,1:32,wavetype,'plot');
% cw2 = cwt(x,1:32,wavetype,'scal');

figure(4)

[c,l]=wavedec(x,dwtlevel,wavetype);
% [c,l] = wavedec(s,3,'db1');
[cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
% Compute and reshape DWT to compare with CWT.
len = length(c);
cfd=zeros(dwtlevel,len);

for k=1:dwtlevel
    d = detcoef(c,l,k);
    cfd(k,:) = imresize(d, size(c), 'nearest');
end
cfd=cfd(:);
I=find(abs(cfd) <sqrt(eps));
cfd(I)=zeros(size(I));
cfd=reshape(cfd,dwtlevel,len);
% Plot DWT.
subplot(311); plot(x); title('Analyzed signal.');
set(gca,'xlim',[0 inf]);
subplot(312); 
image(flipud(wcodemat(cfd,255,'row')));
colormap(pink(255));
set(gca,'yticklabel',[]);
title('Discrete Transform,absolute coefficients');
ylabel('Level');
% Compute CWT and compare with DWT
subplot(313);
ccfs=cwt(x,1:32,wavetype,'plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');