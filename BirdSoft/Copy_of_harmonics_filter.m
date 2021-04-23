function varargout = harmonics_filter(varargin)
% HARMONICS_FILTER MATLAB code for harmonics_filter.fig
%      HARMONICS_FILTER, by itself, creates a new HARMONICS_FILTER or raises the existing
%      singleton*.
%
%      H = HARMONICS_FILTER returns the handle to a new HARMONICS_FILTER or the handle to
%      the existing singleton*.
%
%      HARMONICS_FILTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HARMONICS_FILTER.M with the given input arguments.
%
%      HARMONICS_FILTER('Property','Value',...) creates a new HARMONICS_FILTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before harmonics_filter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to harmonics_filter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help harmonics_filter

% Last Modified by GUIDE v2.5 15-Nov-2017 17:10:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @harmonics_filter_OpeningFcn, ...
                   'gui_OutputFcn',  @harmonics_filter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function handles = resetHandles(handles)
if(isfield(handles , 'detect'))
    handles = rmfield(handles , 'detect');
end
if(isfield(handles , 'envpermed'))
    handles = rmfield(handles , 'envpermed');
end
if(isfield(handles , 'envnonpermed'))
    handles = rmfield(handles , 'envnonpermed');
end
if(isfield(handles , 'SNR_est'))
    handles = rmfield(handles , 'SNR_est');
end
if(isfield(handles , 'xx_env'))
    handles = rmfield(handles , 'xx_env');
end
if(isfield(handles , 'xx_filter'))
    handles = rmfield(handles , 'xx_filter');
end
if(isfield(handles , 'yin'))
    handles = rmfield(handles , 'yin');
end
if(isfield(handles , 'periodictime'))
    handles = rmfield(handles , 'periodictime');
end
if(isfield(handles , 'nonperiodictime'))
    handles = rmfield(handles , 'nonperiodictime');
end
if(isfield(handles , 'xxds'))
    handles = rmfield(handles , 'xxds');
end
if(isfield(handles , 'xx'))
    handles = rmfield(handles , 'xx');
end
if(isfield(handles , 'time'))
    handles = rmfield(handles , 'time');s
end
if(isfield(handles , 'fs'))
    handles = rmfield(handles , 'fs');
end
if(isfield(handles , 'output'))
    handles = rmfield(handles , 'output');
end




% --- Executes just before harmonics_filter is made visible.
function harmonics_filter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to harmonics_filter (see VARARGIN)

% Choose default command line output for harmonics_filter
handles.output = hObject;
handles.xx = varargin{1};
handles.fs = varargin{2};
filename = varargin{3};
yin = varargin{4};
harmonics = varargin{5};
handles.text_filename.String = filename;

[detect , info] = longtrill_syllable_detection(handles.xx , handles.fs , yin , harmonics)

[xxds , dstime , new_fs] = downsample2(handles.xx , handles.fs);


min_segment_time = 0.018;
[periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
    pertime(yin.dips , yin.time , min_segment_time , yin.f0,[2000 , 3500]);

[minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin);


handles.time = dstime;
handles.fs = new_fs;
handles.xxds = xxds;  
handles.periodictime = segments2logical(persegs , handles.time , handles.fs);
handles.nonperiodictime = segments2logical(nonpersegs , handles.time , handles.fs);
handles.yin = yin;


fbounds = [minf0 , maxf0; minf1 , maxf1 ; minf2, maxf2];
axes = [handles.axes_f0 , handles.axes_f1 , handles.axes_f2];

for filternum=1:3
    handles = process_and_plot_filter(fbounds(filternum,2) , fbounds(filternum,1) , handles , filternum);
    
    %plot
    plot(axes(filternum) , handles.time, handles.xx_filter{filternum} , handles.time, handles.xx_env{filternum});
    % figure(filternum)
    % spectrogram(xx_f, 512 , 256 , [] , fs , 'yaxis')
end


axes = handles.axes_resampled;
envpermed = handles.envpermed;
envnonpermed = handles.envnonpermed;
SNR_est = handles.SNR_est;
periodictime = handles.periodictime;
xx_filter = handles.xx_filter;

c = calc_env_coeffs(envpermed , envnonpermed, SNR_est);
[detect , movmed , wenv , low , energy] = longtrill_final_detect(handles.xx_env , envpermed , envnonpermed ,...
    SNR_est, periodictime , new_fs , c , xx_filter , yin);


%GUI and plot
update_edittextfields(handles , minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2);
set_env_coeffs_in_gui(c , handles);
%plot and handle gui
highpass_plot(1000 , handles);
hold(axes , 'on');
plot(axes , dstime , wenv ,dstime , energy, 'c' , dstime , low , 'r' , dstime , movmed , 'g' , dstime , detect*0.5 ,'m');
hold(axes , 'off');


xlim_binder( [handles.axes_resampled , handles.axes_f0 , handles.axes_f1 , handles.axes_f2 , handles.axes_spect] );
handles.axes_resampled.XLim = [0 , length(xxds)/new_fs];
% UIWAIT makes harmonics_filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.detect = detect;
% Update handles structure
guidata(hObject, handles);


function [detect , movmedenvper , wenv , low , energy] = longtrill_final_detect(xx_env , envpermed , envnonpermed , SNR_est , periodictime , fs , c , xx_filter , yin)
%xx_env: cell with 3 spots , each holding envelope of f0 and the first 2
%harmonics
% envpermed: 3-array containing median of envelope for periodic parts of
% signal
%envnonpermed: 3-array containing median of envelope for non-periodic parts
%of signal
%SNR_est: SNR estimation of each of the envelopes' bandwidths
%periodictime: vector of logicals for periodic samples
%fs : sample rate
%c: optional argument, weight coefficients for each envelope. if not
%supplied they are calculated inside the function

if(nargin<7)
    c = calc_env_coeffs(envpermed , envnonpermed, SNR_est);
    disp('calculating coefficients');
end

wenv= weighted_env(xx_env , c);

%high and low bound calculations
high = sum(envpermed .* c);
low = sum(envnonpermed .* c);
% p = 0.15 + 0.05*(sum(c>0));
p=0.2;
if(sum(c>0)==3)
    p=0.25;
end
low = (1-p)*low + p*high;


weighted_bands = weighted_env(xx_filter , c);
energy = hamming_energy(weighted_bands , fs ,  0.01 , 1);
movmedenvper = energy;

%median of periodic envelope
% movmedenvper = wenv;
movmedenvper(~periodictime) = nan;
[~ , badprecent] = mean_yinprobs_trill(yin);
windur = 0.02;
if(badprecent > 0.02)
    movmedenvper = medfilt1(movmedenvper , floor(0.2*fs) , 'omitnan');
    low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
else
    %tested for "fast long trills" , or trills with high lobes between
    %syllables
    disp(['fast trill , ' , num2str(badprecent)])
    movmedenvper = prcfilt(movmedenvper , floor(0.20*fs) , 75 , 'omitnan');
    low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
    
%     movmedenvper = prcfilt(movmedenvper , floor(0.10*fs) , 85 , 'omitnan');
%     low = medfilt1(energy , floor(windur*4*fs) , 'omitnan');
end
movmedenvper(movmedenvper==0) = nan;

high = movmedenvper;

% low = wenv;
% low(periodictime) = nan;
% low = medfilt1(low, floor(0.15*fs) , 'omitnan');
% low = low*ones(size(wenv));


% segs = logical2segments(periodictime , fs);
% segs_between = nan(size(segs) + [0,1]);
% segs_between(1 ,2:end) = segs(2 , 1:end) + 1/fs;
% segs_between(2 ,1:end-1) = segs(1 , 1:end) - 1/fs;
% segs_between(1,1) = 0;
% segs_between(2,end) = (length(wenv)-1)/fs;

% low = nan(size(periodictime));
% for i=1:length(segs_between)
%     between = segments2logical( segs_between(:,i),[0,(length(wenv)-1)/fs],fs);
%     low(between) = prctile(wenv(between) , 30);
% end
% 
% for i=1:length(segs)
%     pstart = floor(segs(1,i)*fs);
%     pend = ceil(segs(2,i)*fs)+1;
%     L = pend-pstart+2;
%     t = 1/L:1/L:(L-1)/L;
%     vstart = low(pstart-1);
%     vend = low(pend+1);
%     low(pstart:pend) = vstart * (1-t) + vend*t;
% end



[detect , low] = syllable_detection(weighted_bands , wenv, low , high ,periodictime , fs , energy);
% detect = cut_env_tails(weighted_bands , fs, energy , detect);
 detect = syllable_detection_postproc(weighted_bands , fs , yin, energy ,high, low , detect);






function [xx_f , env , menvper , envpermed , envnonpermed , SNR_est  , low , high] = process_filtered_signal(xx , fs , lower , higher , periodictime , nonperiodictime)

%filtering the signal to desired spectrum band
wn = [lower , higher]/(fs/2);
wn(end) = min(wn(end) , 0.9);
[b,a] = butter(11 , wn , 'bandpass');
% [b , ~] = fir1(72 , wn , 'bandpass'); a=1;
xx_f = filtfilt(b , a , xx);

%calculating envelope
env = calc_envelope(xx_f);

%median of periodic envelope
menvper = env;
menvper(~periodictime) = nan;
menvper = medfilt1(menvper , floor(0.20*fs) , 'omitnan');

%SNR estimate
envpermed = median(env(periodictime)); %average periodic energy
envnonpermed = median(env(nonperiodictime)); %average nonperiodic energy
SNR_est = 20*log10(envpermed/envnonpermed-1);

%high and low bound calculations
high = envpermed;
low = envnonpermed + (envpermed - envnonpermed) * 0.2;

%median of periodic envelope
movmedenvper = env;
movmedenvper(~periodictime) = nan;
movmedenvper = medfilt1(movmedenvper , floor(0.20*fs) , 'omitnan');
movmedenvper(movmedenvper==0) = nan;

high = movmedenvper;
low = low *ones(size(env));


% energy = hamming_energy(xx_f, fs ,  0.02 , 1);
% 
% detect = syllable_detection(xx ,env , low , high , periodictime ,fs , energy);


function detect = cut_env_tails(xx , fs, energy , detect)

segs = logical2segments(detect , fs);

de = [0 ; energy ; 0];
de = 0.5*fs*(de(3:end)-de(1:end-2)); %centralized difference scheme

segcount = size(segs , 2);
for i=1:segcount
    istart = floor(segs(1,i)*fs)+1;
    iend = floor(segs(2,i)*fs)+1;
    e = energy(istart:iend);
    
%     figure(5)
%     plot(istart:iend , e , 'b')
%     hold on
    
    slopethresh = 15;
    while(abs(de(istart)) < slopethresh)
%         de(istart)
        istart=istart+1;
    end
    while(abs(de(iend)) < slopethresh) , iend=iend-1; end
    
%     plot(istart:iend,energy(istart:iend) , 'g')
%     hold off
    
end

function detect = syllable_detection_postproc(xx , fs , yin, energy , high, low , detect)
%look for undetected syllables in the middle of trill
detectsegs = logical2segments(detect, fs);
seglen = detectsegs(2,:) - detectsegs(1,:);
rhythm = detectsegs(1,2:end) - detectsegs(1,1:end-1);
meanr = mean(rhythm);
stdr = std(rhythm);
gaps = abs(rhythm-meanr)> stdr*2;
for i=1:length(gaps)
    if(gaps(i)>0)
        t = (detectsegs(1,i)+detectsegs(1,i+1))*0.5;
        meanl = mean(seglen(i:i+1));
        newseg = [t , t+meanl];
        newsegyin = yin.time>newseg(1) & yin.time<newseg(2);
        if(mean(yin.dips(newsegyin)) < 0.2 && std(yin.dips(newsegyin)) < 0.15)
            disp(['new syllable spotted:' , num2str(newseg(1)) , ':', num2str(newseg(2))]);
            indices = floor([detectsegs(2,i),detectsegs(1,i+1)].*fs);
            indices = indices(1):indices(2);
            detectsyl = syllable_detection(xx(indices) , xx(indices) , low(indices) , high(indices)/4 , [] , [] , energy(indices));
            detect(indices) = detectsyl;
        end
    end
    
end


%look for undetected syllables at the end of trill
tryflag = true;
while(tryflag)
    t = detectsegs(1,end)+rhythm(end); %approximate start of last missing syllable
    newseg = [t , t+seglen(end)];
    newsegyin = yin.time>newseg(1) & yin.time<newseg(2);
    if(mean(yin.dips(newsegyin)) < 0.3 && std(yin.dips(newsegyin)) < 0.3)
        disp(['new syllable spotted:' , num2str(newseg(1)) , ':', num2str(newseg(2))]);
        indices = floor([detectsegs(2,end),detectsegs(2,end)+meanr+2*stdr].*fs);
        indices(2) = min(indices(2) , length(xx));
        indices = indices(1):indices(2);
        highnan = find(isnan(high(indices)));
        if(~isempty(highnan)) %if high(indices) contains nan
            high(indices(highnan)) = high(indices(highnan(1))-1);
        end
            
        detectsyl = syllable_detection(xx(indices) , xx(indices) , low(indices) , high(indices)/4 , [] , [] , energy(indices));
        detect(indices) = detectsyl;
        syltime = indices(detectsyl>0);
        if(isempty(syltime)) , tryflag = false; 
        else
            syltime = [syltime(1) ; syltime(end)]/fs;
            detectsegs = [detectsegs , syltime];
        end
    else
        tryflag = false;
    end
end


function [detect , lowbound] = syllable_detection(xx , env , lowbound, high , periodictime , fs , energy)
if(size(env) ~= size(high))
    high = high';
end

abovemed = energy > high;


%peak detection
% [peaks , ipeaks] = extrema(env , 1 , 0);
[peaks , ipeaks] = extrema(energy , 1 , 0);

%give periodic peaks a detection bonus
% adjusted_high = high - (high-low)/3*periodictime(ipeaks)';
% adjusted_high = menvper;


% ipeaks = ipeaks(peaks > adjusted_high);

ipeaks = ipeaks(abovemed(ipeaks));


detect = zeros(size(env));


% energy considerations
% windur = 0.02;
% step = 1;
% energy = hamming_energy(xx , fs ,  windur , step);


for i=1:length(ipeaks)
    istart=ipeaks(i);
    iend = ipeaks(i);
    while(istart>1 && energy(istart)>=lowbound(istart))
        istart=istart-1;
    end
    while(iend<length(env) && energy(iend)>=lowbound(iend))
        iend = iend+1;
    end
    
    
    if(lowbound(istart) > lowbound(iend))
        tightbound = lowbound(istart);
        while(energy(iend) < tightbound)
            iend = iend-1;
        end
    else
        tightbound = lowbound(iend);
        while(energy(istart) < tightbound)
            istart = istart+1;
        end
    end
        
    if(istart~=iend)
        detect(istart:iend) = 1;
    end
end



%using area of lobe
% for i=1:length(ipeaks)
% %     disp(['--peak: ' , num2str(ipeaks(i))]);
%     istart=ipeaks(i)-1;
%     iend=ipeaks(i)+1;
%     
%     larea = env(ipeaks(i)); %lobe area 
%     flag = 1;
%     while(flag)
%         flag=0;
% %         disp(['start - ' , num2str(istart ),' : ' ,num2str(env(istart)/larea)]);
% %         disp(['end - ' , num2str(iend ),' : ' ,num2str(env(iend)/larea)]);
%         slope = env(iend+1)-env(iend);
%         if(env(iend)/larea > 10^-3 || slope < -0.01)
%             iend = iend+1;
%             larea = larea + env(iend);
%             flag=1;
%         end
%         slope = env(istart-1)-env(istart);
%         if(env(istart)/larea > 10^-3 || slope < -0.01)
%             istart = istart-1;
%             larea = larea + env(istart);
%             flag=1;
%         end
%     end
% %     figure(5)
% %     plot((istart:iend)/fs , env(istart:iend));
% %     pause(0.1)
%     
%     if(iend-istart>3)
%         detect(istart:iend)=1;
%     end
% end



%using given lower bound
%changed on 14/11/17
% for i=1:length(ipeaks)
%     istart=ipeaks(i);
%     iend = ipeaks(i);
%     while(istart>1 && env(istart)>=low(istart))
%         istart=istart-1;
%     end
%     while(iend<length(env) && env(iend)>=low(iend))
%         iend = iend+1;
%     end
%     
%     if(istart~=iend)
%         detect(istart:iend) = 1;
%     end
% end


function env = calc_envelope(xx)

% env = envelope(xx);
env = envelope(xx , 200, 'peaks');

%changed on 3/11/17
% b = fir1(156 , 0.05);
% delay = mean(grpdelay(b));

% env = filter(b , 1 , env);

% env = medfilt1(env , 90);

% b = fir1(340 , 0.01);
% delay = delay + mean(grpdelay(b));

% env = filter(b , 1 , env);

b = fir1(156 , 0.02);
delay = mean(grpdelay(b));
env = filter(b , 1 , env);

% figure(100)
% freqz(b,1)

env = [env(delay+1:end); zeros(delay,1)];







function env_out = weighted_env(envs , coeffs)
env_out = coeffs(1) * envs{1} + coeffs(2) * envs{2} + coeffs(3) * envs{3};

function c = calc_env_coeffs(emedper , emednonper , snr_est)

for i=1:length(snr_est)
    if(isreal(snr_est(i))) , snr_est(i) = real(snr_est(i)); 
    else snr_est(i) = 0;
    end
end
snr_est(snr_est<0) = 0;
snr_factor = snr_est;
snr_factor(snr_factor<2) = 0;
%changed on 6/11
% snr_factor(snr_factor~=0) = 1;
snr_factor = snr_factor/snr_factor(1); %snr of 30 or above is cosidered fully accurate
snr_factor

emedper = emedper(1)./emedper;
c = emedper .* snr_factor;
c=[1;0;0];



% --- Outputs from this function are returned to the command line.
function varargout = harmonics_filter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = logical2segments(handles.detect , handles.fs);


function [dsxx , dstime , dsfs] = downsample2(xx , fs)
downsample_factor = 2;
dsfs = fs/downsample_factor;
dsxx= resample(xx , 1 , downsample_factor);
dstime = (0:length(dsxx)-1)/dsfs;

function [minf , maxf] = bw_boundaries (fvect , selection , maxbw , p)
%fvect = vector of frequencies
%selection = vector of 1's and 0's (indices chosen)
%p = (optional) precentile above minimum 

if (nargin < 4) , p=0; end

perf = fvect(selection);
minf = prctile(perf , p);
maxf = prctile(perf , 100-p);

[minf , maxf] = clip_to_min_bw(minf , maxf , maxbw);

function [minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin)


[minf0 , maxf0] = bw_boundaries (yin.f0 , periodictimeyin, 800, 10);
[minf1 , maxf1] = bw_boundaries (harmonics(:,1) , periodictimeyin, 800, 0);
[minf2 , maxf2] = bw_boundaries (harmonics(:,2) , periodictimeyin, 800, 0);
if(isnan(minf2)) , minf2 = 3*minf0; end
if(isnan(maxf2)) , maxf2 = 3*maxf0; end
% if(maxf0 > 3300)
%     maxf0=3300;
%     if(minf0>2000)
%         minf0 = 2000;
%     end
% end

function update_edittextfields(handles , minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2)
handles.edit_lower0.String = num2str(floor(minf0));
handles.edit_higher0.String = num2str(ceil(maxf0));
handles.edit_lower1.String = num2str(floor(minf1));
handles.edit_higher1.String = num2str(ceil(maxf1));
handles.edit_lower2.String = num2str(floor(minf2));
handles.edit_higher2.String = num2str(ceil(maxf2));

function highpass_plot(cutoff , handles)
xxds = handles.xxds;
fs = handles.fs;
yin = handles.yin;
b = fir1(56 , cutoff/fs*2 , 'high');
delay = mean(grpdelay(b,1));
xx_plot = filter(b, 1 , xxds);
xx_plot = [xx_plot(delay+1:end) ; zeros(delay , 1)];

plot(handles.axes_resampled , handles.time , xx_plot);

axes(handles.axes_spect)
spectrogram(xx_plot , 256 , 200 , [] , handles.fs , 'yaxis')
xlabel('');
ylabel('');
colorbar off

A = max(xxds);
plot_yin(yin.f0/1000 , yin.dips , yin.time , gca)
hold(gca , 'on');
plot(gca, handles.time , A*handles.periodictime ,'g' , handles.time , A*handles.nonperiodictime , 'r', 'LineWidth' , 2);
hold(gca , 'off');







function handles = process_and_plot_filter(higher , lower , handles , filternum)    
% Not a callback!
if(~isfield(handles , 'xx_filter')) , handles.xx_filter = []; end
if(~isfield(handles , 'xx_env')) , handles.xx_env = []; end
if(~isfield(handles , 'SNR_est')) ,handles.SNR_est = nan(3,1); end
if(~isfield(handles , 'envpermed')) ,  handles.envpermed = nan(3,1); end
if(~isfield(handles , 'envnonpermed')) , handles.envnonpermed = nan(3,1); end

fs = handles.fs;
xx = handles.xxds;
periodictime = handles.periodictime;
nonperiodictime = handles.nonperiodictime;



[xx_f , env ,menvper, envpermed , envnonpermed , SNR_est , low , high] = ...
    process_filtered_signal(xx , fs , lower , higher , periodictime , nonperiodictime);


handles.xx_filter{filternum} = xx_f;
handles.xx_env{filternum} = env;
handles.SNR_est(filternum) = SNR_est;
handles.envpermed(filternum) = envpermed;
handles.envnonpermed(filternum) = envnonpermed;

disp([num2str(filternum-1) ') "energy":' num2str(sum(env(periodictime))) , ' max:' , num2str(prctile(env(periodictime),95))...
    ' envpermed:' , num2str(envpermed)]);

update_snr_textlabel(SNR_est , filternum , handles);






function set_env_coeffs_in_gui(c , handles)
handles.edit_coeff0.String = num2str(c(1));
handles.edit_coeff1.String = num2str(c(2));
handles.edit_coeff2.String = num2str(c(3));


function c = get_env_coeffs(handles)
c0 = str2double(handles.edit_coeff0.String);
c1 = str2double(handles.edit_coeff1.String);
c2 = str2double(handles.edit_coeff2.String);
c = [c0 , c1 , c2];
    


function edit_lower0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower0 as text
%        str2double(get(hObject,'String')) returns contents of edit_lower0 as a double


% --- Executes during object creation, after setting all properties.
function edit_lower0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_higher0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_higher0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_higher0 as text
%        str2double(get(hObject,'String')) returns contents of edit_higher0 as a double


% --- Executes during object creation, after setting all properties.
function edit_higher0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_higher0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_coeff0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coeff0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coeff0 as text
%        str2double(get(hObject,'String')) returns contents of edit_coeff0 as a double


% --- Executes during object creation, after setting all properties.
function edit_coeff0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coeff0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_coeff1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coeff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coeff1 as text
%        str2double(get(hObject,'String')) returns contents of edit_coeff1 as a double


% --- Executes during object creation, after setting all properties.
function edit_coeff1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coeff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_coeff2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coeff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coeff2 as text
%        str2double(get(hObject,'String')) returns contents of edit_coeff2 as a double


% --- Executes during object creation, after setting all properties.
function edit_coeff2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coeff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_filter0.
function pushbutton_filter0_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_filter0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lower = str2double(handles.edit_lower0.String);
higher = str2double(handles.edit_higher0.String);
handles = process_and_plot_filter(higher , lower , handles , 1);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_filter1.
function pushbutton_filter1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_filter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lower = str2double(handles.edit_lower1.String);
higher = str2double(handles.edit_higher1.String);
handles = process_and_plot_filter(higher , lower , handles , 2);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_filter2.
function pushbutton_filter2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_filter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lower = str2double(handles.edit_lower2.String);
higher = str2double(handles.edit_higher2.String);
handles = process_and_plot_filter(higher , lower , handles , 3);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_soundf0.
function pushbutton_soundf0_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_soundf0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
play_filtered_signal(1 , handles);

% --- Executes on button press in pushbutton_soundf1.
function pushbutton_soundf1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_soundf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
play_filtered_signal(2 , handles);

% --- Executes on button press in pushbutton_soundf2.
function pushbutton_soundf2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_soundf2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
play_filtered_signal(3 , handles);

function play_filtered_signal(filternum , handles)
switch (filternum)
    case 0
        ax = handles.axes_resampled;
    case 1
        ax = handles.axes_f0;
    case 2
        ax = handles.axes_f1;
    case 3
        ax = handles.axes_f2;
end

xx_toplay = handles.xx_filter{filternum};

xlim = ax.XLim;
xlim_samples = 1+floor(xlim * handles.fs);
xlim_samples(1) = max(xlim_samples(1) , 1);
xlim_samples(2) = min(xlim_samples(2) , length(xx_toplay));

soundsc(xx_toplay(xlim_samples(1):xlim_samples(2)) , handles.fs);

function edit_lower1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower1 as text
%        str2double(get(hObject,'String')) returns contents of edit_lower1 as a double


% --- Executes during object creation, after setting all properties.
function edit_lower1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_higher1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_higher1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_higher1 as text
%        str2double(get(hObject,'String')) returns contents of edit_higher1 as a double


% --- Executes during object creation, after setting all properties.
function edit_higher1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_higher1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lower2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower2 as text
%        str2double(get(hObject,'String')) returns contents of edit_lower2 as a double


% --- Executes during object creation, after setting all properties.
function edit_lower2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_higher2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_higher2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_higher2 as text
%        str2double(get(hObject,'String')) returns contents of edit_higher2 as a double


% --- Executes during object creation, after setting all properties.
function edit_higher2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_higher2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_snr_textlabel(snrval , filternum , handles)

switch (filternum)
    case 1
        textlabel = handles.text_snr0;
    case 2
        textlabel = handles.text_snr1;
    case 3
        textlabel = handles.text_snr2;
end
textlabel.String = ['SNR: ' ,num2str(snrval)];


% --- Executes on button press in pushbutton_plotenv.
function pushbutton_plotenv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotenv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


axes = handles.axes_resampled;
time = handles.time;
fs = handles.fs;
envpermed = handles.envpermed;
envnonpermed = handles.envnonpermed;
SNR_est = handles.SNR_est;
periodictime = handles.periodictime;
xx_filter = handles.xx_filter;

c = get_env_coeffs(handles)';
% c = calc_env_coeffs(envpermed , envnonpermed, SNR_est);
[detect , movmed , wenv , low] = longtrill_final_detect(handles.xx_env , envpermed , envnonpermed ,...
    SNR_est, periodictime , fs , c , xx_filter);

%plot and handle gui
set_env_coeffs_in_gui(c , handles);
hold(axes , 'on');
plot(axes , time , wenv , time , low , 'r' , time , movmed , 'g' , time , detect*sum(envpermed) ,'m');
hold(axes , 'off');


% --- Executes on button press in pushbutton_energy.
function pushbutton_energy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xx_filter = handles.xx_filter;
fs = handles.fs;
time = handles.time;
windur = 0.03;
step = nan;
energy = hamming_energy(xx_filter{1} , fs ,  windur , step);
lowbound = prcfilt(energy , floor(windur*5*fs) , 30 , 'omitnan');

ax = handles.axes_resampled;
axes(ax);
hold on;
plot(ax , time , energy , time , lowbound);
hold off;


function energy = hamming_energy(xx , fs ,  windur , step)
energy = zeros(size(xx));
winsize = floor(windur*fs);
win = hamming(winsize);
xx = xx(:);
L = length(xx);
transient = ceil(winsize/2);
xx = [zeros(transient , 1) ; xx ; zeros(transient , 1)];
for i=1:L
    xxwin = xx(i:i+winsize-1);
    e = win .* xxwin;
    e = sum(e.^2);
    energy(i) = e;
end