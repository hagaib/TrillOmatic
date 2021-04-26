function varargout = trillOmatic(varargin)
% TRILLOMATIC MATLAB code for trillOmatic.fig
%      TRILLOMATIC, by itself, creates a new TRILLOMATIC or raises the existing
%      singleton*.
%
%      H = TRILLOMATIC returns the handle to a new TRILLOMATIC or the handle to
%      the existing singleton*.
%
%      TRILLOMATIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRILLOMATIC.M with the given input arguments.
%
%      TRILLOMATIC('Property','Value',...) creates a new TRILLOMATIC or raises the
%      existing singleton*.  Starting from the left, propeerty value pairs are
%      applied to the GUI before trillOmatic_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trillOmatic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trillOmatic

% Last Modified by GUIDE v2.5 01-Aug-2020 13:25:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trillOmatic_OpeningFcn, ...
                   'gui_OutputFcn',  @trillOmatic_OutputFcn, ...
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




% --- Executes just before trillOmatic is made visible.
function trillOmatic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trillOmatic (see VARARGIN)

% Choose default command line output for trillOmatic
handles.output = hObject;

handles.filenames = get_list_csv('filenames.csv');
fs = 44100;
p.winsize = floor(0.01*fs);
p.overlap = floor(p.winsize/2);
p.wintype = 2;
p.nfft = 2^10;
p.bw = [0 , 10];

handles.spect_params  = p;

handles.global_filter_gui_values = [0, fs/2];

% Update handles structure
guidata(hObject, handles);

% xlim_binder( [handles.axes_time , handles.axes_spect] );
linkaxes( [handles.axes_time , handles.axes_spect , handles.axes_teo] , 'x');

set(handles.uitoolbar1.Children , 'State' , 'off');


% UIWAIT makes trillOmatic wait for user response (see UIRESUME)
% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.
function varargout = trillOmatic_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function handles = gather_harmonic_info (handles)
xx = handles.xx;
fs = handles.fs;
timexx = handles.time;

if(~isfield(handles , 'yin'))
    yin = yin_wrapper(xx , fs);
    handles.yin = yin;
else
    yin = handles.yin;
end

[harmonics , energies]= harmonics_of_signal(xx , fs , yin , 2);
SNR = SNR_estimation(xx , fs , timexx , yin);
disp(['SNR =' , num2str(SNR)]);

handles.harmenergy = energies;
handles.harmonics = harmonics;


function h = plot_harmonic_content(axes ,yin , harmonics)

%Plotting yin and harmonics on Spectrogram
h = plot_yin(yin.f0/1000 , yin.dips , yin.time , axes);
set(h(end),'AutoUpdate','off')
hold (axes , 'on');
% h2 = plot(yin.time, harmonics(:,1)/1000 , 'g' , yin.time, harmonics(:,2)/1000 , 'g');
% h = [h ; h2];
hold (axes , 'off');
for i=1:length(h)
    set(h(i) , 'HitTest' , 'off');
end


function SNR = SNR_estimation(xx , fs , time , yin)

min_segment_time = 0.018;
[periodictime , nonperiodictime , persegs , nonpersegs] = pertime(yin.dips , yin.time , min_segment_time);
Eavper = mean(xx(segments2logical(persegs , time , fs)).^2);
xxnonper = xx(segments2logical(nonpersegs , time , fs));
if(isempty(xxnonper))
    Eavnonper = 0;
else
    Eavnonper = mean(xxnonper.^2);
end

SNR = 10*log10(Eavper/Eavnonper-1);



function handles = clear_prev_file_data(handles)
if(isfield(handles , 'harmonics'))
    handles = rmfield(handles , 'harmonics');
end
if(isfield(handles , 'yin'))
    handles = rmfield(handles , 'yin');
end
if(isfield(handles , 'xx'))
    handles = rmfield(handles , 'xx');
end
if(isfield(handles , 'fs'))
    handles = rmfield(handles , 'fs');
end
if(isfield(handles , 'time'))
    handles = rmfield(handles , 'time');
end

handles.t_indicators = {};
handles.s_indicators = {};



function handles = plot_axes(handles)
%plot a wavefile and the appropriate spectrogram 
xxplot = handles.xx;
if(~handles.checkbox_toggle_baseline.Value)
    xxplot = xxplot - handles.bb;
end


plot_TEO(handles);


handles.hxx = plot(handles.axes_time , handles.time , xxplot);
xlim([0 , handles.time(end)])

handles.spect = plot_spectrogram(handles);
xlim([0 , handles.time(end)])

handles.axes_time.Tag = 'axes_time';
handles.axes_spect.Tag = 'axes_spect';
set(handles.hxx , 'HitTest' , 'off');
handles.axes_time.ButtonDownFcn = @axes_time_ButtonDownFcn;
handles.axes_spect.ButtonDownFcn = @axes_spect_ButtonDownFcn;
% handles.hxx.ButtonDownFcn = @axes_time_ButtonDownFcn;


function plot_TEO(handles)
teo = handles.teo;

plot(handles.axes_teo ,handles.time ,  teo);

function spect = plot_spectrogram(handles)
    p = handles.spect_params;
    m = matfile('wintypes.mat');
    win_list = m.funcs;
    win_f = win_list{p.wintype};
    axes(handles.axes_spect)
    spectrogram( handles.xx , win_f(p.winsize) , p.overlap , p.nfft , handles.fs , 'yaxis');
    [s , f , t , psd] = spectrogram( handles.xx , win_f(p.winsize) , p.overlap , p.nfft , handles.fs , 'yaxis');
    spect.s = s;
    spect.f = f;
    spect.t = t;
    spect.p = psd;
%     surf(t , f/1000 , abs(s) , 'EdgeColor' , 'None' );
%     axis xy; axis tight; colormap(jet); view(0,90);
    colormap gray
    colorbar off;
    %     colorbar('Location' , 'southoutside')

    handles.axes_spect.XLim = handles.axes_time.XLim;
    handles.axes_spect.YLim = p.bw;
    
    % setting color limits according to Slider
%     precent = 0.99*handles.slider_spect_color.Value;
%     clim = color_lim_spectrogram(psd , precent);
    [clim , precent] = auto_lower_lim_spectrogram(spect);
    line([0, handles.time(end)], handles.global_filter_gui_values(1), 'Color', 'g', 'LineWidth', 1);
    line([0, handles.time(end)], handles.global_filter_gui_values(2), 'Color', 'g', 'LineWidth', 1);
    caxis(handles.axes_spect , clim);
    handles.slider_spect_color.Value = precent;
    
    children = handles.axes_spect.Children;
    for i=1:length(children)
        set(children(i) , 'HitTest' , 'off');
    end
    



function handles = load_new_file(hObject, handles , filename)
%load a new sound file and save it in handles object
fullname = filename;
if(iscell(filename))
    fullname = [filename{1} , filename{2}];
    filename = filename{2};
end
[xx , fs] = audioread(fullname);
% xx = resample(xx , 44.1*10^3, fs); fs = 44.1*10^3; 
xx = xx(:,1);
handles.fs = fs;
handles.xx_raw = xx;
handles.xx = xx;
handles.time = ((0:(length(xx)-1))/fs)';
handles.filename = filename;
handles.text_filename.String = filename;
handles.bb = baseline(xx , fs , [] , 1500 , 'fir1'); 




% --- Executes on selection change in listbox_filenames.
function listbox_filenames_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_filenames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_filenames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_filenames


% --- Executes during object creation, after setting all properties.
function listbox_filenames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_filenames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

filenames = get_list_csv('filenames.csv');
numerals = (1:length(filenames))';
filenames = strcat(num2str(numerals), {') '} , filenames);

set(hObject , 'String' , filenames);

% Update handles structure
guidata(hObject, handles);



%%
% segment indicator methods
%%
function trillsegs = check_for_dana_benchmark(fullfilename)
fullfilename = strrep(fullfilename , '.wav' , '.xlsx');
try
    rawdata = xlsread(fullfilename , 'spectrogram');
    timestamps = rawdata(:,2);

    trillsegs = zeros(2 , length(timestamps) );
    k=1;
    while(k*2-1<length(timestamps) && ~isnan(timestamps(k*2-1)))
        trillsegs(:,k) = timestamps(2*k-1:2*k);
        k=k+1;
    end
    trillsegs = trillsegs(:,1:k-1);
catch
    trillsegs = nan;
end

function plot_benchmark_indicators(handles)
% plots the benchmarks annotated by Dana as part of her work
segments = handles.benchmarksegs;
N = size(segments,2);
time_max = handles.time_indicator_height;
ax = handles.axes_time;
axes(ax);
time_line_len = [-time_max , time_max];
time_indicators = cell(2 , N);
spect_indicators = cell(2 , N);

% plot lines
for i=1:N
    t = segments(1,i);
    l = line([t , t] , time_line_len , 'Color' , 'k' , 'LineWidth' , 2);
    time_indicators{1,i} = l;
    
    t = segments(2,i);
    l = line( [t , t] , time_line_len , 'Color' , 'k' , 'LineWidth' , 2);
    time_indicators{2,i} = l;
end

ax = handles.axes_spect;
axes(ax);
for i=1:N
    t = segments(1,i);
    l = line([t , t] , [0 , 10] , 'Color' , 'k' , 'LineWidth' , 2);
    spect_indicators{1,i} = l;
    
    t = segments(2,i);
    l = line( [t , t] , [0 , 10] , 'Color' , 'k' , 'LineWidth' , 2);
    spect_indicators{2,i} = l;
end

function [time_indicators , spect_indicators , handles] = plot_segment_indicators(axesh , segments , time_max , handles)
%plots on vector axesh draggable segment indicators

if(isfield(handles , 't_indicators'))
    indt = handles.t_indicators;
    inds = handles.s_indicators;
    delete([indt{:}]);
    delete([inds{:}]);
end

N = size(segments,2);

ax = axesh(1);
axes(ax);
time_line_len = [-time_max , time_max];
time_indicators = cell(2 , N);
spect_indicators = cell(2 , N);

% plot lines
for i=1:N
    t = segments(1,i);
    l = line([t , t] , time_line_len , 'Color' , 'b' , 'LineWidth' , 2 ...
                 , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    time_indicators{1,i} = l;
    
    t = segments(2,i);
    l = line( [t , t] , time_line_len , 'Color' , 'r' , 'LineWidth' , 2 ...
                 , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    time_indicators{2,i} = l;
end

ax = axesh(2);
axes(ax);
for i=1:N
    t = segments(1,i);
    l = line([t , t] , [0 , 10] , 'Color' , 'b' , 'LineWidth' , 2 ...
                 , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    spect_indicators{1,i} = l;
    
    t = segments(2,i);
    l = line( [t , t] , [0 , 10] , 'Color' , 'r' , 'LineWidth' , 2 ...
                 , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    spect_indicators{2,i} = l;
end

%save in handles
handles.t_indicators = time_indicators;
handles.s_indicators = spect_indicators;


function make_indicators_draggable(time_indicators , spect_indicators , segments , ax)
N = size(segments,2);
    %update data and give drag ability
for i = 1:N
    for j=1:2
        lt = time_indicators{j , i};
        ls = spect_indicators{j , i};
        lt.UserData = ls;
        ls.UserData = lt;
        
        if(j==2)
            xmin = segments(1 , i);
            try
                xmax = segments(1 , i+1);
            catch
                xmax = inf;
            end
        else
            xmax = segments(2 , i);
            try
                xmin = segments(2 , i-1);
            catch
                xmin = 0;
            end
        end
        if getappdata(lt, 'constraint_type') == 'h'
            setappdata(lt, 'constraint_parameters', [xmin xmax])
        else
            draggable(lt , 'h' , [xmin xmax] , @move_twin , 'endfcn' , @(h)end_move(h , ax));
        end
        if getappdata(ls, 'constraint_type') == 'h'
            setappdata(ls, 'constraint_parameters', [xmin xmax])
        else
            draggable(ls , 'h' , [xmin xmax] , @move_twin , 'endfcn' , @(h)end_move(h , ax)); 
        end
    end
end
    

function move_twin(h)
twin = h.UserData(1);
twin.XData = h.XData;
h.Tag = 'move';


function end_move(h , axhandle)
handles = guidata(axhandle);
if ~strcmp(h.Tag, 'move') 
    button = questdlg('Warning!! Indicators will be deleted. Are you sure??!' , 'Delete Syllable Indicator' , 'Cancel');
    if(isequal(button , 'Yes'))
       %delete indicators
       handles = delete_indicators(h , axhandle);
    end
else
    h.Tag = '';
    handles = update_segs_from_indicators(handles);
    guidata(axhandle , handles);
end
make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , axhandle);


function handles = update_segs_from_indicators(handles)
indicators = handles.t_indicators;
segs = zeros(size(indicators));
for i=1:size(indicators , 2)
    seg = [indicators{:,i}];
    segs(1,i) = seg(1).XData(1);
    segs(2,i) = seg(2).XData(1);
end
handles.trillsegs = segs;

function line_click_callback(hObject, eventdata, ax)
button = questdlg('Warning!! Indicators will be deleted. Are you sure??!' , 'Delete Syllable Indicator' , 'Cancel');
if(isequal(button , 'Yes'))
    %delete indicators
    delete_indicators(hObject , ax);
end

function handles = delete_indicators(ind , hObject)
handles = guidata(hObject);
tinds = handles.t_indicators;
sinds = handles.s_indicators;

xdata = ind.XData(1);
N = length(tinds);
for i=1:N
    if(xdata == tinds{1 , i}.XData(1) || ...
        xdata == tinds{2 , i}.XData(1)) 
        disp(['syllable' , num2str(i)]);
        
        delete([tinds{:, i}]);
        delete([sinds{:, i}]);
        
        tinds(: , i:N-1) = tinds(: , i+1:N);
        sinds(: , i:N-1) = sinds(: , i+1:N);
        tinds = tinds(: , 1:N-1);
        sinds = sinds(: , 1:N-1);
        
        handles.t_indicators = tinds;
        handles.s_indicators = sinds;
        handles = update_segs_from_indicators(handles);
        guidata(hObject , handles);
        break;
    end
end


function handles = create_new_indicator_pair(xPos , handles)
tinds = handles.t_indicators;
sinds = handles.s_indicators;

addflag = true;
breakflag = false;
for i=1:size(tinds , 2)
    if(xPos > tinds{1,i}.XData(1) && xPos < tinds{2,i}.XData(1))
        addflag = false;
        f = msgbox('Cannot add syllable inside another syllable!', 'Error');
        break;
    end
    if(xPos < tinds{2 , i}.XData(1))
        breakflag = true;
        break;
    end
end

if(~breakflag)
    i=i+1;
end

delta_t = 0.005;
if(addflag)
    ax = handles.axes_time;
    axes(ax);
    time_line_len = [-handles.time_indicator_height , handles.time_indicator_height];
    t = xPos - delta_t;
    lt1 = line([t , t] , time_line_len , 'Color' , 'b' , 'LineWidth' , 2); %...
        % , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    t = xPos + delta_t;
    lt2 = line( [t , t] , time_line_len , 'Color' , 'r' , 'LineWidth' , 2); %...
        % , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
                
    tinds = [tinds(:,1:i-1) , {lt1 ; lt2} , tinds(:,i:end)];

    axes(handles.axes_spect);
    t = xPos - delta_t;
    ls1 = line([t , t] , [0 , 10] , 'Color' , 'b' , 'LineWidth' , 2); % ...
       % , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));
    t = xPos + delta_t;
    ls2 = line([t , t] , [0 , 10] , 'Color' , 'r' , 'LineWidth' , 2); % ...
       % , 'ButtonDownFcn' , @(hObject , eventdata)line_click_callback(hObject , eventdata , ax));

    sinds = [sinds(:,1:i-1) , {ls1 ; ls2} , sinds(:,i:end)];

end
handles.t_indicators = tinds;
handles.s_indicators = sinds;

make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , handles.axes_time);

%%
%ButtonDown Callbacks
%%

% --- Executes on button press in pushbutton_go.
function pushbutton_go_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = clear_prev_file_data(handles);
filename = handles.filenames{handles.listbox_filenames.Value};
handles = load_new_file(hObject , handles , filename);
handles.teo = zeros(size(handles.time));
handles = plot_axes(handles);

handles = gather_harmonic_info(handles);
handles = delete_yinhandles(handles);
handles.yinhandles = plot_harmonic_content(handles.axes_spect , handles.yin , handles.harmonics);

set(handles.checkbox_toggle_yin , 'Value' , 1)

handles.time_indicator_height = max(abs(handles.xx - handles.bb));
% Update handles structure
guidata(hObject, handles);

assignin('base' , 'handles' , handles);


% --- Executes on button press in pushbutton_openfile.
function pushbutton_openfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = '';
if(isfield(handles , 'path'))
    path = handles.path;
end
[filename , path] = uigetfile({'*.wav' ; '*.*'} , 'Open wave file' , path);
if(filename==0)
    return;
end
handles.path = path;
handles = clear_prev_file_data(handles);
handles = load_new_file(hObject , handles , {path , filename});
handles.teo = zeros(size(handles.time));
handles = plot_axes(handles);

handles = gather_harmonic_info(handles);
handles = delete_yinhandles(handles);
handles.yinhandles = plot_harmonic_content(handles.axes_spect , handles.yin , handles.harmonics);


set(handles.checkbox_toggle_yin , 'Value' , 1)

handles.time_indicator_height = max(abs(handles.xx - handles.bb));

trillsegs = check_for_dana_benchmark([path , filename]);
if(~isnan(trillsegs))
    handles.benchmarksegs = trillsegs;
    plot_benchmark_indicators(handles);
end
% Update handles structure
guidata(hObject, handles);

assignin('base' , 'handles' , handles);


% --- Executes on button press in pushbutton_spect.
function pushbutton_spect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% openfig('spectrogram_options.fig' , 'reuse')
params = spectrogram_options(handles.spect_params);
if(~isempty(params) )
    handles.spect_params = params;
end
plot_spectrogram(handles);
handles = delete_yinhandles(handles);
set(handles.checkbox_toggle_yin , 'Value' , 0)

if(isfield(handles, 'trillsegs') && isfield(handles, 'time_indicator_height'))
[~, ~, handles] = ...
    plot_segment_indicators([handles.axes_time , handles.axes_spect] , handles.trillsegs , handles.time_indicator_height , handles);
end

guidata(hObject , handles);  


% --- Executes on button press in pushbutton_sound.
function pushbutton_sound_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Sound is playing?')
xlim = handles.axes_time.XLim;
player = playsound(handles.xx , handles.fs , xlim);
handles.player = player;
guidata(hObject , handles);


function player = playsound(xx , fs , t)

samples = 1 + floor(t * fs);
samples(1) = max(samples(1) , 1);
samples(2) = min(samples(2) , length(xx) );

xx = xx( samples(1):samples(2) );

player = audioplayer(xx/max(abs(xx)), fs);
play(player);


% --- Executes on mouse press over axes background.
function axes_time_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

button = questdlg('Adding a new Syllable!!!. Are you sure??!' , 'Add Syllable' , 'Cancel');
if(isequal(button , 'Yes'))
    handles = create_new_indicator_pair(hObject.CurrentPoint(1) , handles);
    handles = update_segs_from_indicators(handles);
    make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , handles.axes_time);
    guidata(hObject , handles);
end
assignin('base' , 'handles' , handles);

% --- Executes on mouse press over axes background.
function axes_spect_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

button = questdlg('Adding a new Syllable!!!. Are you sure??!' , 'Add Syllable' , 'Cancel');
if(isequal(button , 'Yes'))
    handles = create_new_indicator_pair(hObject.CurrentPoint(1) , handles);
    handles = update_segs_from_indicators(handles);
    make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , handles.axes_time);
    guidata(hObject , handles);
end


% --- Executes on button press in checkbox_toggle_yin.
function checkbox_toggle_yin_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_toggle_yin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_toggle_yin
value = get(hObject,'Value');

if(~isfield(handles , 'yinhandles'))
    pushbutton_yin_Callback(hObject, eventdata, handles);
else
    h = handles.yinhandles;
    for i=1:length(h)
        if(value ==1)
            set(h(i) , 'Visible' , 'on');
        else
            set(h(i) , 'Visible' , 'off');
        end
    end
end


% --- Executes on button press in pushbutton_reset_plot.
function pushbutton_reset_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.checkbox_toggle_baseline ,'Value');
if(isfield(handles , 'hxx'))
    hxx = handles.hxx;
    delete(hxx);
    axes(handles.axes_time)
    if(value == 1)
        hxx = plot(handles.time , handles.xx , 'b');
    else
        hxx = plot(handles.time , handles.xx-handles.bb , 'b');
    end
    handles.hxx = hxx;
    guidata(hObject , handles);
    
    handles.axes_time.XLim = handles.axes_spect.XLim;
end


% --- Executes on button press in checkbox_toggle_baseline.
function checkbox_toggle_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_toggle_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_toggle_baseline
value = get(hObject,'Value');
if(isfield(handles , 'hxx'))
    hxx = handles.hxx;
    delete(hxx);
    axes(handles.axes_time)
    hold on
    if(value == 1)
        hxx = plot(handles.time , handles.xx , 'b');
    else
        hxx = plot(handles.time , handles.xx-handles.bb , 'b');
    end
    hold off
    handles.hxx = hxx;
    guidata(hObject , handles);
end


function handles = delete_yinhandles(handles)
if(isfield(handles, 'yinhandles'))
    for i=1:length(handles.yinhandles)
        delete(handles.yinhandles(i))
    end
    handles = rmfield(handles , 'yinhandles');
end


% --- Executes on button press in pushbutton_datatable.
function pushbutton_datatable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_datatable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_segs_from_indicators(handles);
params = extract_trill_parameters(handles.xx , handles.fs , handles.yin , handles.trillsegs , handles.trillenv);
handles.trillparams = params;
trill_data(handles.trillparams , handles.trillsegs , handles.filename);
guidata(hObject , handles);






%%
%% Toolbar Methods
%
%%
% --------------------------------------------------------------------
function uitoggletool_zoomin_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isequal(hObject.State ,'on'))
    h = zoom;
    h.Direction = 'in';
    h.Enable = 'on';
else
    zoom off
end
untoggle_all(hObject, handles);

% --------------------------------------------------------------------
function uitoggletool_zoomout_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isequal(hObject.State ,'on'))
    h = zoom;
    h.Direction = 'out';
    h.Enable = 'on';
else
    zoom off
end
untoggle_all(hObject, handles);

% --------------------------------------------------------------------
function uitoggletool_pan_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
untoggle_all(hObject, handles);
if(isequal(hObject.State ,'on'))
    pan on
else
    pan off
end

 
function untoggle_all(hObject , handles)
tools = handles.uitoolbar1.Children;
for i=1:length(tools)
    if(isprop(tools(i) , 'State'))
        if(~isequal(tools(i) , hObject))
            tools(i).State = 'off';
        end
    end
end
pan off
if(~isequal(hObject.Tag , 'uitoggletool_zoomin') && ~isequal(hObject.Tag , 'uitoggletool_zoomout'))
    zoom off
end


% --- Executes on slider movement.
function slider_spect_color_Callback(hObject, eventdata, handles)
% hObject    handle to slider_spect_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = hObject.Value;
v = 0.99*v;

ax = handles.axes_spect;
clim = color_lim_spectrogram(handles.spect.p , v);
caxis(ax , clim);

function clim = color_lim_spectrogram(psd , precent)
logpsd = 10*log10(psd);
minlogpsd = min(min(logpsd));
maxlogpsd = max(max(logpsd));
lerp = (1-precent)*minlogpsd + precent*maxlogpsd;
clim = [lerp , maxlogpsd];

function [clim , precent] = auto_lower_lim_spectrogram(spect)
logpsd = 10*log10(spect.p);


band = spect.f>1800 & spect.f<3500;
logpsd_start = logpsd(band , spect.t<0.1);
noisepsd = mean(max(logpsd_start));
maxlogpsd = max(max(logpsd));
minlogpsd = min(min(logpsd));

clim = [noisepsd, maxlogpsd];
precent = (noisepsd - minlogpsd) / (maxlogpsd - minlogpsd);



% --- Executes during object creation, after setting all properties.
function slider_spect_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_spect_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_indicators.
function checkbox_indicators_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_indicators (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_indicators
if(isfield(handles , 't_indicators'))
    t_indicators = handles.t_indicators;
    state = handles.t_indicators{1}.Visible;
    if(isequal(state ,'on')) , state = 'off'; else, state = 'on'; end
        
    for i=1:size(t_indicators , 2)
        for j=1:2
            l = t_indicators{j,i};
            l.Visible = state;
        end
    end
    
    if(isfield(handles , 's_indicators'))
        s_indicators = handles.s_indicators;
        for i=1:size(s_indicators , 2)
            for j=1:2
                l = s_indicators{j,i};
                l.Visible = state;
            end
        end
    end
end

%% Automatic Segmentation algorithms
% --- Executes on button press in pushbutton_ymfa.
function pushbutton_ymfa_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ymfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[trill_f0band , tpeaks] = area_of_interest4( [] , false , handles.xx , handles.fs);
trill_f0band
[detect , segs , env , smoothteo] = longtrill_syllable_detection41(handles.xx , handles.fs , handles.yin, [], trill_f0band);
handles.trillsegs = segs;
handles.trillenv = env;
handles.teo = smoothteo;

[time_indicators , spect_indicators , handles] = ...
    plot_segment_indicators([handles.axes_time , handles.axes_spect] , segs , handles.time_indicator_height , handles);

make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , handles.axes_time);
plot_TEO(handles);

guidata(hObject , handles);


% --- Executes on button press in pushbutton_aoidetect.
function pushbutton_aoidetect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_aoidetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[trill_f0band , tpeaks] = area_of_interest4( [] , false , handles.xx , handles.fs);
aoidetect = tpeaks;
trill_f0band
[detect , segs , env] = longtrill_syllable_detectionAOI(handles.xx , handles.fs , handles.yin, aoidetect, trill_f0band , false);
% [detect , segs , env] = longtrill_syllable_detectionAOI(handles.xx , handles.fs , handles.yin, aoidetect, [1800 , 3600] , false);
handles.trillsegs = segs;
handles.trillenv = env;

[time_indicators , spect_indicators , handles] = ...
    plot_segment_indicators([handles.axes_time , handles.axes_spect] , segs , handles.time_indicator_height , handles);

make_indicators_draggable(handles.t_indicators , handles.s_indicators , handles.trillsegs , handles.axes_time);
plot_TEO(handles);

guidata(hObject , handles);


%% band width filters

function bw_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bw_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_high_edit as text
%        str2double(get(hObject,'String')) returns contents of bw_high_edit as a double


% --- Executes during object creation, after setting all properties.
function bw_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bw_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bw_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_low_edit as text
%        str2double(get(hObject,'String')) returns contents of bw_low_edit as a double


% --- Executes during object creation, after setting all properties.
function bw_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bp_filter_checkbox.
function bp_filter_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to bp_filter_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bp_filter_checkbox
is_on = get(hObject, 'Value');
handles.bp_filter_on = is_on;
if is_on
    f_low = str2double(handles.bw_low_edit.String) * 1000;
    f_high = str2double(handles.bw_high_edit.String) * 1000;
    xxfilt = filter_signal(handles.xx_raw, handles.fs, f_low, f_high); 
    handles.xx = xxfilt;
    handles.global_filter_gui_values = [f_low, f_high]; 
else
    handles.xx = handles.xx_raw;
end
guidata(hObject, handles);
yin_toggle_value = get(handles.checkbox_toggle_yin , 'Value');
if yin_toggle_value
    pushbutton_yin_Callback(hObject, eventdata, handles);
end
plot_axes(handles);


function xxfilt = filter_signal(xx , fs , flow , fhigh)
%filtering the signal to desired spectrum band
if flow == 0
    if fhigh == fs/2; return; end
    % lowpass
    n = cheb2ord(fhigh/fs*2 , fhigh*1.02/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, fhigh/fs*2 , 'low');
elseif fhigh >= fs/2
    % highpass
    n = cheb2ord(flow/fs*2 , flow*0.98/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, flow/fs*2 , 'high');
else
    n = cheb2ord([flow, fhigh]/fs*2 , [flow*0.98 , fhigh*1.02]/fs*2 , 0.5 , 40);
    [z,p,k] = cheby2(n , 40, [flow , fhigh]/fs*2 , 'bandpass');
end
sos = zp2sos(z,p,k);
xxfilt = sosfilt(sos , xx);

    
   