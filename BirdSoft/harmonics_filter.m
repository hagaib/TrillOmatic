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

[detect , info] = longtrill_syllable_detection(handles.xx , handles.fs , yin , harmonics);



handles.time = info.time;
handles.fs = info.fs;
handles.xxds = info.xxds;  
handles.periodictime = info.periodictime;
handles.nonperiodictime = info.nonperiodictime;
handles.yin = yin;
handles.detect = detect;
handles.info = info;
dstime = info.time;


axes = [handles.axes_f0 , handles.axes_f1 , handles.axes_f2];

for filternum=1:3
    
    %plot
    plot(axes(filternum) , info.time, info.fbands(filternum).xx , info.time, info.fbands(filternum).env);
    % figure(filternum)
    % spectrogram(xx_f, 512 , 256 , [] , fs , 'yaxis')
end


axes = handles.axes_resampled;

%GUI and plot
update_edittextfields(handles , info.fbounds);
set_env_coeffs_in_gui(info.c , handles);
%plot and handle gui
highpass_plot(1000 , handles);
hold(axes , 'on');
plot(axes , dstime , info.wenv ,dstime , info.energy*100, 'c' , dstime , info.lowthresh*100 , 'r' , dstime , info.highthresh*100 , 'g' , dstime , detect*0.5 ,'m');
hold(axes , 'off');
legend(axes , {'signal' , 'Energy' , 'Envelope' , 'lower threshold' , 'high threshold' , 'detection'});


xlim_binder( [handles.axes_resampled , handles.axes_f0 , handles.axes_f1 , handles.axes_f2 , handles.axes_spect] );
handles.axes_resampled.XLim = [0 , length(info.xxds)/info.fs];
% UIWAIT makes harmonics_filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% Update handles structure
guidata(hObject, handles);












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















% --- Outputs from this function are returned to the command line.
function varargout = harmonics_filter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = logical2segments(handles.detect , handles.fs);
varargout{3} = handles.info;




function update_edittextfields(handles ,fbounds)
handles.edit_lower0.String = num2str(floor(fbounds(1,1)));
handles.edit_higher0.String = num2str(ceil(fbounds(1,2)));
handles.edit_lower1.String = num2str(floor(fbounds(2,1)));
handles.edit_higher1.String = num2str(ceil(fbounds(2,2)));
handles.edit_lower2.String = num2str(floor(fbounds(3,1)));
handles.edit_higher2.String = num2str(ceil(fbounds(3,2)));

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


