function varargout = birdsong_model_0(varargin)
% BIRDSONG_MODEL0 MATLAB code for birdsong_model0.fig
%      BIRDSONG_MODEL0, by itself, creates a new BIRDSONG_MODEL0 or raises the existing
%      singleton*.
%
%      H = BIRDSONG_MODEL0 returns the handle to a new BIRDSONG_MODEL0 or the handle to
%      the existing singleton*.
%
%      BIRDSONG_MODEL0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIRDSONG_MODEL0.M with the given input arguments.
%
%      BIRDSONG_MODEL0('Property','Value',...) creates a new BIRDSONG_MODEL0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before birdsong_model0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to birdsong_model0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help birdsong_model0

% Last Modified by GUIDE v2.5 24-Jul-2017 16:26:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @birdsong_model_0_OpeningFcn, ...
                   'gui_OutputFcn',  @birdsong_model_0_OutputFcn, ...
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


% --- Executes just before birdsong_model0 is made visible.
function birdsong_model_0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to birdsong_model0 (see VARARGIN)



set(hObject , 'toolbar' , 'figure');

% Choose default command line output for birdsong_model0
handles.output = hObject;
filename = get(handles.filename_editText, 'String');

% Update handles structure
guidata(hObject, handles);

% [xx,fs]=audioread(filename);
% handles.fs = fs;
% handles.xx = xx;
handles.filename = filename;

% load_file(hObject , handles,handles.filename);
% handles = guidata(hObject);

handles.xx_color = 'b';
handles.yy_color = 'r';

% axes1_plot(handles);

handles.stepSize = str2double(get(handles.stepSize_editText, 'String'));
handles.winDur = str2double(get(handles.winDur_editText, 'String'));
handles.noiseCutoff = str2double(get(handles.noiseCutoff_editText, 'String'));

tag = get(handles.baseline_buttongroup.SelectedObject , 'Tag');
handles.baselineType = baselineTypeForRadioButtonTag(tag);

handles.fadeType = 'lpc-hamming';
handles.pitchEval = 'yinbird';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes birdsong_model0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = birdsong_model_0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function filename_editText_Callback(hObject, eventdata, handles)
% hObject    handle to filename_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_editText as text
%        str2double(get(hObject,'String')) returns contents of filename_editText as a double
filename = get(hObject, 'String');
if(~strcmp(filename , handles.filename))
    load_file(hObject , handles , filename);
end
handles = guidata(hObject);
axes1_plot(handles);
axes4_plot(handles);

function load_file(hObject , handles , filename)
[xx,fs]=audioread(filename);
handles.fs = fs;
[rows,cols] = size(xx);
if (cols>rows)
    xx = xx';
    tmp = cols;
    cols = rows;
    rows = tmp;
end
if cols==2 , xx = 0.5*(xx(:,1) + xx(:,2)); end %convert stereo to mono
handles.xx = xx;
handles.filename = filename;    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function filename_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadFile_pushbutton.
function loadFile_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFile_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,inputPathName] = uigetfile('../sounds/*.wav','Select an input wav file');

if (filename ~= 0)
    set(handles.filename_editText, 'String',[inputPathName , filename]);
    load_file(hObject , handles , [inputPathName , filename]);
    handles = guidata(hObject);
    axes1_plot(handles);
    axes4_plot(handles);
end




% --- Executes on button press in synthesize_button.
function synthesize_button_Callback(hObject, eventdata, handles)
% hObject    handle to synthesize_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[y , freq , base , amp , freq_time] = chirpeval_m0(handles.xx , handles.fs , handles.winDur ,  ...
            handles.stepSize, handles.noiseCutoff , 1800 , handles.fadeType , handles.baselineType , ...
            handles.pitchEval);
handles.yy = y;
handles.baseline = base;
guidata(hObject, handles);


axes23_plot(handles);
play_synth_sound(handles);

timescale = 1;
if(length(handles.xx) <= handles.fs)
    %less than 1 sec -> scale in miliseconds
    timescale = 1000;
end

win_samples = floor(handles.winDur*handles.fs);
step_samples = floor(handles.stepSize * handles.fs);
overlap_samples = win_samples - step_samples;

freq(freq==0) = NaN;
axes(handles.axes4);
spectrogram(handles.xx,win_samples,overlap_samples,[],handles.fs,'yaxis');
hold on
plot(freq_time*timescale,freq/1000,'b');
hold off
axis([0 inf 0 10]);



function axes1_plot(handles)

time = (0:length(handles.xx)-1)/handles.fs;
axes(handles.axes1);
plot(time , handles.xx , handles.xx_color);
xlim([0,inf]);

function axes23_plot(handles)
time = (0:length(handles.xx)-1)'/handles.fs;

yy_plot = zeros(size(time));

if (get(handles.chirpPlot_checkbox , 'Value') == 1.0)
    yy_plot = handles.yy; 
end

if (get(handles.baselinePlot_checkbox , 'Value') == 1.0)
    yy_plot = yy_plot + handles.baseline;
end

axes(handles.axes2);
plot(time , yy_plot , handles.yy_color);
xlim([0,inf]);

axes(handles.axes3);
plot(time , handles.xx , handles.xx_color);
hold on;
plot(time , yy_plot , handles.yy_color);
hold off;
xlim([0,inf]);

function axes4_plot(handles)

win_samples = floor(handles.winDur*handles.fs);
step_samples = floor(handles.stepSize * handles.fs);
overlap_samples = win_samples - step_samples;
axes(handles.axes4);
spectrogram(handles.xx,win_samples,overlap_samples,[],handles.fs,'yaxis');
axis([0 inf 0 10]);


function noiseCutoff_editText_Callback(hObject, eventdata, handles)
% hObject    handle to noiseCutoff_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noiseCutoff_editText as text
%        str2double(get(hObject,'String')) returns contents of noiseCutoff_editText as a double

s = get(hObject , 'String');
v = str2double(s);
handles.noiseCutoff = v;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function noiseCutoff_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseCutoff_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winDur_editText_Callback(hObject, eventdata, handles)
% hObject    handle to winDur_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winDur_editText as text
%        str2double(get(hObject,'String')) returns contents of winDur_editText as a double
s = get(hObject , 'String');
v = str2double(s);
handles.winDur = v;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function winDur_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winDur_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepSize_editText_Callback(hObject, eventdata, handles)
% hObject    handle to stepSize_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSize_editText as text
%        str2double(get(hObject,'String')) returns contents of stepSize_editText as a double
s = get(hObject , 'String');
v = str2double(s);
handles.stepSize = v;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stepSize_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSize_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in baseline_checkbox.
function baseline_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to baseline_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of baseline_checkbox

v = get(hObject , 'Value');
radioGroup = handles.baseline_buttongroup;


if (v == 1.0)
    set(radioGroup, 'Visible' , 'on');
    tag = get(radioGroup.SelectedObject , 'Tag');
    handles.baselineType = baselineTypeForRadioButtonTag(tag);
    
elseif (v == 0)
    set(radioGroup, 'Visible' , 'off'); 
    handles.baselineType = 'non';
end
guidata(hObject, handles);

% --- Executes when selected object is changed in baseline_buttongroup.
function baseline_buttongroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in baseline_buttongroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag = get(hObject , 'Tag');
handles.baselineType = baselineTypeForRadioButtonTag(tag);
guidata(hObject, handles);

function str = baselineTypeForRadioButtonTag(tag)
    str = 'non';
    switch tag
        case 'movingAverage_radiobutton'
            str = 'moving average';
        case 'idealLowPass_radiobutton'
            str = 'ideal low pass';
        case 'ampAverage_radiobutton'
            str = 'amplitude average sampling';
    end
    

% --- Executes when selected object is changed in fade_buttongroup.
function fade_buttongroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in fade_buttongroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag = get(hObject , 'Tag');
str = '';

switch (tag)
    case 'fade1_radiobutton'
        str = 'lpc-hamming';
    case 'fade2_radiobutton'
        str = 'lpc-poly';
    case 'fade3_radiobutton'
        str = 'exponential';
    case 'fade4_radiobutton'
        str = 'hermite';
end
handles.fadeType = str;
guidata(hObject , handles);


% --- Executes on button press in segment_checkbox.
function segment_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to segment_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of segment_checkbox
isCheck = get(hObject , 'Value');
if (isCheck == 1.0)
    set(handles.segment_panel, 'Visible' , 'on');
else
    set(handles.segment_panel, 'Visible' , 'off');
    load_file(hObject , handles,handles.filename);
    handles = guidata(hObject);
    axes1_plot(handles);
    axes4_plot(handles);
end
guidata(hObject , handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over startSegment_editText.
function startSegment_editText_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to startSegment_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
[x,y] = ginput(1);
set(hObject , 'String' , num2str(x));


% --- Executes during object creation, after setting all properties.
function startSegment_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startSegment_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function endSegment_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endSegment_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over endSegment_editText.
function endSegment_editText_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to endSegment_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
[x,y] = ginput(1);
set(hObject , 'String' , num2str(x));



% --- Executes on button press in segmentation_pushbutton.
function segmentation_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startTime = str2double(get(handles.startSegment_editText , 'String'));
endTime = str2double(get(handles.endSegment_editText , 'String'));

startIndex = floor(startTime * handles.fs);
startIndex = max(startIndex , 0);
endIndex = floor(endTime * handles.fs);
endIndex = min(endIndex , length(handles.xx));

handles.xx = handles.xx(startIndex:endIndex);
guidata(hObject , handles);
axes1_plot(handles);
axes4_plot(handles);


% --- Executes on button press in play1_button.
function play1_button_Callback(hObject, eventdata, handles)
% hObject    handle to play1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

soundsc(handles.xx , handles.fs);

% --- Executes on button press in play2_pushbutton.
function play2_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to play2_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
play_synth_sound(handles);

function play_synth_sound(handles)
yy_sound = zeros(size(handles.yy));

if (get(handles.chirpPlot_checkbox , 'Value') == 1.0)
    yy_sound = handles.yy; 
end

if (get(handles.baselinePlot_checkbox , 'Value') == 1.0)
    yy_sound = yy_sound + handles.baseline;
end

soundsc(yy_sound , handles.fs);


% --- Executes on button press in chirpPlot_checkbox.
function chirpPlot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to chirpPlot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chirpPlot_checkbox
axes23_plot(handles);


% --- Executes on button press in baselinePlot_checkbox.
function baselinePlot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to baselinePlot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of baselinePlot_checkbox
axes23_plot(handles);


% --- Executes when selected object is changed in pitch_buttongroup.
function pitch_buttongroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in pitch_buttongroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag = get(hObject , 'Tag');
str = '';

switch (tag)
    case 'yinbird_radiobutton'
        str = 'yinbird';
    case 'zcr_radiobutton'
        str = 'zcr';
end
handles.pitchEval = str;
guidata(hObject , handles);
