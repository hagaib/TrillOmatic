function varargout = spectrogram_options(varargin)
% SPECTROGRAM_OPTIONS MATLAB code for spectrogram_options.fig
%      SPECTROGRAM_OPTIONS, by itself, creates a new SPECTROGRAM_OPTIONS or raises the existing
%      singleton*.
%
%      H = SPECTROGRAM_OPTIONS returns the handle to a new SPECTROGRAM_OPTIONS or the handle to
%      the existing singleton*.
%
%      SPECTROGRAM_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTROGRAM_OPTIONS.M with the given input arguments.
%
%      SPECTROGRAM_OPTIONS('Property','Value',...) creates a new SPECTROGRAM_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spectrogram_options_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spectrogram_options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spectrogram_options

% Last Modified by GUIDE v2.5 09-Jan-2018 17:36:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spectrogram_options_OpeningFcn, ...
                   'gui_OutputFcn',  @spectrogram_options_OutputFcn, ...
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


% --- Executes just before spectrogram_options is made visible.
function spectrogram_options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spectrogram_options (see VARARGIN)

% Choose default command line output for spectrogram_options
handles.output = hObject;


if(~isempty(varargin))
    spect_params = varargin{1};
    handles.edit_winsize.String = spect_params.winsize;
    handles.popupmenu_wintype.Value =  spect_params.wintype;
    handles.edit_overlap.String = spect_params.overlap;
    handles.edit_bw_bottom.String = num2str(spect_params.bw(1));
    handles.edit_bw_top.String = num2str(spect_params.bw(2));
    if(~isempty(spect_params.nfft) )
        handles.checkbox_nfft.Value = 1;
        handles.edit_nfft.String = spect_params.nfft;
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spectrogram_options wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isequal(get(hObject , 'waitstatus') , 'waiting') )
    uiresume(hObject);
else
    delete(hObject);
end

% --- Outputs from this function are returned to the command line.
function varargout = spectrogram_options_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( ~isempty(handles) )
    if ( isequal(handles.output , 'OK') )
        spect_params.winsize  = str2double(handles.edit_winsize.String);
        spect_params.wintype = handles.popupmenu_wintype.Value;
        spect_params.overlap = str2double(handles.edit_overlap.String);
        spect_params.nfft = [];
        
        %bandwidth
        bottom = str2double(handles.edit_bw_bottom.String);
        top = str2double(handles.edit_bw_top.String);
        spect_params.bw = [bottom , top];
        
        if (handles.checkbox_nfft.Value)
            spect_params.nfft = str2double(handles.edit_nfft.String);
        end
        varargout{1} = spect_params;
    else
        varargout{1} = [];
    end
    
    delete(handles.figure1)
else
    varargout{1} = [];
end

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('OK')
handles.output = 'OK';
guidata(hObject, handles);
uiresume(handles.figure1)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('cancel')
handles.output = 'cancel';
guidata(hObject, handles);
uiresume(hObject.Parent)



function edit_winsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_winsize as text
%        str2double(get(hObject,'String')) returns contents of edit_winsize as a double


% --- Executes during object creation, after setting all properties.
function edit_winsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_winsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_overlap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_overlap as text
%        str2double(get(hObject,'String')) returns contents of edit_overlap as a double


% --- Executes during object creation, after setting all properties.
function edit_overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nfft_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nfft as text
%        str2double(get(hObject,'String')) returns contents of edit_nfft as a double


% --- Executes during object creation, after setting all properties.
function edit_nfft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_nfft.
function checkbox_nfft_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nfft


% --- Executes on selection change in popupmenu_wintype.
function popupmenu_wintype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_wintype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_wintype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_wintype


% --- Executes during object creation, after setting all properties.
function popupmenu_wintype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_wintype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

m = matfile('wintypes.mat');
set(hObject , 'String' , m.names);



function edit_bw_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bw_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bw_bottom as text
%        str2double(get(hObject,'String')) returns contents of edit_bw_bottom as a double


% --- Executes during object creation, after setting all properties.
function edit_bw_bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bw_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bw_top_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bw_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bw_top as text
%        str2double(get(hObject,'String')) returns contents of edit_bw_top as a double


% --- Executes during object creation, after setting all properties.
function edit_bw_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bw_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
