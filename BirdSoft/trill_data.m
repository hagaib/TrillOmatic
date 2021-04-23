function varargout = trill_data(varargin)
% TRILL_DATA MATLAB code for trill_data.fig
%      TRILL_DATA, by itself, creates a new TRILL_DATA or raises the existing
%      singleton*.
%
%      H = TRILL_DATA returns the handle to a new TRILL_DATA or the handle to
%      the existing singleton*.
%
%      TRILL_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRILL_DATA.M with the given input arguments.
%
%      TRILL_DATA('Property','Value',...) creates a new TRILL_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trill_data_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trill_data_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trill_data

% Last Modified by GUIDE v2.5 25-Apr-2018 18:07:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trill_data_OpeningFcn, ...
                   'gui_OutputFcn',  @trill_data_OutputFcn, ...
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


% --- Executes just before trill_data is made visible.
function trill_data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trill_data (see VARARGIN)

% Choose default command line output for trill_data
handles.output = hObject;
handles.params = varargin{1};
handles.segments = varargin{2};
filename = varargin{3};


handles.edit_filename.String = strrep(filename , '.wav' , '.csv');
handles.edit_title.String = strrep(filename , '.wav' , '');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trill_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


fill_table(handles);


% --- Outputs from this function are returned to the command line.
function varargout = trill_data_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function fill_table(handles)

table1 = handles.uitable1;
table2 = handles.uitable2;
params = handles.params;
segs = handles.segments;

syl_titles = ...
{'Syllable Count' , ...
'Trill Duration' , ...
'Syllable Rate' , ...
'Syllable Duration' , ...
'Interval Duration' , ...
'PRI' , ...
'Syllable Bandwidth' , ...
'Frequency Slope (Syllable)' ,... 
'Max Amplitude' , ...
'Attack Time' , ...
'Attack Relative Time' ,... 
'Release Time' , ...
'Release Relative Time' ,... 
'Max Band' ,...
'Min Band' , ...
'f0 Start Max' , ...
'f0 End Max' , ...
'Delta Max' , ...
'Max Slope' , ...
'f0 Start Min' , ...
'f0 End Min' , ...
'Delta Min' , ...
'Min Slope'};

syl_data = {
    params.syllable_count , ...
    params.trill_dur, ...
    1/median((params.syllable_dur(1:end-1) + params.gaps_dur)) , ...
    median(params.syllable_dur) , ...
    median(params.gaps_dur) , ...
    median(params.syllable_dur(1:end-1) + params.gaps_dur) , ...
    median(params.syllable_bw) , ...
    median(params.syllable_bw ./ params.syllable_dur) , ...
    '????' , ...
    params.trill_attack_dur , ...
    params.trill_attack_dur / params.trill_dur , ...
    params.trill_release_dur , ...
    params.trill_release_dur / params.trill_dur , ...
    max(params.syllable_f0max - params.syllable_f0min) , ...
    min(params.syllable_f0max - params.syllable_f0min) , ...
    params.syllable_f0max(1) , ...
    params.syllable_f0max(end) , ...
    params.syllable_f0max(1) - params.syllable_f0max(end) , ...
    (params.syllable_f0max(1) - params.syllable_f0max(end))/params.trill_dur , ...
    params.syllable_f0min(1) , ...
    params.syllable_f0min(end) , ...
    params.syllable_f0min(1) - params.syllable_f0min(end) , ...
    (params.syllable_f0min(1) - params.syllable_f0min(end))/params.trill_dur};

table1.ColumnName = syl_titles;
table1.Data = syl_data;


data2 = [segs' , ...
         params.syllable_f0max' , ...
         params.syllable_f0min' , ...
         params.syllable_bw' , ...
         (params.syllable_bw ./ params.syllable_dur)' , ...
         params.syllable_dur' , ...
         [params.gaps_dur , nan]' , ...
         [(params.syllable_dur(1:end-1) + params.gaps_dur)' ; nan]];

table2.ColumnName = {'Start time' , 'End Time' , 'Max f0' , 'Min f0' , ...
            'Bandwidth' , 'f slope' , 'Duration' , 'Interval Duration' , 'Period Time'};
table2.Data = data2;



% NewData = mat2cell(table2.Data,ones(size(table2.Data,1),1),ones(size(table2.Data,2),1));

% write_list_csv('excelfile.xls' , NewData);


% --- Executes on button press in pushbutton_excel.
function pushbutton_excel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_excel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = handles.edit_filename.String;
filenamecsv = strrep(filename ,'.xls' , '.csv');
table1 = handles.uitable1;
table2 = handles.uitable2;

data1 = [table1.ColumnName' ; cellfun(@num2str, table1.Data, 'UniformOutput', false)];
data2 = [table2.ColumnName' ; num2cell(table2.Data)];
title = handles.edit_title.String;


write_list_csv(filenamecsv , {title}, 'w');
write_list_csv(filenamecsv , data2 , 'a');
write_list_csv(filenamecsv , {'-'} , 'a');
write_list_csv(filenamecsv , data1 , 'a');

% xlswrite(filename , table1.ColumnName , 'A2');
% xlswrite(filename , [table1.Data{:}'] , 'A3');



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_title_Callback(hObject, eventdata, handles)
% hObject    handle to edit_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_title as text
%        str2double(get(hObject,'String')) returns contents of edit_title as a double


% --- Executes during object creation, after setting all properties.
function edit_title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
