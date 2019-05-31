function varargout = PreprocessGUI(varargin)
% PREPROCESSGUI MATLAB code for PreprocessGUI.fig
%      PREPROCESSGUI, by itself, creates a new PREPROCESSGUI or raises the existing
%      singleton*.
%
%      H = PREPROCESSGUI returns the handle to a new PREPROCESSGUI or the handle to
%      the existing singleton*.
%
%      PREPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSGUI.M with the given input arguments.
%
%      PREPROCESSGUI('Property','Value',...) creates a new PREPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PreprocessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PreprocessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PreprocessGUI

% Last Modified by GUIDE v2.5 14-May-2019 14:45:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreprocessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PreprocessGUI_OutputFcn, ...
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


% --- Executes just before PreprocessGUI is made visible.
function PreprocessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreprocessGUI (see VARARGIN)

% Choose default command line output for PreprocessGUI
handles.output = hObject;
warning('off', 'MATLAB:unknownObjectNowStruct');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PreprocessGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PreprocessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1. select folder
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folder = uigetdir;
files=dir(fullfile(handles.folder, 'AD0_*.mat'));
for i = 1:length(files)
    figures{i} = files(i).name;
end
set(handles.listbox1, 'String', figures);
guidata(hObject,handles);



% --- Executes on button press in pushbutton2. accept
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listacc = cellstr(get(handles.listbox1, 'String'));
selectedacc = get(handles.listbox1, 'Value');
choice_listbox1 = listacc(selectedacc);
updated_listbox2 = cellstr(get(handles.listbox2, 'String'));
newmenuacc = [choice_listbox1; updated_listbox2];
set(handles.listbox2, 'String', newmenuacc);
current2 = get(handles.listbox1, 'String');
newitems2 = current2;
newitems2(selectedacc) = [];
set(handles.listbox1, 'String', newitems2);



% --- Executes on button press in pushbutton3. reject
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listrej = cellstr(get(handles.listbox1, 'String'));
selectedrej = get(handles.listbox1, 'Value');
choice_listbox4r = listrej(selectedrej);
updated_listbox3 = cellstr(get(handles.listbox3, 'String'));
newmenurej = [choice_listbox4r; updated_listbox3];
set(handles.listbox3, 'String', newmenurej);
current1 = get(handles.listbox1, 'String');
newitems1 = current1;
newitems1(selectedrej) = [];
set(handles.listbox1, 'String', newitems1);


% --- Executes on selection change in listbox1. all traces
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

list = get(handles.listbox1, 'string');
selected = get(handles.listbox1, 'value');
load(fullfile(handles.folder, cell2mat(list(selected))));
name = list{selected};
I = eval(name(1:end-4));
axes(handles.axes1)
cla;
plot(I.rawdata, 'Color', [0.3010, 0.7450, 0.9330])
title(strcat(name(1:end-4), ' Full Trace'), 'Interpreter', 'none');
xlabel('Frames (1/10000 s)');
ylabel('mV');
ylim auto
xlim([0 length(I.rawdata)]);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.accepted
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3. rejected
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5. acc back
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listaccr = cellstr(get(handles.listbox2, 'String'));
selectedaccr = get(handles.listbox2, 'Value');
choice_listbox2r = listaccr(selectedaccr);
updated_listbox1 = cellstr(get(handles.listbox1, 'String'));
newmenuaccr = [choice_listbox2r; updated_listbox1];
set(handles.listbox1, 'String', newmenuaccr);
current2 = get(handles.listbox2, 'String');
newitems2 = current2;
newitems2(selectedaccr) = [];
set(handles.listbox2, 'String', newitems2);


% --- Executes on button press in pushbutton6. rejback
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listrej = cellstr(get(handles.listbox1, 'String'));
selectedrej = get(handles.listbox1, 'Value');
choice_listbox1r = listrej(selectedrej);
updated_listbox3 = cellstr(get(handles.listbox3, 'String'));
newmenurej = [choice_listbox1r; updated_listbox3];
set(handles.listbox3, 'String', newmenurej);
current1 = get(handles.listbox1, 'String');
newitems1 = current1;
newitems1(selectedrej) = [];
set(handles.listbox1, 'String', newitems1);


% --- Executes on button press in pushbutton7. save
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(get(handles.listbox1, 'String')) == 0
    f = warndlg('You have not processed all traces','Warning');
else
    opts.Interpreter = "none";
    opts.Default = "No";
    answer = questdlg('You are about to move your rejected traces -- Continue?', 'Move Traces',  'Yes', 'No',  opts); 
end 

if strcmp(answer, 'Yes')
    allacc = get(handles.listbox3, 'String');
    mkdir(fullfile(handles.folder, 'Rejected Traces'));
    for i = 1:length(allacc)-1
        filestring = string(allacc(i));
        movefile(fullfile(handles.folder, filestring), fullfile(handles.folder, 'Rejected Traces'));
    end
else
    return;
end
warning('on', 'MATLAB:unknownObjectNowStruct');
warning('on', 'MATLAB:MKDIR:DirectoryExists');


