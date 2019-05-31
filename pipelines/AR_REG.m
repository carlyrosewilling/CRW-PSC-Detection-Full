function varargout = AR_REG(varargin)
% AR_REG MATLAB code for AR_REG.fig
%      AR_REG, by itself, creates a new AR_REG or raises the existing
%      singleton*.
%
%      H = AR_REG returns the handle to a new AR_REG or the handle to
%      the existing singleton*.
%
%      AR_REG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AR_REG.M with the given input arguments.
%
%      AR_REG('Property','Value',...) creates a new AR_REG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AR_REG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AR_REG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AR_REG

% Last Modified by GUIDE v2.5 08-Mar-2019 15:49:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AR_REG_OpeningFcn, ...
                   'gui_OutputFcn',  @AR_REG_OutputFcn, ...
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


% --- Executes just before AR_REG is made visible.
function AR_REG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AR_REG (see VARARGIN)

% Choose default command line output for AR_REG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AR_REG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AR_REG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select. Select
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folder = uigetdir;
files=dir(fullfile(handles.folder, '*-proc.mat'));
for i = 1:length(files)
    figures{i} = files(i).name;
end
set(handles.listbox1, 'String', figures);
guidata(hObject,handles);


% --- Executes on selection change in listbox1. Traces
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
I = eval(name(1:end-9));
axes(handles.axes1)
cla;
plot(I.peristim1, 'Color', [0.3010, 0.7450, 0.9330]);
hold on;
scatter(I.event_times1, repmat(mean(I.peristim1),1,length(I.event_times1)), 50, [1, 0.33, 0.64]);
title(strcat('Peristimulus 1 ', I.params.traces_file(1:end-4), ' Data with Events'), 'Interpreter', 'none');
xlabel('Frames (1/10000 s)');
ylabel('mV');
ylim auto
xlim([0 length(I.peristim1)]);
    
axes(handles.axes2);
cla;
plot(zscore(I.filtered_trace1), 'Color', [0.3010, 0.7450, 0.9330]);
hold on;
scatter((I.event_times1-15), repmat(mean(zscore(I.filtered_trace1)),1,length(I.event_times1)), 50, [1, 0.33, 0.64]);
title(strcat('Filtered Trace ', I.params.traces_file(1:end-4), ' Data with Events'), 'Interpreter', 'none');
xlabel('Frames (1/10000 s)');
ylabel('mV');
ylim auto
xlim([0 length(I.filtered_trace1)]);




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


% --- Executes on button press in accept. Accept
function accept_Callback(hObject, eventdata, handles)
% hObject    handle to accept (see GCBO)
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



% --- Executes on button press in reject. 
function reject_Callback(hObject, eventdata, handles)
% hObject    handle to reject (see GCBO)
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


% --- Executes on button press in sendaccback. send accepted back
function sendaccback_Callback(hObject, eventdata, handles)
% hObject    handle to sendaccback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listaccr = cellstr(get(handles.listbox2, 'String'));
selectedaccr = get(handles.listbox2, 'Value');
choice_listbox2r = listaccr(selectedaccr);
updated_listbox11 = cellstr(get(handles.listbox1, 'String'));
newmenuaccr = [choice_listbox2r; updated_listbox11];
set(handles.listbox1, 'String', newmenuaccr);
current2 = get(handles.listbox2, 'String');
newitems2 = current2;
newitems2(selectedaccr) = [];
set(handles.listbox2, 'String', newitems2);



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


% --- Executes on button press in sendrejback. sent rej back
function sendrejback_Callback(hObject, eventdata, handles)
% hObject    handle to sendrejback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listrejr = cellstr(get(handles.listbox3, 'String'));
selectedrejr = get(handles.listbox3, 'Value');
choice_listbox3r = listrejr(selectedrejr);
updated_listbox1 = cellstr(get(handles.listbox1, 'String'));
newmenurejr = [choice_listbox3r; updated_listbox1];
set(handles.listbox1, 'String', newmenurejr);
current3 = get(handles.listbox3, 'String');
newitems3 = current3;
newitems3(selectedrejr) = [];
set(handles.listbox3, 'String', newitems3);


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


% --- Executes on button press in save. save
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(get(handles.listbox1, 'String')) == 0
    f = warndlg('You have not processed all traces','Warning');
end 

allacc = get(handles.listbox2, 'String');
selected2 = get(handles.listbox2, 'value');
load(fullfile(handles.folder, cell2mat(allacc(selected2))));
name2 = allacc{selected2};
I2 = eval(name2(1:end-9));
epoch = I2(1).params.epoch;
cell = I2(1).params.cell;
T = array2table(allacc);
writetable(T,fullfile(handles.folder, strcat('Accepted_Traces_cell', cell, '_epoch', epoch, '.xlsx')));
