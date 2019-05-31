function varargout = Accept_Reject_GUI(varargin)

%Written by CRW, 4 Oct 2018
    %Last Updated 10 Oct 2018
    
% ACCEPT_REJECT_GUI MATLAB code for Accept_Reject_GUI.fig
%      ACCEPT_REJECT_GUI, by itself, creates a new ACCEPT_REJECT_GUI or raises the existing
%      singleton*.
%
%      H = ACCEPT_REJECT_GUI returns the handle to a new ACCEPT_REJECT_GUI or the handle to
%      the existing singleton*.
%
%      ACCEPT_REJECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACCEPT_REJECT_GUI.M with the given input arguments.
%
%      ACCEPT_REJECT_GUI('Property','Value',...) creates a new ACCEPT_REJECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Accept_Reject_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Accept_Reject_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Accept_Reject_GUI

% Last Modified by GUIDE v2.5 10-Oct-2018 13:52:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Accept_Reject_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Accept_Reject_GUI_OutputFcn, ...
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


% --- Executes just before Accept_Reject_GUI is made visible.
function Accept_Reject_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Accept_Reject_GUI (see VARARGIN)

% Choose default command line output for Accept_Reject_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Accept_Reject_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Accept_Reject_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1. Accept
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listacc = cellstr(get(handles.listbox4, 'String'));
selectedacc = get(handles.listbox4, 'Value');
choice_listbox4 = listacc(selectedacc);
updated_listbox2 = cellstr(get(handles.listbox2, 'String'));
newmenuacc = [choice_listbox4; updated_listbox2];
set(handles.listbox2, 'String', newmenuacc);
current2 = get(handles.listbox4, 'String');
newitems2 = current2;
newitems2(selectedacc) = [];
set(handles.listbox4, 'String', newitems2);





% --- Executes on button press in pushbutton2. Reject
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listrej = cellstr(get(handles.listbox4, 'String'));
selectedrej = get(handles.listbox4, 'Value');
choice_listbox4r = listrej(selectedrej);
updated_listbox3 = cellstr(get(handles.listbox3, 'String'));
newmenurej = [choice_listbox4r; updated_listbox3];
set(handles.listbox3, 'String', newmenurej);
current1 = get(handles.listbox4, 'String');
newitems1 = current1;
newitems1(selectedrej) = [];
set(handles.listbox4, 'String', newitems1);



% --- Executes on button press in pushbutton6. Select Folder
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folder = uigetdir;
files=dir(fullfile(handles.folder, '*-proc.mat'));
for i = 1:length(files)
    figures{i} = files(i).name;
end
set(handles.listbox4, 'String', figures);
guidata(hObject,handles);


% --- Executes on selection change in listbox2. Accepted
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


% --- Executes on selection change in listbox3. Rejected
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


% --- Executes on selection change in listbox4. All traces
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listbox4, 'string');
selected = get(handles.listbox4, 'value');
load(fullfile(handles.folder, cell2mat(list(selected))));
name = list{selected};
I = eval(name(1:end-9));
axes(handles.axes1)
cla;
    if isfield(I, 'besseled') == 1;
        plot(I.besseled, 'Color', [0.3010, 0.7450, 0.9330]);
        hold on;
        scatter(I.event_times, repmat(mean(I.besseled),1,length(I.event_times)), 50, [1, 0.33, 0.64]); %plots events ontop of unfiltered trace
        title(strcat('Unfiltered  ', I.params.traces_file(1:end-4), ' Data with Events'), 'Interpreter', 'none');
        xlabel('Frames (1/10000 s)');
        ylabel('mV');
        ylim auto
        xlim([0 length(I.besseled)]);
    else 
        plot(I.raw_trace, 'Color', [0.3010, 0.7450, 0.9330]); 
        hold on;
        scatter(I.event_times, repmat(mean(I.raw_trace),1,length(I.event_times)), 50, [1, 0.33, 0.64]); %plots events ontop of unfiltered trace
        title(strcat('Unfiltered ', ' ', I.params.traces_file(1:end-4), ' ', ' Data with Events'), 'Interpreter', 'none');
        ylabel('mV');
        ylim auto
        xlim([0 length(I.raw_trace)]);
    end
    
axes(handles.axes4);
cla;
    plot(zscore(I.filtered_trace), 'Color', [0.3010, 0.7450, 0.9330]); %plots filtered trace
    hold on;
    scatter(I.event_times, repmat(mean(zscore(I.filtered_trace)), 1, length(I.event_times)), 40, [1, 0.33, 0.64], 'filled');
    hold on;
    hline = refline([0 I.params.init_method.threshold]);
    hline.Color = [0.88 0.13 0.54];
    hline.LineWidth = 1;
    hline.DisplayName = 'Threshold for Spike Detection (2.35 Std)';
    ylabel('ZScore of Filtered Trace');
    ylim auto
    xlim([0 length(I.filtered_trace)]);
    
    %scatter times vs amps 84%	9%	41%
axes(handles.axes5);
cla;
    scatter(I.event_times, I.event_amp, 40, [.84 .09 .41], 'filled');
    xlabel('Frames (1/10000 s)');
    ylabel('Amplitude of Event (mV)');
    ylim auto
    xlim([0 length(I.filtered_trace)]);



% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listaccr = cellstr(get(handles.listbox2, 'String'));
selectedaccr = get(handles.listbox2, 'Value');
choice_listbox2r = listaccr(selectedaccr);
updated_listbox41 = cellstr(get(handles.listbox4, 'String'));
newmenuaccr = [choice_listbox2r; updated_listbox41];
set(handles.listbox4, 'String', newmenuaccr);
current2 = get(handles.listbox2, 'String');
newitems2 = current2;
newitems2(selectedaccr) = [];
set(handles.listbox2, 'String', newitems2);



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listrejr = cellstr(get(handles.listbox3, 'String'));
selectedrejr = get(handles.listbox3, 'Value');
choice_listbox3r = listrejr(selectedrejr);
updated_listbox4 = cellstr(get(handles.listbox4, 'String'));
newmenurejr = [choice_listbox3r; updated_listbox4];
set(handles.listbox4, 'String', newmenurejr);
current3 = get(handles.listbox3, 'String');
newitems3 = current3;
newitems3(selectedrejr) = [];
set(handles.listbox3, 'String', newitems3);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(get(handles.listbox1, 'String')) == 0
    f = warndlg('You have not processed all traces','Warning');
end 

allacc = get(handles.listbox3, 'String');
selected2 = get(handles.listbox3, 'value');
load(fullfile(handles.folder, cell2mat(allacc(selected2))));
name2 = allacc{selected2};
I2 = eval(name2(1:end-9));
epoch = I2(1).params.epoch;
cell = I2(1).params.cell;
T = array2table(allacc);
writetable(T,fullfile(handles.folder, strcat('Accepted_Traces_cell', cell, '_epoch', epoch, '.xlsx')));
