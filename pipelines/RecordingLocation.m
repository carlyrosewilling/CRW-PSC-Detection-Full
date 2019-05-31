function varargout = RecordingLocation(varargin)
% RECORDINGLOCATION MATLAB code for RecordingLocation.fig
%      RECORDINGLOCATION, by itself, creates a new RECORDINGLOCATION or raises the existing
%      singleton*.
%
%      H = ATLASREDO returns the handle to a new RECORDINGLOCATION or the handle to
%      the existing singleton*.
%
%      RECORDINGLOCATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECORDINGLOCATIONO.M with the given input arguments.
%
%      RECORDINGLOCATION('Property','Value',...) creates a new RECORDINGLOCATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RecordingLocation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RecordingLocation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help atlasredo

% Last Modified by GUIDE v2.5 04-Jan-2019 12:55:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RecordingLocation_OpeningFcn, ...
                   'gui_OutputFcn',  @RecordingLocation_OutputFcn, ...
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


% --- Executes just before atlasredo is made visible.
function RecordingLocation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to atlasredo (see VARARGIN)

% Choose default command line output for atlasredo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes atlasredo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RecordingLocation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectorigin.
function selectorigin_Callback(hObject, eventdata, handles)
% hObject    handle to selectorigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[xcoordinate, ycoordinate] = ginput;
coordinate = [xcoordinate, ycoordinate];
set(handles.selectorigin, 'UserData', coordinate);
guidata(hObject,handles)


function x_Callback(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double



% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_Callback(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y as text
%        str2double(get(hObject,'String')) returns contents of y as a double


% --- Executes during object creation, after setting all properties.
function y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in slicemenu.
function slicemenu_Callback(hObject, eventdata, handles)
% hObject    handle to slicemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns slicemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from slicemenu
files=dir(fullfile('//Volumes/Carly Rose/2 - Code/1 - MATLAB/CRW-PSC-Detection-master/Atlas Tiffs', '*.tiff'));
for i = 1:length(files)
    slices{i} = files(i).name;
end
handles.slicemenu.String = slices;
guidata(hObject, handles);
contents = cellstr(get(hObject, 'String'));
selected = contents{get(hObject, 'Value')};
SelectedSlice = imread(fullfile('//Volumes/Carly Rose/2 - Code/1 - MATLAB/CRW-PSC-Detection-master/Atlas Tiffs', contents{get(hObject, 'Value')}));
SelectedSlice = SelectedSlice(:,:,1:3);
axes(handles.axes1);
hold off;
imshow(SelectedSlice, 'InitialMagnification', 'fit');
set(handles.slicemenu, 'UserData', SelectedSlice)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function slicemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveimage.
function saveimage_Callback(hObject, eventdata, handles)
% hObject    handle to saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folder = uigetdir;
fullcoordinates = get(handles.drawlocation, 'UserData');
axes(handles.axes1)
hold off;
Selectedslice = get(handles.slicemenu, 'UserData');
Slicewithmarker2 = insertShape(Selectedslice, 'Circle', [fullcoordinates(1) fullcoordinates(2) 50], 'Color', [139 10 80] , 'LineWidth', 20);
imshow(Slicewithmarker2, 'InitialMagnification', 'fit');
cell = get(handles.edit3, 'String');
imwrite(Slicewithmarker2, fullfile(handles.folder,strcat('cell_', cell, '_atlas.png')));
%save(fullfile(handles.folder, strcat('cell_', cell, '_locationcoordinates.mat')), fullcoordinates);


% --- Executes on button press in drawlocation.
function drawlocation_Callback(hObject, eventdata, handles)
% hObject    handle to drawlocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
yinput = get(handles.y,'String');
coordinates = get(handles.selectorigin, 'UserData');
xcoor = coordinates(1);
ycoor = coordinates(2);
ylocation = ycoor - (str2num(yinput)/1000)*775;
xinput = get(handles.x,'String');
xlocation = xcoor + (str2num(xinput)/1000)*775;
fullcoordinates = [xlocation ylocation 50];
axes(handles.axes1)
hold off;
Selectedslice = get(handles.slicemenu, 'UserData');
Slicewithmarker = insertShape(Selectedslice, 'Circle', [xlocation ylocation 50], 'Color', [139 10 80] , 'LineWidth', 20);
imshow(Slicewithmarker, 'InitialMagnification', 'fit');
set(handles.drawlocation, 'UserData', fullcoordinates);
guidata(hObject,handles);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
