function varargout = kitAskUser(varargin)
% KITASKUSER MATLAB code for kitAskUser.fig
%      KITASKUSER, by itself, creates a new KITASKUSER or raises the existing
%      singleton*.
%
%      H = KITASKUSER returns the handle to a new KITASKUSER or the handle to
%      the existing singleton*.
%
%      KITASKUSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KITASKUSER.M with the given input arguments.
%
%      KITASKUSER('Property','Value',...) creates a new KITASKUSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kitAskUser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kitAskUser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kitAskUser

% Last Modified by GUIDE v2.5 20-May-2015 14:11:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kitAskUser_OpeningFcn, ...
                   'gui_OutputFcn',  @kitAskUser_OutputFcn, ...
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


% --- Executes just before kitAskUser is made visible.
function kitAskUser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kitAskUser (see VARARGIN)

% Output is the job struct which is passed in as 1st parameter..
handles.output = varargin{1};

% Populate with current options.
opts = handles.output.options;

set(handles.minSpotsPerFrame,'string',num2str(opts.minSpotsPerFrame));
set(handles.betterBackground,'value',opts.betterBackground);
set(handles.robustStats,'value',opts.robustStats);
set(handles.maxCloseGap,'string',num2str(opts.maxCloseGap));
set(handles.minSearchRadiusInliers,'string',num2str(opts.minSearchRadius(1)));
set(handles.minSearchRadiusUnaligned,'string',num2str(opts.minSearchRadius(2)));
set(handles.minSearchRadiusLagging,'string',num2str(opts.minSearchRadius(3)));
set(handles.maxSearchRadiusInliers,'string',num2str(opts.maxSearchRadius(1)));
set(handles.maxSearchRadiusUnaligned,'string',num2str(opts.maxSearchRadius(2)));
set(handles.maxSearchRadiusLagging,'string',num2str(opts.maxSearchRadius(3)));
set(handles.coordSysChannel,'string',num2str(opts.coordSystemChannel));
set(handles.coordSys,'value',coordSysToIndex(opts.coordSystem));
set(handles.maxAlignAngle,'string',num2str(opts.maxSisterAlignmentAngle));
set(handles.maxSisterSep,'string',num2str(opts.maxSisterSeparation));
set(handles.minTrackOverlap,'string',num2str(opts.minSisterTrackOverlap));
set(handles.minSpotsPerFrame','string',num2str(opts.minSpotsPerFrame));
set(handles.trackMethod1,'value',coordModeToIndex(opts.coordMode{1}));
set(handles.trackMethod2,'value',coordModeToIndex(opts.coordMode{2}));
set(handles.trackMethod3,'value',coordModeToIndex(opts.coordMode{3}));
set(handles.spotMethod1,'value',spotModeToIndex(opts.spotMode{1}));
set(handles.spotMethod2,'value',spotModeToIndex(opts.spotMode{2}));
set(handles.spotMethod3,'value',spotModeToIndex(opts.spotMode{3}));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kitAskUser wait for user response (see UIRESUME)
uiwait(handles.figure1);

function idx = coordModeToIndex(coordMode)
switch coordMode
    case 'centroid'
        idx = 1;
    case 'gaussian'
        idx = 2;
    case 'none'
        idx = 3;
end

function idx = coordSysToIndex(coordSys)
switch coordSys
    case 'plate'
        idx = 1;
    case 'poles'
        idx = 2;
    case 'image'
        idx = 3;
end

function idx = spotModeToIndex(spotMode)
switch spotMode
    case 'histcut'
        idx = 1;
    case 'wavelet'
        idx = 2;
    case 'none'
        idx = 3;
end


% --- Outputs from this function are returned to the command line.
function varargout = kitAskUser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update options.
opts = handles.output.options;

opts.minSpotsPerFrame = str2double(get(handles.minSpotsPerFrame,'string'));
opts.betterBackground = get(handles.betterBackground,'value');
opts.robustStats = get(handles.robustStats,'value');
opts.maxCloseGap = str2double(get(handles.maxCloseGap,'string'));
opts.minSearchRadius(1) = str2double(get(handles.minSearchRadiusInliers,'string'));
opts.minSearchRadius(2) = str2double(get(handles.minSearchRadiusUnaligned,'string'));
opts.minSearchRadius(3) = str2double(get(handles.minSearchRadiusLagging,'string'));
opts.maxSearchRadius(1) = str2double(get(handles.maxSearchRadiusInliers,'string'));
opts.maxSearchRadius(2) = str2double(get(handles.maxSearchRadiusUnaligned,'string'));
opts.maxSearchRadius(3) = str2double(get(handles.maxSearchRadiusLagging,'string'));
opts.coordSystemChannel = str2double(get(handles.coordSysChannel,'string'));
coordSysModes = get(handles.coordSys,'string');
opts.coordSystem = lower(coordSysModes{get(handles.coordSys,'value')});
opts.maxSisterAlignmentAngle = str2double(get(handles.maxAlignAngle,'string'));
opts.maxSisterSeparation = str2double(get(handles.maxSisterSep,'string'));
opts.minSisterTrackOverlap = str2double(get(handles.minTrackOverlap,'string'));
opts.minSpotsPerFrame = str2double(get(handles.minSpotsPerFrame','string'));
trackingModes = get(handles.trackMethod1,'string');
opts.coordMode{1} = lower(trackingModes{get(handles.trackMethod1,'value')});
opts.coordMode{2} = lower(trackingModes{get(handles.trackMethod2,'value')});
opts.coordMode{3} = lower(trackingModes{get(handles.trackMethod3,'value')});
spotModes = get(handles.spotMethod1,'string');
opts.spotMode{1} = lower(spotModes{get(handles.spotMethod1,'value')});
opts.spotMode{2} = lower(spotModes{get(handles.spotMethod2,'value')});
opts.spotMode{3} = lower(spotModes{get(handles.spotMethod3,'value')});

handles.output.options = opts;
handles.output.filename = get(handles.jobsetName,'string');
varargout{1} = handles.output;
% Close GUI.
close(hObject);



function movieDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to movieDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movieDirectory as text
%        str2double(get(hObject,'String')) returns contents of movieDirectory as a double


% --- Executes during object creation, after setting all properties.
function movieDirectory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movieDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in movies.
function movies_Callback(hObject, eventdata, handles)
% hObject    handle to movies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns movies contents as cell array
%        contents{get(hObject,'Value')} returns selected item from movies


% --- Executes during object creation, after setting all properties.
function movies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Min',0.0,'Max',2.0);


% --- Executes on button press in trackNowBtn.
function trackNowBtn_Callback(hObject, eventdata, handles)
% hObject    handle to trackNowBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update job struct.
handles.output.movieDirectory = get(handles.movieDirectory, 'String');
% Selected files.
movieFiles = get(handles.movies, 'String');
selFiles = get(handles.movies, 'Value');
handles.output.movieFiles = movieFiles(selFiles);

% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in selectDirBtn.
function selectDirBtn_Callback(hObject, eventdata, handles)
% hObject    handle to selectDirBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirName = uigetdir([], 'Select directory tree containing movies');
if ~isempty(dirName)
    set(handles.movieDirectory, 'String', dirName);
    % Find movie files.
    movieFiles = kitFindFiles(dirName, kitSupportedFormats());
    % Strip search directory from filenames.
    for i=1:length(movieFiles)
      movieFiles{i} = strrep(movieFiles{i},[dirName filesep],'');
    end
    set(handles.movies, 'String', movieFiles);
end

% --- Executes on selection change in trackMethod1.
function trackMethod1_Callback(hObject, eventdata, handles)
% hObject    handle to trackMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trackMethod1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trackMethod1


% --- Executes during object creation, after setting all properties.
function trackMethod1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in linkingMethod.
function linkingMethod_Callback(hObject, eventdata, handles)
% hObject    handle to linkingMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns linkingMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from linkingMethod


% --- Executes during object creation, after setting all properties.
function linkingMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linkingMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ktChannel_Callback(hObject, eventdata, handles)
% hObject    handle to ktChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ktChannel as text
%        str2double(get(hObject,'String')) returns contents of ktChannel as a double


% --- Executes during object creation, after setting all properties.
function ktChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ktChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in robustStats.
function robustStats_Callback(hObject, eventdata, handles)
% hObject    handle to robustStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of robustStats



function maxCloseGap_Callback(hObject, eventdata, handles)
% hObject    handle to maxCloseGap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxCloseGap as text
%        str2double(get(hObject,'String')) returns contents of maxCloseGap as a double


% --- Executes during object creation, after setting all properties.
function maxCloseGap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxCloseGap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minSearchRadiusInliers_Callback(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusInliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSearchRadiusInliers as text
%        str2double(get(hObject,'String')) returns contents of minSearchRadiusInliers as a double


% --- Executes during object creation, after setting all properties.
function minSearchRadiusInliers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusInliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSearchRadiusInliers_Callback(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusInliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSearchRadiusInliers as text
%        str2double(get(hObject,'String')) returns contents of maxSearchRadiusInliers as a double


% --- Executes during object creation, after setting all properties.
function maxSearchRadiusInliers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusInliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coordSysChannel_Callback(hObject, eventdata, handles)
% hObject    handle to coordSysChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coordSysChannel as text
%        str2double(get(hObject,'String')) returns contents of coordSysChannel as a double


% --- Executes during object creation, after setting all properties.
function coordSysChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coordSysChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in coordSys.
function coordSys_Callback(hObject, eventdata, handles)
% hObject    handle to coordSys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coordSys contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coordSys


% --- Executes during object creation, after setting all properties.
function coordSys_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coordSys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxCloseGapLagging_Callback(hObject, eventdata, handles)
% hObject    handle to maxCloseGapLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxCloseGapLagging as text
%        str2double(get(hObject,'String')) returns contents of maxCloseGapLagging as a double


% --- Executes during object creation, after setting all properties.
function maxCloseGapLagging_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxCloseGapLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minSearchRadiusLagging_Callback(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSearchRadiusLagging as text
%        str2double(get(hObject,'String')) returns contents of minSearchRadiusLagging as a double


% --- Executes during object creation, after setting all properties.
function minSearchRadiusLagging_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSearchRadiusLagging_Callback(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSearchRadiusLagging as text
%        str2double(get(hObject,'String')) returns contents of maxSearchRadiusLagging as a double


% --- Executes during object creation, after setting all properties.
function maxSearchRadiusLagging_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusLagging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSearchRadiusUnaligned_Callback(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSearchRadiusUnaligned as text
%        str2double(get(hObject,'String')) returns contents of maxSearchRadiusUnaligned as a double


% --- Executes during object creation, after setting all properties.
function maxSearchRadiusUnaligned_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSearchRadiusUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minSearchRadiusUnaligned_Callback(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSearchRadiusUnaligned as text
%        str2double(get(hObject,'String')) returns contents of minSearchRadiusUnaligned as a double


% --- Executes during object creation, after setting all properties.
function minSearchRadiusUnaligned_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSearchRadiusUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxCloseGapUnaligned_Callback(hObject, eventdata, handles)
% hObject    handle to maxCloseGapUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxCloseGapUnaligned as text
%        str2double(get(hObject,'String')) returns contents of maxCloseGapUnaligned as a double


% --- Executes during object creation, after setting all properties.
function maxCloseGapUnaligned_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxCloseGapUnaligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxAlignAngle_Callback(hObject, eventdata, handles)
% hObject    handle to maxAlignAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxAlignAngle as text
%        str2double(get(hObject,'String')) returns contents of maxAlignAngle as a double


% --- Executes during object creation, after setting all properties.
function maxAlignAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxAlignAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxSisterSep_Callback(hObject, eventdata, handles)
% hObject    handle to maxSisterSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxSisterSep as text
%        str2double(get(hObject,'String')) returns contents of maxSisterSep as a double


% --- Executes during object creation, after setting all properties.
function maxSisterSep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSisterSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in betterBackground.
function betterBackground_Callback(hObject, eventdata, handles)
% hObject    handle to betterBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of betterBackground



function minTrackOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to minTrackOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minTrackOverlap as text
%        str2double(get(hObject,'String')) returns contents of minTrackOverlap as a double


% --- Executes during object creation, after setting all properties.
function minTrackOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minTrackOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minSpotsPerFrame_Callback(hObject, eventdata, handles)
% hObject    handle to minSpotsPerFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSpotsPerFrame as text
%        str2double(get(hObject,'String')) returns contents of minSpotsPerFrame as a double


% --- Executes during object creation, after setting all properties.
function minSpotsPerFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSpotsPerFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trackMethod3.
function trackMethod3_Callback(hObject, eventdata, handles)
% hObject    handle to trackMethod3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trackMethod3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trackMethod3


% --- Executes during object creation, after setting all properties.
function trackMethod3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackMethod3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trackMethod2.
function trackMethod2_Callback(hObject, eventdata, handles)
% hObject    handle to trackMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trackMethod2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trackMethod2


% --- Executes during object creation, after setting all properties.
function trackMethod2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jobsetName_Callback(hObject, eventdata, handles)
% hObject    handle to jobsetName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jobsetName as text
%        str2double(get(hObject,'String')) returns contents of jobsetName as a double


% --- Executes during object creation, after setting all properties.
function jobsetName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jobsetName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spotMethod3.
function spotMethod3_Callback(hObject, eventdata, handles)
% hObject    handle to spotMethod3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spotMethod3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spotMethod3


% --- Executes during object creation, after setting all properties.
function spotMethod3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotMethod3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spotMethod2.
function spotMethod2_Callback(hObject, eventdata, handles)
% hObject    handle to spotMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spotMethod2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spotMethod2


% --- Executes during object creation, after setting all properties.
function spotMethod2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spotMethod1.
function spotMethod1_Callback(hObject, eventdata, handles)
% hObject    handle to spotMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spotMethod1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spotMethod1


% --- Executes during object creation, after setting all properties.
function spotMethod1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
