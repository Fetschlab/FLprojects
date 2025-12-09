function varargout = rigBgui_advsettings(varargin)
% RIGBGUI_ADVSETTINGS MATLAB code for rigBgui_advsettings.fig
%      RIGBGUI_ADVSETTINGS, by itself, creates a new RIGBGUI_ADVSETTINGS or raises the existing
%      singleton*.
%
%      H = RIGBGUI_ADVSETTINGS returns the handle to a new RIGBGUI_ADVSETTINGS or the handle to
%      the existing singleton*.
%
%      RIGBGUI_ADVSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIGBGUI_ADVSETTINGS.M with the given input arguments.
%
%      RIGBGUI_ADVSETTINGS('Property','Value',...) creates a new RIGBGUI_ADVSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rigBgui_advsettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rigBgui_advsettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rigBgui_advsettings

% Last Modified by GUIDE v2.5 25-May-2021 10:59:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rigBgui_advsettings_OpeningFcn, ...
                   'gui_OutputFcn',  @rigBgui_advsettings_OutputFcn, ...
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


% --- Executes just before rigBgui_advsettings is made visible.
function rigBgui_advsettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rigBgui_advsettings (see VARARGIN)

% Choose default command line output for rigBgui_advsettings
handles.output = hObject;

% grab data from main GUI and settingsStruct
paradigm = getappdata(0,'paradigm_value');
subject  = getappdata(0,'subject_value');

set(handles.subject,'String',['Subject: ' subject]); % fill in subject field

load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);

% identify handles (GUI fields) to update, then grab their values from the
% local settingsStruct (location stored in TooltipString)

% pass in hObject, otherwise we'll reset the values in the main GUI too!
checkboxes = findobj(hObject,'Style','checkbox');
for n = 1:length(checkboxes)
    try
        eval(['loadedValue = ' checkboxes(n).TooltipString ';']);
        checkboxes(n).Value = loadedValue;
    catch
        checkboxes(n).Value = 0;
    end
end

editboxes = findobj(hObject,'Style','edit');
for n = 1:length(editboxes)
    try
        eval(['loadedValue = ' editboxes(n).TooltipString ';']);
        editboxes(n).Value = loadedValue;
        if strcmp(editboxes(n).Tag,'headingMinMaxNstep')
            loadedValue(abs(loadedValue)<1e-10) = 0; % handles the eps cases
            unsignedNonzero = unique(abs(loadedValue(loadedValue~=0)));
            editboxes(n).String = [num2str(min(unsignedNonzero)) ' ' num2str(max(unsignedNonzero)) ' ' num2str(length(unsignedNonzero))];
            % grab the extras, if any
            largestValInd = find(loadedValue==unsignedNonzero(end));
            extra = 0;
            while largestValInd + extra < length(loadedValue)
                extra = extra + 1;
                editboxes(n).String = ([editboxes(n).String ' ' num2str(loadedValue(largestValInd+extra))]);
            end
        else
            editboxes(n).String = num2str(loadedValue);
        end
    catch
        editboxes(n).Value = NaN;
        editboxes(n).String = 'NaN';
    end
end
guidata(hObject,handles);

% Update handles structure
guidata(hObject, handles);

viewdist = handles.viewdist.Value;
heightcm = handles.heightcm.Value;
set(handles.fovy_text, 'String',sprintf('Field of view: %.2f',2*atand((heightcm/2)/viewdist)))

% UIWAIT makes rigBgui_advsettings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rigBgui_advsettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN CreateFcns

% General
function viewdist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interocularDist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function requiredFixFraction_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corrloop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Motion parameters
function headingDur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function headingAmpl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function headingSigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Visual display/stimulus

function widthcm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function heightcm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fovy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clipNear_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clipFar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function visclipFar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eyeHeight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bgColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Dot features
function dotNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dotSizeinWorld_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dotColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minDotLifetime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxDotLifetime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function apertureDiam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% END CreateFcns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)

paradigm = getappdata(0,'paradigm_value');
subject  = getappdata(0,'subject_value');

% load the last saved settingsStruct...
load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
% ...but then overwrite loaded values using GUI data.
% this is just the reverse of the process from loadSettings

% here we don't pass in hObject, so that we save ALL the settings, from the
% main and advanced GUIs
checkboxes = findobj('Style','checkbox');
for n = 1:length(checkboxes)
    try
        eval([checkboxes(n).TooltipString ' = ' num2str(checkboxes(n).Value) ';']);
    catch
        eval([checkboxes(n).TooltipString ' = 0;']);
    end
end

editboxes = findobj('Style','edit');
for n = 1:length(editboxes)
    try
        eval([editboxes(n).TooltipString ' = [' num2str(editboxes(n).Value) '];']);
    catch
        editboxes(n).Value = NaN;
    end
end

% now save it back to the .mat file
save(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat'],'settingsStruct');

d = dialog('Position',[800 500 250 150],'Name','Done');
txt = uicontrol('Parent',d,...
       'Style','text',...
       'Position',[20 70 210 40],...
       'String','Settings saved!');
btn = uicontrol('Parent',d,...
       'Position',[95 60 60 25],...
       'String','OK',...
       'Callback','delete(gcf)');   
   
%%%%%%%%%%%%%%%%%%%% CHECK BOXES %%%%%%%%%%%%%%%%%%%%

function stereomode_Callback(hObject, eventdata, handles)

function isDots_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%% EDIT BOXES %%%%%%%%%%%%%%%%%%%%

% General

function viewdist_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

viewdist = handles.viewdist.Value;
heightcm = handles.heightcm.Value;
set(handles.fovy_text, 'String',sprintf('Field of view: %.2f',2*atand((heightcm/2)/viewdist)))
guidata(hObject, handles);

function interocularDist_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function requiredFixFraction_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function corrloop_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

% Motion parameters

function headingDur_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingSigma_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingAmpl_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);


% Visual stimulus parameters

function widthcm_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function heightcm_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

viewdist = handles.viewdist.Value;
heightcm = handles.heightcm.Value;
set(handles.fovy_text, 'String',sprintf('Field of view:\t%.2f',2*atand((heightcm/2)/viewdist)))
guidata(hObject, handles);

function clipNear_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function clipFar_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);


function visclipFar_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function eyeHeight_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function bgColor_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

% dot features
function dotNumber_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function dotSizeinWorld_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function dotColor_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function minDotLifetime_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function maxDotLifetime_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function apertureDiam_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);
