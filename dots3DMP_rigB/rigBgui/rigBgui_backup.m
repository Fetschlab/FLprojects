function varargout = rigBgui(varargin)
% RIGBGUI MATLAB code for rigBgui.fig
%      RIGBGUI, by itself, creates a new RIGBGUI or raises the existing
%      singleton*.
%
%      H = RIGBGUI returns the handle to a new RIGBGUI or the handle to
%      the existing singleton*.
%
%      RIGBGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIGBGUI.M with the given input arguments.
%
%      RIGBGUI('Property','Value',...) creates a new RIGBGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rigBgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rigBgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rigBgui

% Last Modified by GUIDE v2.5 15-May-2019 09:28:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rigBgui_OpeningFcn, ...
                   'gui_OutputFcn',  @rigBgui_OutputFcn, ...
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


% --- Executes just before rigBgui is made visible.
function rigBgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rigBgui (see VARARGIN)

msgbox('Don''t forget to connect to Motion Platform in Ubuntu menu bar (double-arrow icon)');

% Choose default command line output for rigBgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes rigBgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rigBgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the deactivate flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to deactivate the data.

sca; clear global mp;


% --- Executes on selection change in paradigm.
function paradigm_Callback(hObject, eventdata, handles)
% hObject    handle to paradigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns paradigm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from paradigm

evalin(


% --- Executes during object creation, after setting all properties.
function paradigm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paradigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in subject.
function subject_Callback(hObject, eventdata, handles)
% hObject    handle to subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in activate.
function activate_Callback(hObject, eventdata, handles)
% hObject    handle to activate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mpStartup;


% --- Executes on button press in deactivate.
function deactivate_Callback(hObject, eventdata, handles)
% hObject    handle to deactivate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mp
mp.deactivate


% --- Executes on button press in park.
function park_Callback(hObject, eventdata, handles)
% hObject    handle to park (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mp
if mp.parked==0
    parkToLoad_gui;
    mp.parked = 1;
end


% --- Executes on button press in returnToZero.
function returnToZero_Callback(hObject, eventdata, handles)
% hObject    handle to returnToZero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mp
if size(mp, 1)==1
    mp.returnToZero; WaitSecs(0.2); mp.returnToZero;
    % returnToZero sets parked=0
else
    error('MP not connected/activated');
    mp.parked=0;
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paradigm = handles.paradigm.String{handles.paradigm.Value};
subject = handles.subject.String{handles.subject.Value};
load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
eval(['p = pldaps(@' paradigm '_setup, ''' subject ''', settingsStruct);']);
p.run;


% % % % --- Executes on button press in stop.
% % % function stop_Callback(hObject, eventdata, handles)
% % % % hObject    handle to stop (see GCBO)
% % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % handles    structure with handles and user data (see GUIDATA)
% % % 
% % % % can't execute code in keyboardCmd 'q' because no access to p
% % % 
% % % % % this also doesn't work (robot keypresses not detected by KbQueueCheck
% % % % 
% % % % %Initialize the java engine 
% % % % import java.awt.*;
% % % % import java.awt.event.*;
% % % % %Create a Robot-object to do the key-pressing
% % % % rob=Robot;
% % % % for n = 1:5
% % % %     commandwindow;
% % % %     rob.keyPress(KeyEvent.VK_Q);
% % % %     rob.keyRelease(KeyEvent.VK_Q);
% % % %     WaitSecs(0.1);
% % % % end


% --- Executes on button press in clearTemp.
function clearTemp_Callback(hObject, eventdata, handles)
% hObject    handle to clearTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = dialog('Position',[800 500 250 150],'Name','CAREFUL');
txt = uicontrol('Parent',d,...
       'Style','text',...
       'Position',[20 80 210 40],...
       'String','Clear TEMP folder? Are you sure?');
btn1 = uicontrol('Parent',d,...
       'Position',[60 20 60 25],...
       'String','Yes',...
       'Callback','delete /home/fetschlab/data/TEMP/*; delete(gcf);');
btn2 = uicontrol('Parent',d,...
       'Position',[125 20 60 25],...
       'String','No',...
       'Callback','delete(gcf)');


% --- Executes on button press in copyData.
function copyData_Callback(hObject, eventdata, handles)
% hObject    handle to copyData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subject = handles.subject.String{handles.subject.Value};
command = ['sudo scp /home/fetschlab/data/' subject '* fetschlab@172.30.3.33:/volume1/data/' subject];
system(command,'-echo');
% NEXT: add this as a dialogue after run ends


% --- Executes on button press in checkboxPass.
function checkboxPass_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxPass
settingsStruct.pldaps.pass = get(hObject,'Value');
set(hObject,'Value',settingsStruct.pldaps.pass);




function editPass_Callback(hObject, eventdata, handles)
% hObject    handle to editPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPass as text
%        str2double(get(hObject,'String')) returns contents of editPass as a double
settingsStruct.pldaps.pass = str2double(get(hObject,'String'));
set(hObject,'String',num2str(settingsStruct.pldaps.pass));


% --- Executes during object creation, after setting all properties.
function editPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
