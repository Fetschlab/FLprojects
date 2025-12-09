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

% Last Modified by GUIDE v2.5 30-May-2019 13:48:17

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

h = msgbox('Don''t forget to connect to Motion Platform in Ubuntu menu bar (double-arrow icon)');
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)+300 pos(3) pos(4)]);

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



% CF: Stack up all the CreateFcns here to make reading the rest easier

function paradigm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function subject_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function amountRewardLowConf_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function amountRewardHighConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function penaltyErrorLowConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function penaltyErrorHighConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function probOneConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function whichOneConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function confTask_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function headingModality_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function headingDelta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function headingMinMaxNstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function headingMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function headingNstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function numReps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function textFeedback_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dotCoh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function numSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% END CreateFcns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in paradigm.
function paradigm_Callback(hObject, eventdata, handles)

% --- Executes on selection change in subject.
function subject_Callback(hObject, eventdata, handles)

% --- Executes on button press in activate.
function activate_Callback(hObject, eventdata, handles)
mpStartup;

% --- Executes on button press in deactivate.
function deactivate_Callback(hObject, eventdata, handles)
global mp
mp.deactivate

% --- Executes on button press in park.
function park_Callback(hObject, eventdata, handles)
global mp
if mp.parked==0
    parkToLoad_gui;
    mp.parked = 1;
end

% --- Executes on button press in returnToZero.
function returnToZero_Callback(hObject, eventdata, handles)

global mp
if size(mp, 1)==1
    mp.returnToZero; WaitSecs(0.2); mp.returnToZero;
    % returnToZero sets parked=0
else
    error('MP not connected/activated');
    mp.parked=0;
end


% --- Executes on button press in clearTemp.
function clearTemp_Callback(hObject, eventdata, handles)

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

subject = handles.subject.String{handles.subject.Value};
command = ['sudo scp /home/fetschlab/data/' subject '* fetschlab@172.30.3.33:/volume1/data/' subject];
system(command,'-echo');
% NEXT: add this as a dialogue after run ends


% --- Executes on button press in loadSettings.
function loadSettings_Callback(hObject, eventdata, handles)

paradigm = handles.paradigm.String{handles.paradigm.Value};
subject = handles.subject.String{handles.subject.Value};
load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);

% identify handles (GUI fields) to update, then grab their values from the
% local settingsStruct (location stored in TooltipString)
checkboxes = findobj('Style','checkbox');
for n = 1:length(checkboxes)
    try
        eval(['loadedValue = ' checkboxes(n).TooltipString ';']);
        checkboxes(n).Value = loadedValue;
    catch
        checkboxes(n).Value = 0;
    end
end

editboxes = findobj('Style','edit');
for n = 1:length(editboxes)
    try
        eval(['loadedValue = ' editboxes(n).TooltipString ';']);
        editboxes(n).Value = loadedValue;
        if strcmp(editboxes(n).Tag,'headingMinMaxNstep')
            unsignedNonzero = unique(abs(loadedValue(loadedValue~=0)));
            editboxes(n).String = [num2str(min(unsignedNonzero)) ' ' num2str(max(unsignedNonzero)) ' ' num2str(length(unsignedNonzero))];
        else        
            editboxes(n).String = num2str(loadedValue);
        end
    catch
        editboxes(n).Value = NaN;
        editboxes(n).String = 'NaN';
    end
end
% Update handles structure
% % % handles.settingsStruct = settingsStruct;
guidata(hObject,handles);


% --- Executes on button press in saveSettings.
function saveSettings_Callback(hObject, eventdata, handles) %#ok<*INUSL>

paradigm = handles.paradigm.String{handles.paradigm.Value};
subject = handles.subject.String{handles.subject.Value};

% load the last saved settingsStruct...
load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
% ...but then overwrite loaded values using GUI data.
% this is just the reverse of the process from loadSettings
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
    
% % %     % seems like values are not loaded unless editbox is used, even after
% % %     % load settings???
% % %     if ~strcmp(paradigm,'SaccadeTraining') && all(isnan(editboxes(n).Value));
% % %         keyboard
% % %     end
end

% now save it back to the .mat file
save(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat'],'settingsStruct');

% % why the f- does this not work
% try
% d = dialog('Position',[800 500 250 150],'Name','CAREFUL');
% txt = uicontrol('Parent',d,...
%        'Style','text',...
%        'Position',[20 80 210 40],...
%        'String','Overwrite settings in .mat file? Are you sure?');
% btn3 = uicontrol('Parent',d,...
%        'Position',[60 20 60 25],...
%        'String','Yes',...
%        'Callback','save([''/home/fetschlab/PLDAPSFL/settingsStruct_'' paradigm ''_'' subject ''.mat''],''settingsStruct''); delete(gcf);');
% btn4 = uicontrol('Parent',d,...
%        'Position',[125 20 60 25],...
%        'String','No',...
%        'Callback','delete(gcf)');
% catch me
%     keyboard
% end

d = dialog('Position',[800 500 250 150],'Name','Done');
txt = uicontrol('Parent',d,...
       'Style','text',...
       'Position',[20 70 210 40],...
       'String','Settings saved!');
btn = uicontrol('Parent',d,...
       'Position',[95 60 60 25],...
       'String','OK',...
       'Callback','delete(gcf)');
   

   
   % --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

sca;
% reconnect/activate, or else occasionally platform won't move
mpStartup; pause(0.2);

paradigm = handles.paradigm.String{handles.paradigm.Value};
subject = handles.subject.String{handles.subject.Value};

% load the last saved settingsStruct...
load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
% ...but then overwrite loaded values using GUI data.
% this is just the reverse of the process from loadSettings
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
    
% % %     % seems like values are not loaded unless editbox is used, even after
% % %     % load settings???
% % %     if ~strcmp(paradigm,'SaccadeTraining') && all(isnan(editboxes(n).Value));
% % %         keyboard
% % %     end
    
end

eval(['p = pldaps(@' paradigm '_setup, ''' subject ''', settingsStruct);']);
p.run;



%%%%%%%%%%%%%%%%%%%% EDIT BOXES %%%%%%%%%%%%%%%%%%%%

function amountRewardLowConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.reward.amountRewardLowConf = str2num(get(hObject,'String')); %#ok<*ST2NM>
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function amountRewardHighConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.reward.amountRewardHighConf = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function penaltyErrorLowConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.ITI.penaltyErrorLowConf = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function penaltyErrorHighConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.ITI.penaltyErrorHighConf = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function probOneConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.behavior.probOneConf = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function whichOneConf_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.behavior.whichOneConf = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

% --- Executes on button press in confTask.
function confTask_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingModality_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.stimulus.headingModality = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingDelta_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.stimulus.headingDelta = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingMinMaxNstep_Callback(hObject, eventdata, handles)
try
    hdgVect = str2num(get(hObject,'String'));
    hdg = logspace(log10(hdgVect(1)),log10(hdgVect(2)),hdgVect(3));
% % %     handles.settingsStruct.stimulus.headingTheta = [-fliplr(hdg) 0 hdg];
    set(hObject,'Value',[-fliplr(hdg) 0 hdg]);
catch
    error('Heading Range must be specified as [min max nstep]');
end
guidata(hObject, handles);

function numReps_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.stimulus.numReps = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function textFeedback_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.behavior.textFeedback = str2num(get(hObject,'String'));
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function dotCoh_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function numSteps_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%% CHECKBOXES %%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in pass.
function pass_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.pldaps.pass = get(hObject,'Value');
% % % guidata(hObject, handles);

% --- Executes on button press in dummyMP.
function dummyMP_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.pldaps.dummyMP = get(hObject,'Value');
% % % guidata(hObject, handles);

% --- Executes on button press in overlay.
function overlay_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.display.useOverlay = get(hObject,'Value');
% % % guidata(hObject, handles);

% --- Executes on button press in imu.
function imu_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.pldaps.useIMU = get(hObject,'Value');
% % % guidata(hObject, handles);

% --- Executes on button press in nexonar.
function nexonar_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.pldaps.useNexonar = get(hObject,'Value');
% % % guidata(hObject, handles);

% --- Executes on button press in RTtask.
function RTtask_Callback(hObject, eventdata, handles)
% % % handles.settingsStruct.behavior.RTtask = get(hObject,'Value');
% % % guidata(hObject, handles);


