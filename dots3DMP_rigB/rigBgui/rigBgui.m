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

% Last Modified by GUIDE v2.5 05-Sep-2024 15:56:55

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
set(h,'Position',[pos(1) pos(2)+200 pos(3) pos(4)]);

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
function rewardErrorLowConf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rewardIncrement_CreateFcn(hObject, eventdata, handles)
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
function headingPhi_CreateFcn(hObject, eventdata, handles)
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
function hdgStair_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function contTarg_CreateFcn(hObject, eventdata, handles)
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
function headingStepSD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fixReward_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fixRewardProb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function targetR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function targetTheta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'Backmodality = 1groundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function probOfMemorySaccade_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function timeTargDisappears_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function audioFeedback_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fixationPosition_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function colorCorrTarg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function randomizeReward_CreateFcn(hObject, eventdata, handles)
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
mpStartup(handles.dummyMP.Value);

% --- Executes on button press in deactivate.
function deactivate_Callback(hObject, eventdata, handles)
global mp
if size(mp, 1)==0
    mpStartup(handles.dummyMP.Value);
end
mp.deactivate

% --- Executes on button press in park.
function park_Callback(hObject, eventdata, handles)
global mp
if size(mp, 1)==0
    mpStartup(handles.dummyMP.Value);
end
ishuman = strcmp(handles.subject.String{handles.subject.Value},'human');

if mp.parked==0
    parkToLoad_gui(ishuman);
    mp.parked = 1;
end

% --- Executes on button press in returnToZero.
function returnToZero_Callback(hObject, eventdata, handles)
global mp
if size(mp, 1)==0
    mpStartup(handles.dummyMP.Value);
end
mp.returnToZero; WaitSecs(0.1); mp.returnToZero;

function mountPosition_Callback(hObject, eventdata, handles)
global mp
if size(mp, 1)==0
    mpStartup(handles.dummyMP.Value);
end
ishuman = strcmp(handles.subject.String{handles.subject.Value},'human');

if mp.mountPosition==0
    mountPosition_gui(ishuman);
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
command = ['sudo scp -P 3333 /home/fetschlab/data/' subject '* fetschlab@172.30.3.33:/var/services/homes/fetschlab/data/' subject];
%system(command,'-echo');

% trying to make a confirmation button, like the delete temp version
% above, but can't reference 'command' or 'subject' in the callback...
d = dialog('Position',[800 500 250 150],'Name','Confirm copy');
txt = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 80 210 40],...
    'String',['Copy all ' subject ' files from /data to NAS? Are you sure?']);

btn1 = uicontrol('Parent',d,...
    'Position',[60 20 60 25],...
    'String','Yes',...
    'Callback', @(hObject, eventdata) copyDatabtn1_Callback(hObject, eventdata, handles));
btn2 = uicontrol('Parent',d,...
    'Position',[125 20 60 25],...
    'String','No',...
    'Callback','delete(gcf);');
    
% --- Executes on button press in Yes in copyData.  
function copyDatabtn1_Callback(hObject, eventdata, handles)
subject = handles.subject.String{handles.subject.Value};
command = ['sudo scp -P 3333 /home/fetschlab/data/' subject '* fetschlab@172.30.3.33:/var/services/homes/fetschlab/data/' subject];
delete(gcf);
system(command,'-echo');



% --- Executes on button press in advanced_settings.
function advanced_settings_Callback(hObject, eventdata, handles)
setappdata(0,'subject_value',handles.subject.String{handles.subject.Value});
setappdata(0,'paradigm_value',handles.paradigm.String{handles.paradigm.Value});
rigBgui_advsettings;

% --- Executes on button press in advanced_trialTable.
function advanced_trialTable_Callback(hObject, eventdata, handles)
rigBgui_advTrialTable;



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
        %                 if strcmp(editboxes(n).Tag,'headingMinMaxNstep')
        %                     keyboard
        %                 end
        if strcmp(editboxes(n).Tag,'headingMinMaxNstep') && any(strcmp(handles.paradigm.String{handles.paradigm.Value},{'dots3DMP','dots3DMPtuning','motion2DMP'}))
            loadedValue(abs(loadedValue)<1e-10) = 0; % handles the eps cases
            unsignedNonzero = unique(abs(loadedValue(loadedValue~=0)));
            editboxes(n).String = [num2str(min(unsignedNonzero)) ' ' num2str(max(unsignedNonzero)) ' ' num2str(length(unsignedNonzero))];
            % grab the extras, if any
            largestValInd = find(loadedValue==unsignedNonzero(end), 1, 'last' );
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


% --- Executes on button press in saveSettings.
function saveSettings_Callback(hObject, eventdata, handles) %#ok<*INUSL>

paradigm = handles.paradigm.String{handles.paradigm.Value};
subject = handles.subject.String{handles.subject.Value};

% load the last saved settingsStruct...
try
    load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
catch
    disp('Can''t find settings struct with this name, creating new one')
end

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

% compute number of trials based on selected conditions
function trCount_Callback(hObject, eventdata, handles)

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
    catch exception
        disp(exception.message);
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
guidata(hObject,handles);

hdgs    = sort(handles.headingMinMaxNstep.Value);
mods    = sort(handles.headingModality.Value);
cohs    = sort(handles.dotCoh.Value);
deltas  = sort(handles.headingDelta.Value);
numReps = sort(handles.numReps.Value);


% display heading angles in text
set(handles.headings_text,'String',sprintf('%.1f  ',hdgs))

% and display total number of trial conditions

% compute conditions (trial types)
runByGui = 1;
[~,hdg] = generateCondList_3DMP([],runByGui,hdgs,mods,cohs,deltas,numReps);

nTr = length(hdg);

textstr = sprintf('%d conds\nx %d reps \n= %d TOTAL',nTr,numReps,nTr*numReps);
set(handles.trialcount,'String',textstr);


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

% clear PTB screens just in case one was leftover
sca;

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
end

if settingsStruct.behavior.audioFeedback > 0
    settingsStruct.sound.use = 1;
    settingsStruct.sound.usePsychPortAudio = 1;
end % otherwise, leave existing setting alone

if settingsStruct.mouse.useAsEyepos==1
    settingsStruct.mouse.use = 1;
    settingsStruct.pldaps.draw.eyepos.use = 1; % no point in mouse eyepos if can't see it
end % otherwise, leave existing setting alone

% store the settings and trajectories in a protocol file, if desired
if settingsStruct.pldaps.saveProtFile==1
    if strcmp(paradigm,'dots3DMP') || strcmp(paradigm,'motion2DMP') || strcmp(paradigm,'dots3DMPtuning')
        settingsStruct.pldaps.protFilename = ['/home/fetschlab/dots3DMP/protFiles/' paradigm ...
            '_cond_' handles.headingModality.String ...
            '_hdg_' handles.headingMinMaxNstep.String ...
            '_coh_' handles.dotCoh.String ...
            '_delta_' handles.headingDelta.String ...
            '_nSteps_' handles.numSteps.String ...
            '_stdev_' handles.headingStepSD.String ...
            '_hdgAmpl' num2str(settingsStruct.stimulus.headingAmpl) ...
            '_hdgDur' num2str(settingsStruct.stimulus.headingDur) ...
            '_hdgSigma' num2str(settingsStruct.stimulus.headingSigma) ...
            '.mat'];
    else
        %         settingsStruct.pldaps.protFilename = ['/home/fetschlab/dots3DMP/protFiles/' paradigm ...
        %             '_cond_' handles.headingModality.String ...
        %             '_hdg_' handles.headingMinMaxNstep.String ...
        %             '_coh_' handles.dotCoh.String ...
        %             '_delta_' handles.headingDelta.String ...
        %             '_nSteps_' handles.numSteps.String ...
        %             '_stdev_' handles.headingStepSD.String ...
        %             '.mat'];
        
        % just save super basic protFileName for e.g. mapping paradigms, will be overwritten
        settingsStruct.pldaps.protFilename = ['/home/fetschlab/dots3DMP/protFiles/' paradigm '.mat'];
    end
else
    settingsStruct.pldaps.protFilename = [];
end

if strcmp(subject,'human') && (strcmp(paradigm,'dots3DMP') || strcmp(paradigm,'motion2DMP') || strcmp(paradigm,'TargetDots3DMP'))
    fprintf('\n******\n');
    subjectTag = input('Enter 3-letter subject code: ','s');
else
    subjectTag = [];
end

% trCount_Callback(hObject, eventdata, handles); % run this automatically?

eval(['p = pldaps(@' paradigm '_setup, ''' subject subjectTag ''', settingsStruct);']);

p.run;




% % *****************************************************************
% % function Stop_Callback(hObject, eventdata, handles)
% %
% %   disp('STOP not yet implemented');
% %   return;
% %
% % Should be possible, using flags, to change Run button to Stop button when
% % it's pressed, and vice versa.
% % p.run needs to know when this button is pressed to quit the while loop
% % *****************************************************************



% start a task immediately after eye-cal
if strcmp(paradigm,'EyeCalibration') && (settingsStruct.behavior.start3DMPAfterEyeCal || settingsStruct.behavior.startStepsAfterEyeCal)
    % need to essentially repeat the code from loadData callback. can't
    % simply call that function, because we need to determine paradigm from
    % the checkboxes not the 'select' list
    
    settingsStruct.behavior.useEyeCal = 1;
    
    if settingsStruct.behavior.start3DMPAfterEyeCal && ~settingsStruct.behavior.startStepsAfterEyeCal
        paradigm = 'dots3DMP';
    elseif ~settingsStruct.behavior.start3DMPAfterEyeCal && settingsStruct.behavior.startStepsAfterEyeCal
        paradigm = 'motion2DMP';
    else
        error('cannot check both startXAfterEyeCal boxes');
    end
    subject = handles.subject.String{handles.subject.Value};
    load(['/home/fetschlab/PLDAPSFL/settingsStruct_' paradigm '_' subject '.mat']);
    
    
    % load the settings and trajectories from protocol file, if desired
    if settingsStruct.pldaps.saveProtFile==1
        settingsStruct.pldaps.protFilename = ['/home/fetschlab/dots3DMP/protFiles/' paradigm ...
            '_cond_' handles.headingModality.String ...
            '_hdg_' handles.headingMinMaxNstep.String ...
            '_coh_' handles.dotCoh.String ...
            '_delta_' handles.headingDelta.String ...
            '_nSteps_' handles.numSteps.String ...
            '_stdev_' handles.headingStepSD.String ...
            '.mat'];
    else
        settingsStruct.pldaps.protFilename = [];
    end
    
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
    
    eval(['p = pldaps(@dots3DMP_setup, ''' subject ''', settingsStruct);']);
    p.run;
    
end




%%%%%%%%%%%%%%%%%%%% EDIT BOXES %%%%%%%%%%%%%%%%%%%%

function amountRewardLowConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function amountRewardHighConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function rewardErrorLowConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function rewardIncrement_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function penaltyErrorLowConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function penaltyErrorHighConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function probOneConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function whichOneConf_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function confTask_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

% function delay2Conf_Callback(hObject, eventdata, handles)
% set(hObject,'Value',str2num(get(hObject,'String')));
% guidata(hObject, handles);

function headingModality_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingDelta_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingTheta_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingPhi_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

% this one needs special processing to get vector of heading angles
function headingMinMaxNstep_Callback(hObject, eventdata, handles)
if any(strcmp(handles.paradigm.String{handles.paradigm.Value},{'dots3DMP','dots3DMPtuning','motion2DMP'}))
    try
        hdgVect = str2num(get(hObject,'String'));
        hdgVect(hdgVect==0) = eps; % to avoid log(0)
        hdg = logspace(log10(hdgVect(1)),log10(hdgVect(2)),hdgVect(3));
        %hdg = [hdg(1) hdg ];  % changed by YC 20240117, manully increase the trial num of hdg(1) and hdg(2)
        hdg = [-fliplr(hdg) hdg];
        if length(hdgVect)>3 % grab optional extras, including zero
            hdg = [hdg hdgVect(4:end)];
        end
        hdg(abs(hdg)<0.01) = 0; % change back to zero!!
        set(hObject,'Value',hdg);
        trCount_Callback(hObject, eventdata, handles)
    catch
        error('Heading Range must be specified as [min max nstep [optional extras]]');
    end
else
    set(hObject,'Value',str2num(get(hObject,'String')));
end
guidata(hObject, handles);

function numReps_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function hdgStair_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function textFeedback_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function dotCoh_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function numSteps_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function headingStepSD_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function fixReward_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function fixRewardProb_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function targetR_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function targetTheta_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function probOfMemorySaccade_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function timeTargDisappears_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function audioFeedback_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function fixationPosition_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function colorCorrTarg_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

function randomizeReward_Callback(hObject, eventdata, handles)
set(hObject,'Value',str2num(get(hObject,'String')));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%% CHECKBOXES %%%%%%%%%%%%%%%%%%%%


function pass_Callback(hObject, eventdata, handles)

function dummyMP_Callback(hObject, eventdata, handles)

function rewardIncModsIndep_Callback(hObject, eventdata, handles)

function overlay_Callback(hObject, eventdata, handles)

function inRig_Callback(hObject, eventdata, handles)

function imu_Callback(hObject, eventdata, handles)

function nexonar_Callback(hObject, eventdata, handles)

function ripple_Callback(hObject, eventdata, handles)

function RTtask_Callback(hObject, eventdata, handles)

function contTarg_Callback(hObject, eventdata, handles)

function saveProtFile_Callback(hObject, eventdata, handles)

function useEyeCal_Callback(hObject, eventdata, handles)

function start3DMPAfterEyeCal_Callback(hObject, eventdata, handles)

function startStepsAfterEyeCal_Callback(hObject, eventdata, handles)

function mpMask_Callback(hObject, eventdata, handles)

function mouseUseAsEyepos_Callback(hObject, eventdata, handles)

function hotSpotCorrection_Callback(hObject, eventdata, handles)

function mouseStimPos_Callback(hObject, eventdata, handles)

function dupVes_Callback(hObject, eventdata, handles)

function neuropixels_Callback(hObject, eventdata, handles)
