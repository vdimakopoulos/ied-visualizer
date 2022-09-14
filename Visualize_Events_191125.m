function varargout = Visualize_Events_191125(varargin)
% VISUALIZE_EVENTS_191125 MATLAB code for Visualize_Events_191125.fig
%      VISUALIZE_EVENTS_191125, by itself, creates a new VISUALIZE_EVENTS_191125 or raises the existing
%      singleton*.
%
%      H = VISUALIZE_EVENTS_191125 returns the handle to a new VISUALIZE_EVENTS_191125 or the handle to
%      the existing singleton*.
%
%      VISUALIZE_EVENTS_191125('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZE_EVENTS_191125.M with the given input arguments.
%
%      VISUALIZE_EVENTS_191125('Property','Value',...) creates a new VISUALIZE_EVENTS_191125 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Visualize_Events_191125_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Visualize_Events_191125_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Visualize_Events_191125

% Last Modified by GUIDE v2.5 02-Mar-2020 13:02:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Visualize_Events_191125_OpeningFcn, ...
    'gui_OutputFcn',  @Visualize_Events_191125_OutputFcn, ...
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


% --- Executes just before Visualize_Events_191125 is made visible.
function Visualize_Events_191125_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Visualize_Events_180704 (see VARARGIN)

%% Initialize figure
Initialize_Main_Figure(hObject)

%% Read input parameters
InputParams = varargin{1};
handles.slider2.Value = 0;

%% Check inputs % TODO
PlotParams = Check_Input_Params(InputParams);
PlotParams.strSaveImagesFolderPath = InputParams.strSaveImagesFolderName; % TODO
PlotParams.strPlotTitle = InputParams.strPlotTitle; % TODO

handles.MarkText.String = PlotParams.strPlotTitle;
handles.MarkText.FontSize = 13;
handles.MarkText.FontWeight = 'bold';
handles.MarkText.BackgroundColor = 'w';

if(~isfield(PlotParams,'SpikeMarkings'))
    PlotParams.SpikeMarkings = [];
end

%% Put into handles
handles.PlotParams = PlotParams;

%% Plot raw and filtered signals
handles = Plot_Raw_Signal(handles);

%% Take a screenshot of the whole signal
% Initialize screenshot number
handles.nScreenshot = 0;
pushbutton_save_image_Callback(handles.pushbutton_save_image,[],handles)
handles.nScreenshot = 1;

%% Initialize spike markings
if(~isfield(InputParams,'EventPosition'))
    handles.EventPosition = [];
    handles.nEvent = 0;
else
    handles.EventPosition = InputParams.EventPosition;
    handles.nEvent = size(handles.EventPosition,1);
    % If no third column is give, automatically assign group 1
    if(size(handles.EventPosition,2)==2)
        handles.EventPosition = [handles.EventPosition,ones(size(handles.EventPosition,1),1)];
    end
end

%% Save screenshots
if(~isfield(InputParams,'flagSaveScreenshots'))
    handles.flagSaveScreenshots = 0;
else
    handles.flagSaveScreenshots = InputParams.flagSaveScreenshots ;
end

%% Color of markings
handles.strSpikeColors = {'b','m','r','c',[0.9290,0.6940,0.1250],[0.4940,0.1840,0.5560]};

%% Initialize button for marking color
handles.nMarkingGroup = 1;
handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

%%
handles.pushbotton_color_boxes(1) = handles.pushbutton_color_1;
handles.pushbotton_color_boxes(2) = handles.pushbutton_color_2;
handles.pushbotton_color_boxes(3) = handles.pushbutton_color_3;
handles.pushbotton_color_boxes(4) = handles.pushbutton_color_4;
handles.pushbotton_color_boxes(5) = handles.pushbutton_color_5;
handles.pushbotton_color_boxes(6) = handles.pushbutton_color_6;

%% Initialize colors of pushbuttons
for nMarkingGroupColor = 1:length(handles.strSpikeColors)
    handles.pushbotton_color_boxes(nMarkingGroupColor).BackgroundColor = handles.strSpikeColors{nMarkingGroupColor};
end

%% Plot input markings on the signals (if there are any)
for nEvent = 1:handles.nEvent
    nChannelSpike = handles.EventPosition(nEvent,1);
    t_spike = handles.EventPosition(nEvent,2);
    nEventGroup = handles.EventPosition(nEvent,3);
    
    iChannelSpike = find(handles.PlotParams.ListOfChannels_ToPlot==nChannelSpike);
    
    if(~isempty(iChannelSpike))
        [~,ind_t_axis] = min(abs(handles.PlotParams.t-t_spike));
        ind_t_axis_range = (ind_t_axis-2000*0.25):(ind_t_axis+2000*0.25);
        
        hold on
        if(nEventGroup>0)
            handles.SpikePlots{nEvent} = plot(handles.PlotParams.t(ind_t_axis_range),repmat(-(iChannelSpike-1)*handles.PlotParams.YShift(1),length(ind_t_axis_range),1)+...
                handles.PlotParams.YShift(1)*(rand(1)*0.1-0.05),'.m');
            handles.SpikePlots{nEvent}.Color = handles.strSpikeColors{nEventGroup};
        else
            handles.SpikePlots{nEvent} = [];
        end
    else
        handles.SpikePlots{nEvent} = [];
    end
end

%% Set initial xlims % TODO
xlim([0,PlotParams.tWindow(1)])
handles.PlotParams.t = double(handles.PlotParams.t); % TODO
t_max = handles.PlotParams.t(end);
SlideStep = handles.PlotParams.tWindow(1)/t_max;
set(handles.slider2,'SliderStep',[SlideStep,SlideStep])
SliderValue = get(handles.slider2,'Value');
SliderSteps = get(handles.slider2,'SliderStep');
SliderLimits = [SliderValue,SliderValue+SliderSteps(2)];
xlimMax = handles.PlotParams.t(end);
xlimVals = xlimMax*SliderLimits;
set(handles.axes_raw,'XLim',xlimVals)

%%
% set(handles.MarkText,'String',''); % TODO what to do with this

%% Figure position style?
% set(handles.figure1,'PaperUnits','centimeters');
% set(handles.figure1,'PaperPosition',[-10,0,10,15]); %x_width=10cm y_width=15cm
% set(handles.figure1,'PaperPosition',[0,0,10,15]); %x_width=10cm y_width=15cm

%% Save image for whole signal
% strImageFileName = 'Whole_Interval';

% strSaveImagesFolderName = handles.strSaveImagesFolderName;

% set(handles.figure1,'PaperUnits','centimeters');
% set(handles.figure1,'PaperPosition',[-10,0,10,15]); %x_width=10cm y_width=15cm
% set(handles.figure1,'PaperPosition',[0,0,10,15]); %x_width=10cm y_width=15cm

% export_fig(fig,[strSaveImagesFolderName,strImageFileName,'.','png'],['-','png'],'-transparent','-r100')
% saveas(handles.figure1,[strSaveImagesFolderName,strImageFileName],'bmp');

% set(handles.figure1,'PaperUnits','centimeters');
% set(handles.figure1,'PaperPosition',[-10,0,10,15]); %x_width=10cm y_width=15cm

%% If save screenshots mode is selected, save every 15 s window then close
if(handles.flagSaveScreenshots)
    handles.nScreenshot = 1;
    pushbutton_save_image_Callback(handles.pushbutton_save_image,[],handles)
    for nRun = 1:floor(1/handles.slider2.SliderStep(2))
        handles.slider2.Value =  handles.slider2.Value+handles.slider2.SliderStep(2);
        slider2_Callback(handles.slider2,[],handles)
        handles.nScreenshot = nRun+1;
        pushbutton_save_image_Callback(handles.pushbutton_save_image,[],handles)
    end
    %     Visualize_Events_191125_OutputFcn([],[],handles)
    return;
    pushbutton_quit_Callback(handles.pushbutton_quit,[],handles)
end

%% Choose default command line output for Visualize_Events_191125
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% % UIWAIT makes Visualize_Events_180704 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Visualize_Events_191125_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% h = findobj('Tag','figure1');
% nEventValidity = handles.nEventValidity;
EventPosition = handles.EventPosition;

if(~isempty(EventPosition))
%     EventPosition(EventPosition(:,3)==-1,:) = [];
    [~,indSort] = sort(EventPosition(:,2));
    EventPosition = EventPosition(indSort,:);
end

% Get default command line output from handles structure
varargout{1} = EventPosition;

% The figure can be deleted now
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(~, ~, ~)
% hObject    handle to the selected object in uibuttongroup1
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% class = get(handles.uibuttongroup1.SelectedObject,'Tag');
% handles.ClassifiedVisually(handles.current_event) = str2double(class(end));

% guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
% function ReverseButton_CreateFcn(~, ~, ~)
% hObject    handle to ReverseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);


function [ handles ] = Plot_Raw_Signal(handles)

% Plot raw signal
nBand_ToPlot = 1;
% Plot parameters
PlotParams = handles.PlotParams;
strPlotColor = 'k';
% Main axis
axMain = handles.axes_raw;
% Clear axis
axes(axMain)
cla

% YLim for the axis
YLim_Ax = [-(PlotParams.YMargin(2)+length(PlotParams.ListOfChannels_ToPlot)-1),...
    PlotParams.YMargin(1)]*PlotParams.YShift(nBand_ToPlot);

%% Plot signals with a y-shift for each axis
plSingleChannel = zeros(length(PlotParams.ListOfChannels_ToPlot),1);
ranYShift = zeros(length(PlotParams.ListOfChannels_ToPlot),1);
for iChannel_ToPlot = 1:length(PlotParams.ListOfChannels_ToPlot)
    nChannel = PlotParams.ListOfChannels_ToPlot(iChannel_ToPlot);
    
    % Shift the signal
    YShift_SingleChannel = -(iChannel_ToPlot-1)*PlotParams.YShift(nBand_ToPlot);
    % Save for ylabels
    ranYShift(iChannel_ToPlot) = YShift_SingleChannel;
    
    % Plot signal
    plSingleChannel(iChannel_ToPlot) = plot(PlotParams.t,PlotParams.dataAll{nBand_ToPlot}(nChannel,:)+YShift_SingleChannel,strPlotColor);
    hold on
end

% Plot the y-axis to show signal amplitude
pause(1)
yyaxis right
plot(PlotParams.t,NaN(1,length(PlotParams.t)));

% XLims
xlim([0,PlotParams.t(end)])
% XLabel
xlabel('Time (s)')

% YAxis
axMain.YAxis(1).Limits = YLim_Ax;
% Link left and right yaxes
linkprop([axMain.YAxis(1) axMain.YAxis(2)],'Limits');
linkprop([axMain.YAxis(1) axMain.YAxis(2)],'TickValues');
% Colors and font size for right axis
axMain.YAxis(2).Color = 0.4*[1,1,1];
axMain.YAxis(2).FontSize = 8.5;

% Electrode names on the left axis
strLabel = PlotParams.ElectrodeLabels(PlotParams.ListOfChannels_ToPlot);
strLabel = strrep(strLabel,'_','\_');
strLabel = strrep(strLabel,'\\_','\_');
axMain.YAxis(1).TickValues = sort(ranYShift);
axMain.YAxis(2).TickValues = sort(ranYShift);
axMain.YAxis(1).TickLabels = flipud(strLabel(:));

% Grids
grid on
grid minor
% Tick directions and length
set(axMain,'TickDir','out')
axMain.TickLength(1) = 0.001;

% Save axis and plots in handles
handles.axMain = axMain;
handles.pl = plSingleChannel;


% --- Executes on button press in pushbutton_close.
% function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% close(handles.gui)

% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_save_image.
function pushbutton_save_image_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nScreenshot = handles.nScreenshot;
fig = handles.figure1;
if(nScreenshot==0)
    strImageFileName = 'Selected_Channels_Whole_Interval';
else
    strImageFileName = sprintf('Screenshot_%.5d',nScreenshot);
end
% Main folder for screenshots
strSaveImagesFolderName = handles.PlotParams.strSaveImagesFolderPath;

% Paper position mode
fig.PaperPositionMode = 'auto';

% Save image
saveas(fig,[strSaveImagesFolderName,strImageFileName],'png');

% export_fig(fig,[strSaveImagesFolderName,strImageFileName,'.','png'],['-','png'],'-transparent','-r100')
% set(handles.figure1,'PaperUnits','centimeters');
% set(handles.figure1,'PaperPosition',[-10,0,10,15]); %x_width=10cm y_width=15cm
% set(handles.figure1,'PaperPosition',[0,0,10,15]); %x_width=10cm y_width=15cm

% Increment screenshot number
nScreenshot = nScreenshot+1;
handles.nScreenshot = nScreenshot;

guidata(hObject,handles);


% --- Executes on button press in pushbutton_quit.
function pushbutton_quit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SliderValue = get(hObject,'Value');
SliderSteps = get(handles.slider2,'SliderStep');
SliderLimits = [handles.SliderValue,handles.SliderValue+SliderSteps(2)];
xlimMax = handles.PlotParams.t(end);
xlimVals = xlimMax*SliderLimits;
set(handles.axes_raw,'XLim',xlimVals)
xlimValsMiddle = mean(xlimVals);
xlimVals2 = [xlimValsMiddle-0.5,xlimValsMiddle+0.5];
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PushbuttonMark.
function PushbuttonMark_Callback(hObject, eventdata, handles)
% hObject    handle to PushbuttonMark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginputx(1);

axLimX = get(gca,'XLim');
axLimY = get(gca,'YLim');

cond1 = (x<axLimX(1))||(x>axLimX(2));
cond2 = (y<axLimY(1))||(y>axLimY(2));

if(cond1||cond2)
    fprintf('\nOut of loop\n')
    return;
end

y_shiftRaw = handles.PlotParams.YShift(1);

iChannelSpike = floor(-(y-y_shiftRaw/2)/y_shiftRaw)+1;
nChannelSpike = handles.PlotParams.ListOfChannels_ToPlot(iChannelSpike);
t_spike = x;

handles.nEvent = handles.nEvent+1;

if(~isempty(handles.EventPosition))
    handles.EventPosition(handles.nEvent,:) = [nChannelSpike,t_spike,handles.nMarkingGroup];
else
    handles.EventPosition = [nChannelSpike,t_spike,handles.nMarkingGroup];
end
[~,ind_t_axis] = min(abs(handles.PlotParams.t-t_spike));
ind_t_axis_range = (ind_t_axis-2000*0.25):(ind_t_axis+2000*0.25);

hold on
handles.SpikePlots{handles.nEvent} = plot(handles.PlotParams.t(ind_t_axis_range),repmat(-(iChannelSpike-1)*y_shiftRaw,length(ind_t_axis_range),1)+...
    handles.PlotParams.YShift(1)*(rand(1)*0.1-0.05),'.m');
% Mark with the current marking group
handles.SpikePlots{handles.nEvent}.Color = handles.strSpikeColors{handles.nMarkingGroup};

while(1)
    [x,y] = ginputx(1);
    axLimX = get(gca,'XLim');
    axLimY = get(gca,'YLim');
    cond1 = (x<axLimX(1))||(x>axLimX(2));
    cond2 = (y<axLimY(1))||(y>axLimY(2));
    
    if(cond1||cond2)
        fprintf('\nOut of loop\n')
        break;
    end
    
    iChannelSpike = floor(-(y-y_shiftRaw/2)/y_shiftRaw)+1;
    nChannelSpike = handles.PlotParams.ListOfChannels_ToPlot(iChannelSpike);
    
    t_spike = x;
    handles.nEvent = handles.nEvent+1;
    
    handles.EventPosition(handles.nEvent,:) = [nChannelSpike,t_spike,handles.nMarkingGroup];
    
    [~,ind_t_axis] = min(abs(handles.PlotParams.t-t_spike));
    ind_t_axis_range = (ind_t_axis-2000*0.25):(ind_t_axis+2000*0.25);
    
    hold on
    handles.SpikePlots{handles.nEvent} = plot(handles.PlotParams.t(ind_t_axis_range),repmat(-(iChannelSpike-1)*y_shiftRaw,length(ind_t_axis_range),1)+...
        handles.PlotParams.YShift(1)*(rand(1)*0.1-0.05),'.m');
    handles.SpikePlots{handles.nEvent}.Color = handles.strSpikeColors{handles.nMarkingGroup};
    
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginputx(1);

axLimX = get(gca,'XLim');
axLimY = get(gca,'YLim');

cond1 = (x<axLimX(1))||(x>axLimX(2));
cond2 = (y<axLimY(1))||(y>axLimY(2));

if(cond1||cond2)
    fprintf('\nOut of loop\n')
    return;
end

y_shiftRaw = handles.PlotParams.YShift(1);

iChannelSpike = floor(-(y-y_shiftRaw/2)/y_shiftRaw)+1;
nChannelSpike = handles.PlotParams.ListOfChannels_ToPlot(iChannelSpike);

t_spike = x;

[~,ind_t_axis] = min(abs(handles.PlotParams.t-t_spike));
for ii = 1:size(handles.EventPosition,1)
    if(handles.EventPosition(ii,1)==nChannelSpike)
        t_spike_in_array = handles.EventPosition(ii,2);
        [~,ind_t_axis_in_array] = min(abs(handles.PlotParams.t-t_spike_in_array));
        ind_t_axis_range_in_array = (ind_t_axis_in_array-2000*0.25):(ind_t_axis_in_array+2000*0.25);
        cond = ~isempty(find(intersect(ind_t_axis_range_in_array,ind_t_axis),1));
        if(cond)
            if(handles.EventPosition(ii,3)~=-1)
                handles.EventPosition(ii,3) = -1;
                delete(handles.SpikePlots{ii})
            end
        end
    end
end

guidata(hObject,handles);


function Initialize_Main_Figure(hObject)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Visible figure
hObject.Visible = 'on';
% Display menubar and toolbar
hObject.MenuBar = 'figure';
hObject.ToolBar = 'figure';
% Pixels as units
hObject.Units = 'pixels';

% Get display positions
[~,nNumberOfMonitors] = Get_Monitor_Positions_Pixels();

% Display fullscreen figure in chosen monitor, default is monitor 2
nMonitor = min(nNumberOfMonitors,2);

% Make the figure fullscreen in the selected monitor
Make_Fullscreen_Figure_Selected_Monitor(nMonitor,hObject)

guidata(hObject);


function [MonitorPosition, nNumberOfMonitors] = Get_Monitor_Positions_Pixels()

% Get display positions
MonitorPosition = get(0,'MonitorPositions');
% Number of displays
nNumberOfMonitors = size(MonitorPosition,1);


function Make_Fullscreen_Figure_Selected_Monitor(nMonitor, hObject)

% Get display positions
[MonitorPosition,nNumberOfMonitors] = Get_Monitor_Positions_Pixels();

% Do nothing if nMonitor > nNumberOfMonitors
nMonitor = min(nNumberOfMonitors,nMonitor);

% Position of selected display
DisplayPosition = MonitorPosition(nMonitor,:);

% Set figure size to fullscreen in the chosen monitor
set(hObject,'Outerposition',DisplayPosition)

% Make sure that is actually fullscreen
FigureJavaFrame = get(hObject,'JavaFrame');
pause(1);
set(FigureJavaFrame,'Maximized',1);

guidata(hObject);


function [PlotParams] = Check_Input_Params(InputParams)

% Raw signal
if(~isfield(InputParams,'data'))
    error('Raw signal ''data'' must be an input!')
end
% Check size of raw signal
siz = size(InputParams.data);
if(siz(2)<siz(1))
    error('Size of ''data'' must be (number of channels)x(number of samples)!')
end
% Assign input signal
PlotParams.data = InputParams.data;

% Sampling frequency
if(~isfield(InputParams,'fs'))
    error('What is the sampling frequency?')
end
PlotParams.fs = InputParams.fs;

% Time axis
if(isfield(InputParams,'t'))
    if(length(InputParams.t)~=size(PlotParams.data))
        error('Sizes of input signal and time axis do not match!')
    end
    PlotParams.t = InputParams.t;
else
    PlotParams.t = 0:1/PlotParams.fs:(size(PlotParams.data,2)-1)/PlotParams.fs;
end

% Electrode labels
if(isfield(InputParams,'ElectrodeLabels'))
    PlotParams.ElectrodeLabels = InputParams.ElectrodeLabels;
else
    PlotParams.ElectrodeLabels = cellfun(@num2str,num2cell((1:size(InputParams.data,1))',ones(size(InputParams.data,1),1),1),'UniformOutput',0);
end

% Signal bands to plot
if(isfield(InputParams,'SignalBand_ToPlot'))
    if(isnumeric(InputParams.SignalBand_ToPlot))
        PlotParams.SignalBand_ToPlot = InputParams.SignalBand_ToPlot;
    end
else % Set using HFO bands to plot
end

% HFO bands to plot
strInputBandNames = {{'ripple'},{'fastripple','FR'}};
strInputBandNames = cellfun(@lower,strInputBandNames,'UniformOutput',0);
if(isfield(InputParams,'HFOBand_ToPlot'))
    if(isnumeric(InputParams.HFOBand_ToPlot))
        if(find(~ismember(InputParams.HFOBand_ToPlot,[1,2])))
            error('Specified band can take the values: 1 (ripple) or 2 (fast ripple)')
        end
        PlotParams.HFOBand_ToPlot = InputParams.HFOBand_ToPlot;
    elseif(ischar(InputParams.HFOBand_ToPlot))
        InputParams.HFOBand_ToPlot = lower(InputParams.HFOBand_ToPlot);
        InputParams.HFOBand_ToPlot = strrep(InputParams.HFOBand_ToPlot,' ','');
        if(find(~ismember(strrep(lower(InputParams.HFOBand_ToPlot),' ',''),[strInputBandNames{:}])))
            error('Specified band can take the values: ''ripple'',''FR'',''Fast Ripple''')
        end
        switch InputParams.HFOBand_ToPlot
            case strInputBandNames{1}
                PlotParams.HFOBand_ToPlot = 1;
            case strInputBandNames{2}
                PlotParams.HFOBand_ToPlot = 2;
        end
    elseif(iscell(InputParams.HFOBand_ToPlot))
        InputParams.HFOBand_ToPlot = lower(InputParams.HFOBand_ToPlot);
        InputParams.HFOBand_ToPlot = strrep(InputParams.HFOBand_ToPlot,' ','');
        if(find(~ismember(InputParams.HFOBand_ToPlot,[strInputBandNames{:}])))
            error('Specified band can take the values: ''ripple'',''FR'',''Fast Ripple''')
        end
        PlotParams.HFOBand_ToPlot = find(cellfun(@(x) ~isempty(find(x,1)),cellfun(@(x) ismember(x,InputParams.HFOBand_ToPlot),strInputBandNames,'UniformOutput',0)));
    end
else
    % Default is raw signal, signal in ripple band and signal in FR band
    PlotParams.HFOBand_ToPlot = [1,2];
end

% Signal bands to plot
if(~isfield(PlotParams,'SignalBand_ToPlot'))
    PlotParams.SignalBand_ToPlot = [1,PlotParams.HFOBand_ToPlot+1];
end

% Filter coefficients
% Needed only if the filtered signal is not provided for all bands to be plotted
cond1 = isfield(InputParams,'dataFiltered'); % Filtered signal is an input
if(cond1)
    cond2 = isempty(find(cellfun(@isempty,InputParams.dataFiltered(PlotParams.HFOBand_ToPlot)),1)); % Filtered signal given for all bands to plot
    cond3 = isempty(find(~cellfun(@(x) isequal(size(x),size(PlotParams.data)),InputParams.dataFiltered(PlotParams.HFOBand_ToPlot),'UniformOutput',1),1)); % Filtered signal has the same size as raw signal
else
    cond2 = 0;
    cond3 = 0;
end
cond4 = PlotParams.SignalBand_ToPlot==1;
condLoadCoefficients = ~((cond1&&cond2&&cond3)||cond4);
if(condLoadCoefficients)
    % File name to load filter coefficients
    if(isfield(InputParams,'FilterCoeffFilePath'))
        PlotParams.FilterCoeffFilePath = InputParams.FilterCoeffFilePath;
    else
        PlotParams.FilterCoeffFilePath = 'FIR_2kHz';
    end
    if(isfield(InputParams,'FilterCoeff'))
        if(isempty(InputParams.FilterCoeff)) % Load default if empty
            FilterCoeff = load(PlotParams.FilterCoeffFilePath);
            FilterCoeff = FilterCoeff.filter;
        else
            FilterCoeff = InputParams.FilterCoeff;
        end
    else % Load default if empty
        FilterCoeff = load('FIR_2kHz');
        FilterCoeff = FilterCoeff.filter;
    end
    % Check loaded coefficients
    if(isempty(find(~isfield(FilterCoeff,{'Rb','Ra','FRb','FRa'}),1)))
        PlotParams.FilterCoeff{1}.b = FilterCoeff.Rb;
        PlotParams.FilterCoeff{1}.a = FilterCoeff.Ra;
        PlotParams.FilterCoeff{2}.b = FilterCoeff.FRb;
        PlotParams.FilterCoeff{2}.a = FilterCoeff.FRa;
    else
        % Check whether it is in the cell format
        try
            cond1 = isempty(find(~cellfun(@(x) isempty(find(~x,1)),cellfun(@(x) isfield(x,{'a','b'}),FilterCoeff(PlotParams.HFOBand_ToPlot),'UniformOutput',0),'UniformOutput',1),1));
        catch
            error('Provided filter does not have the right format')
        end
        if(~cond1)
            error('Provided filter does not have the right format')
        end
        PlotParams.FilterCoeff = FilterCoeff;
    end
else
    FilterCoeff = [];
end

% Filtered signal
condFilterSignal = zeros(max(PlotParams.HFOBand_ToPlot),1);
if(isfield(InputParams,'dataFiltered')) % Filtered signal is an input
    if(isempty(find(ismember(find(cellfun(@isempty,InputParams.dataFiltered)),PlotParams.HFOBand_ToPlot),1))) % Filtered signal given for all bands to plot
        for nHFOBand = PlotParams.HFOBand_ToPlot
            PlotParams.dataFiltered{nHFOBand} = InputParams.dataFiltered{nHFOBand};
        end
    else % Filtered signal not given for all bands to plot
        for nHFOBand = PlotParams.HFOBand_ToPlot
            if(nHFOBand<=length(InputParams.dataFiltered))
                if(~isempty(InputParams.dataFiltered{nHFOBand})) % Filtered signal given for this band
                    % Check size of filtered signal
                    if(~isequal(size(InputParams.dataFiltered{nHFOBand}),size(PlotParams.data)))
                        error('Sizes of input raw and filtered signals do not match!')
                    end
                    PlotParams.dataFiltered{nHFOBand} = InputParams.dataFiltered{nHFOBand};
                else % Filtered signal is not an input, filter using default or input coefficients
                    condFilterSignal(nHFOBand) = 1;
                end
            else
                condFilterSignal(nHFOBand) = 1;
            end
        end
    end
else % Filter signal
    for nHFOBand = PlotParams.HFOBand_ToPlot
        condFilterSignal(nHFOBand) = 1;
    end
end
% Filter in bands that are not given
for nHFOBand = PlotParams.SignalBand_ToPlot(PlotParams.SignalBand_ToPlot-1>0)
    if(condFilterSignal(nHFOBand))
        b = PlotParams.FilterCoeff{nHFOBand}.b;
        a = PlotParams.FilterCoeff{nHFOBand}.a;
        PlotParams.dataFiltered{nHFOBand} = filtfilt(b,a,PlotParams.data')';
    end
end

% Shift in y-axis for plots
% Default values
yShiftRaw = 1000;
yShiftFiltered = [20,10];
if(isfield(InputParams,'YShift')) % Shift in y-axis given for raw signal and filtered signal
    if(length(InputParams.YShift)>=(max(PlotParams.HFOBand_ToPlot)+1)) % Sufficiently many values provided for yShift
        PlotParams.YShift = InputParams.YShift;
    elseif(length(InputParams.YShift)==(length(PlotParams.HFOBand_ToPlot)+1)) % Sufficienlty many values provided for yShift
        PlotParams.YShift(1) = InputParams.YShift(1);
        for nHFOBand = PlotParams.HFOBand_ToPlot
            PlotParams.YShift(nHFOBand+1) = InputParams.YShift(find(PlotParams.HFOBand_ToPlot==nHFOBand)+1);
        end
        % Missing values for yShift
        for nHFOBand = setdiff(1:length(yShiftFiltered),PlotParams.HFOBand_ToPlot)
            PlotParams.YShift(nHFOBand+1) = yShiftFiltered(nHFOBand);
        end
    elseif(length(InputParams.YShift)<(length(PlotParams.HFOBand_ToPlot)+1)) % Sufficiently many values not provided
        if(isempty(InputParams.YShift)) % No value provided, take the defaults
            InputParams.YShift = [yShiftRaw,yShiftFiltered];
        elseif(length(InputParams.YShift)==1) % Only value for raw signal provided
            PlotParams.YShift(1) = InputParams.YShift(1);
            PlotParams.YShift(2:3) = yShiftFiltered;
        else % Values for some filtered signals provided %% TODO
            PlotParams.YShift(1) = InputParams.YShift(1);
            PlotParams.YShift(PlotParams.SignalBand_ToPlot(1:(length(InputParams.YShift)-1))) = InputParams.YShift(2:end);
            PlotParams.YShift(PlotParams.SignalBand_ToPlot(end)) = yShiftFiltered(end);
        end
    end
else % Shift in y-axis not given
    PlotParams.YShift = [yShiftRaw,yShiftFiltered];
end

% Time window for plots
% Default values
tWindowRaw = 1.6;
tWindowFiltered = [0.3,0.3];
if(isfield(InputParams,'tWindow')) % Shift in y-axis given for raw signal and filtered signal
    if(length(InputParams.tWindow)>=(max(PlotParams.HFOBand_ToPlot)+1)) % Sufficiently many values provided for tWindow
        PlotParams.tWindow = InputParams.tWindow;
    elseif(length(InputParams.tWindow)==(length(PlotParams.HFOBand_ToPlot)+1)) % Sufficienlty many values provided for tWindow
        PlotParams.tWindow(1) = InputParams.tWindow(1);
        for nHFOBand = PlotParams.HFOBand_ToPlot
            PlotParams.tWindow(nHFOBand+1) = InputParams.tWindow(find(PlotParams.HFOBand_ToPlot==nHFOBand)+1);
        end
        % Missing values for tWindow
        for nHFOBand = setdiff(1:length(tWindowFiltered),PlotParams.HFOBand_ToPlot)
            PlotParams.tWindow(nHFOBand+1) = tWindowFiltered(nHFOBand);
        end
    elseif(length(InputParams.tWindow)<(length(PlotParams.HFOBand_ToPlot)+1)) % Sufficiently many values not provided
        if(isempty(InputParams.tWindow)) % No value provided, take the defaults
            InputParams.tWindow = [tWindowRaw,tWindowFiltered];
        elseif(length(InputParams.tWindow)==1) % Only value for raw signal provided
            PlotParams.tWindow(1) = InputParams.tWindow(1);
            PlotParams.tWindow(2:3) = tWindowFiltered;
        else % Values for some filtered signals provided %% TODO
            PlotParams.tWindow(1) = InputParams.tWindow(1);
            PlotParams.tWindow(PlotParams.SignalBand_ToPlot(1:(length(InputParams.tWindow)-1))) = InputParams.tWindow(2:end);
            PlotParams.tWindow(PlotParams.SignalBand_ToPlot(end)) = tWindowFiltered(end);
        end
    end
else % Shift in y-axis not given
    PlotParams.tWindow = [tWindowRaw,tWindowFiltered];
end

% List of channels to plot
if(isfield(InputParams,'ListOfChannels_ToPlot')) % List of channels to plot given
    if(isnumeric(InputParams.ListOfChannels_ToPlot)) % Numbers of channels given
        PlotParams.ListOfChannels_ToPlot = InputParams.ListOfChannels_ToPlot;
    elseif(iscell(InputParams.ListOfChannels_ToPlot)) % Names of channels given
        [~,ind1,ind2] = intersect(PlotParams.ElectrodeLabels(:),InputParams.ListOfChannels_ToPlot(:));
        [~,indSort] = sort(ind2);
        PlotParams.ListOfChannels_ToPlot = ind1(indSort)';
    else % Input not understood, use default
        warning('List of channels to plot not given in correct format')
        PlotParams.ListOfChannels_ToPlot = 1:size(PlotParams.data,1);
    end
else
    PlotParams.ListOfChannels_ToPlot = 1:size(PlotParams.data,1);
end

% Margins in the y-axes above and below
if(isfield(InputParams,'YMargin')) % List of channels to plot given
    PlotParams.YMargin = InputParams.YMargin;
else % Use the default values
    PlotParams.YMargin = [1,1];
end

% Signals in different bands together
PlotParams.dataAll{1} = PlotParams.data;
for nHFOBand = PlotParams.SignalBand_ToPlot(PlotParams.SignalBand_ToPlot-1>0)
    PlotParams.dataAll{nHFOBand+1} = PlotParams.dataFiltered{nHFOBand};
end

% Signal offset % TODO
% Default is removal of mean
strDetrendRawSignalMethod = 'RemoveMean';
if(isfield(InputParams,'DetrendRawSignalMethod'))
    switch(InputParams.DetrendRawSignalMethod)
        case 'RemoveMean'
            PlotParams.DetrendRawSignalMethod = strDetrendRawSignalMethod;
            nBand_ToPlot = 1;
            for nChannel = 1:size(PlotParams.dataAll{nBand_ToPlot},1)
                PlotParams.dataAll{nBand_ToPlot}(nChannel,:) = ...
                    PlotParams.dataAll{nBand_ToPlot}(nChannel,:)-mean(PlotParams.dataAll{nBand_ToPlot}(nChannel,:));
            end
        otherwise
            error('Unknown detrending method')
    end
else
    PlotParams.DetrendRawSignalMethod = strDetrendRawSignalMethod;
    nBand_ToPlot = 1;
    for nChannel = 1:size(PlotParams.dataAll{nBand_ToPlot},1)
        PlotParams.dataAll{nBand_ToPlot}(nChannel,:) = ...
            PlotParams.dataAll{nBand_ToPlot}(nChannel,:)-mean(PlotParams.dataAll{nBand_ToPlot}(nChannel,:));
    end
end

% Plot colors %% TODO default for now  % TODO different color for different event types
PlotParams.nPlotColors_EventType = {[0,0,1];[0,1,1];[1,0,0]}; % Default :'b','c','r'


% --- Executes on button press in pushbutton_changecolor.
function pushbutton_changecolor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_changecolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = handles.nMarkingGroup+1;
handles.nMarkingGroup = mod(handles.nMarkingGroup,length(handles.strSpikeColors));
if(handles.nMarkingGroup==0)
    handles.nMarkingGroup = length(handles.strSpikeColors);
end
handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_change_marking_color.
function pushbutton_change_marking_color_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_change_marking_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginputx(1);

axLimX = get(gca,'XLim');
axLimY = get(gca,'YLim');
cond1 = (x<axLimX(1))||(x>axLimX(2));
cond2 = (y<axLimY(1))||(y>axLimY(2));

if(cond1||cond2)
    fprintf('\nOut of loop\n')
else
    y_shiftRaw = handles.PlotParams.YShift(1);
    
    iChannelSpike = floor(-(y-y_shiftRaw/2)/y_shiftRaw)+1;
    nChannelSpike = handles.PlotParams.ListOfChannels_ToPlot(iChannelSpike);
    
    t_spike = x;
    
    [~,ind_t_axis] = min(abs(handles.PlotParams.t-t_spike));
    for ii = 1:size(handles.EventPosition,1)
        if(handles.EventPosition(ii,1)==nChannelSpike)
            t_spike_in_array = handles.EventPosition(ii,2);
            [~,ind_t_axis_in_array] = min(abs(handles.PlotParams.t-t_spike_in_array));
            ind_t_axis_range_in_array = (ind_t_axis_in_array-2000*0.25):(ind_t_axis_in_array+2000*0.25);
            cond = ~isempty(find(intersect(ind_t_axis_range_in_array,ind_t_axis),1));
            if(cond)
                if(handles.EventPosition(ii,3)~=-1)
                    handles.EventPosition(ii,3) = handles.nMarkingGroup;
                    handles.SpikePlots{ii}.Color = handles.strSpikeColors{handles.nMarkingGroup};
                end
            end
        end
    end
end

guidata(hObject,handles);



% --- Executes on button press in pushbutton_color_1.
function pushbutton_color_1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 1;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_color_2.
function pushbutton_color_2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 2;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_color_3.
function pushbutton_color_3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 3;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_color_4.
function pushbutton_color_4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 4;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_color_5.
function pushbutton_color_5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 5;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_color_6.
function pushbutton_color_6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nMarkingGroup = 6;

handles.pushbutton_changecolor.BackgroundColor = handles.strSpikeColors{handles.nMarkingGroup};
handles.pushbutton_changecolor.String = sprintf('Group %d',handles.nMarkingGroup);

guidata(hObject,handles);
