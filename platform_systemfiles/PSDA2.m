function varargout = PSDA2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PSDA2_OpeningFcn, ...
    'gui_OutputFcn',  @PSDA2_OutputFcn, ...
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

function PSDA2_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

handles.CDataOpen               = double(imread('Open_Lock.jpg'))/255;
handles.CDataClosed             = double(imread('Closed_Lock.jpg'))/255;
handles.toggle1.CData           = double(imread('Legend.jpg'))/255;
handles.runRegular.CData         = double(imread('Play.jpg'))/255;
handles.ax2Limits.CData         = double(imread('Limits.jpg'))/255;
handles.REG_DisplayOptions.CData    = double(imread('Settings.jpg'))/255;
handles.treebutton.CData        = double(imread('tree.jpg'))/255;
handles.deleteButton.CData      = handles.CDataOpen;
handles.runPCE.CData      = double(imread('PlayPCE.jpg'))/255;
handles.PCE_DisplayOptions.CData = double(imread('Settings.jpg'))/255;
plot(handles.ax1,nan,nan,'.-','tag','lambda');
xlabel(handles.ax1,'Displacement(cm)','fontsize',10)
ylabel(handles.ax1,'lambda D','fontsize',10)

%% Retrieve data from SeismicHazard
if nargin==7 % called from SeismicHazard Toolbox
    handles = mat2psda(handles,varargin{:});
end

handles.current = [];
handles.REGOptions.mode  = 1;
handles.REGOptions.param = [1 0 0 50 0  1 0 0];
handles.REGOptions.rnd   = [];
handles.PCEOptions       = [1 0 50 1 1];

guidata(hObject, handles);
% uiwait(handles.FIGpsda);

function varargout = PSDA2_OutputFcn(hObject, eventdata, handles)
varargout{1} = [];

% ----------------  FILE MENU ---------------------------------------------
function File_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>

function ImportSeismicHazard_Callback(hObject, eventdata, handles)

if isfield(handles,'defaultpath_others')
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(handles.defaultpath_others,'*.txt'),'select file');
else
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(pwd,'*.txt'),'select file');
end

if FILTERINDEX==0,return;end

% loads hazard curves
handles.defaultpath_others=pathname;
handles = wipePSDAModel(handles);
handles = mat2psda(handles,pathname,filename);
handles.FIGpsda.Name = sprintf('%s - Probabilistic Slope Displacement Analysis - PSDA',filename);
guidata(hObject,handles)

function Save_Callback(hObject, eventdata, handles)

function Exit_Callback(hObject, eventdata, handles)
close(handles.FIGpsda)

function FIGpsda_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

% ----------------  EDIT MENU ---------------------------------------------

function Edit_Callback(hObject, eventdata, handles)

function EditLogicTree_Callback(hObject, eventdata, handles)
if ~isfield(handles,'sys'), return;end

[handles.T1,handles.T2,handles.T3,handles.Ts_param,handles.ky_param]=PSDA_Logic_tree2(...
    handles.model,...
    handles.sys,...
    handles.Ts_param,...
    handles.ky_param,...
    handles.T1,...
    handles.T2,...
    handles.T3);
[handles.table.Data,handles.IJK]=main_psda(handles.T1,handles.T2,handles.T3);

guidata(hObject,handles)

% ----------------  ANALYSIS PARAMETERS MENU ---------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)

function runRegular_Callback(hObject, eventdata, handles)
if ~isfield(handles,'sys'), return;end
if ~isempty(handles.model)
    handles=runPSDA_regular(handles);
    handles.deleteButton.CData  = handles.CDataClosed;
    plot_PSDA_regular(handles)
end
guidata(hObject,handles)

function runPCE_Callback(hObject, eventdata, handles)
if ~isfield(handles,'sys'), return;end
if ~isempty(handles.modelpce)
    handles=runPSDA_pce(handles);
    handles.deleteButton.CData  = handles.CDataClosed;
    plot_PSDA_pce(handles)
end
guidata(hObject,handles)

function deleteButton_Callback(hObject, eventdata, handles)
delete(findall(handles.ax1,'type','line'));
handles.deleteButton.CData  = handles.CDataOpen;
guidata(hObject,handles)

% ----------------  POP MENUS ---------------------------------------------

function pop_source_Callback(hObject, eventdata, handles)

function pop_source_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pop_site_Callback(hObject, eventdata, handles)
plot_PSDA_regular(handles);
guidata(hObject,handles)

function pop_site_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------  PLOTTING  ---------------------------------------------

function mouseClick(hObject,eventdata,handles)

x        = getAbsCoords(handles.ax1);
d_ptr    = interp1(handles.d,(1:length(handles.d))',x,'nearest','extrap');
d        = handles.d(d_ptr);
ch       = findall(handles.ax1,'tag','line');
ch.XData = d*[1 1];
ch.YData = handles.ax1.YLim;

function [x, y] = getAbsCoords(h_ax)
crd = get(h_ax, 'CurrentPoint');
x = crd(2,1);
y = crd(2,2);

function table_CellSelectionCallback(hObject, eventdata, handles)
if size(eventdata.Indices,1)==1
    handles.current=eventdata.Indices(1);
end
guidata(hObject,handles)

function REG_DisplayOptions_Callback(hObject, eventdata, handles)
Nbranches = size(handles.table.Data,1);
str = cell(Nbranches,1);
for i=1:Nbranches
    str{i}=sprintf('Branch %g',i);
end
handles.REGOptions = DisplacementOptions(str,handles.REGOptions);
plot_PSDA_regular(handles);
guidata(hObject,handles)

function treebutton_Callback(hObject, eventdata, handles)

if ~isfield(handles,'sys')
    return
end


[handles.T1,handles.T2,handles.T3,handles.Ts_param,handles.ky_param]=PSDA_Logic_tree2(...
    handles.model,...
    handles.sys,...
    handles.Ts_param,...
    handles.ky_param,...
    handles.T1,...
    handles.T2,...
    handles.T3);

[handles.table.Data,handles.IJK]=main_psda(handles.T1,handles.T2,handles.T3);
guidata(hObject,handles)

function ax2Limits_Callback(hObject, eventdata, handles)
handles.ax1=ax2Control(handles.ax1);

function ConditionalIntegration_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'on' , hObject.Checked='off';
    case 'off', hObject.Checked='on';
end

function ClearModel_Callback(hObject, eventdata, handles)

handles=wipePSDAModel(handles);
guidata(hObject,handles)

function handles=wipePSDAModel(handles)

if isfield(handles,'sys')
    handles=rmfield(handles,{'sys','model','modelpce','opt','h'});
end

handles.current = [];
handles.REGOptions.mode  = 1;
handles.REGOptions.param = [1 0 0 50 0  1 0 0];
handles.REGOptions.rnd   = [];
handles.PCEOptions       = [1 0 50 1 1];
handles.FIGpsda.Name     = 'Probabilistic Slope Displacement Analysis - PSDA';
handles.table.Data    = cell(0,4);
handles.tablePCE.Data = cell(0,7);
delete(findall(handles.ax1,'type','line'))
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,nan,nan,'.-','tag','lambda');

function Tools_Callback(hObject, eventdata, handles)

function PCE_DisplayOptions_Callback(hObject, eventdata, handles)
if isempty(handles.modelpce)
    return
end
str = {handles.modelpce.id};
handles.PCEOptions = DisplacementPCEOptions(handles.PCEOptions,str);
plot_PSDA_pce(handles);
guidata(hObject,handles)

function toggle1_Callback(hObject, eventdata, handles)

ch = findall(handles.FIGpsda,'type','legend');
if ~isempty(ch)
    switch ch.Visible
        case 'on',  ch.Visible='off';
        case 'off', ch.Visible='on';
    end
end
guidata(hObject,handles);

function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)

switch hObject.String
    case 'Epistemic Uncertainty through Logic Tree'
        handles.table.Enable      = 'on';
        handles.runRegular.Enable = 'on';
        handles.treebutton.Enable = 'on';
        handles.REG_DisplayOptions.Enable = 'on';
        handles.tablePCE.Enable = 'off';
        handles.runPCE.Enable   = 'inactive';
        handles.PCE_DisplayOptions.Enable='inactive';
        plot_PSDA_regular(handles)
        
    case 'Polynomial Chaos Expansion Models'
        handles.table.Enable      = 'off';
        handles.runRegular.Enable = 'inactive';
        handles.treebutton.Enable = 'inactive';
        handles.REG_DisplayOptions.Enable = 'inactive';
        handles.tablePCE.Enable = 'on';
        handles.runPCE.Enable   = 'on';
        handles.PCE_DisplayOptions.Enable='on';
        plot_PSDA_pce(handles)
end

function OpenDrivingFile_Callback(hObject, eventdata, handles)

if ispc
    winopen(handles.sys.filename)
end

function Examples_Callback(hObject, eventdata, handles)

function Ex1_Callback(hObject, eventdata, handles)
delete(findall(handles.FIGpsda,'type','legend'))
delete(findall(handles.ax1,'type','line')); drawnow
handles.ax1.NextPlot        = 'add';
handles.ax1.ColorOrderIndex = 1;

% -------------------------------------------------------------------------
% loads and plots results from independent calculations
% -------------------------------------------------------------------------
load PSDA2_Validation_Ex1_2 D_vector percentile_10_MC_example percentile_50_MC_example percentile_90_MC_example Disp_MC_samples
plot(handles.ax1,D_vector,percentile_10_MC_example)
plot(handles.ax1,D_vector,percentile_50_MC_example)
plot(handles.ax1,D_vector,percentile_90_MC_example)
plot(handles.ax1,D_vector,mean(Disp_MC_samples,1))

% -------------------------------------------------------------------------
% computes calculations directly in the platform to test psda_BT2007_pce
% -------------------------------------------------------------------------
load PSDA2_Validation_Ex1_2 Sa_vector Hazard_MC_samples PC_Coefficients_Sum_Hazards
param.sampling    = 'stable';
param.avgHazard   = 'N';
param.integration = 'mcs0';
realSa            = 100;
realD             = 100;
d                 = D_vector;
Ts_param          = [1.0 0.3 100];
ky_param          = [0.1 0.3 100];
im                = Sa_vector(:);
PCE               = Hazard_MC_samples';
Cz                = PC_Coefficients_Sum_Hazards;
hd                = psda_BT2007_pce(d, Ts_param, ky_param, im, PCE,param,Cz,realSa,realD);
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,d,prctile(hd,10),'.')
plot(handles.ax1,d,prctile(hd,50),'.')
plot(handles.ax1,d,prctile(hd,90),'.')
plot(handles.ax1,d,mean(hd,1),'.')

L=legend(...
    '10th Percentile MCS (Benchmark)', ...
    '50th Percentile MCS (Benchmark)', ...
    '90th Percentile MCS (Benchmark)', ...
    'Mean                MCS (Benchmark)', ...
    '10th Percentile MCS (PSDA2)', ...
    '50th Percentile MCS (PSDA2)', ...
    '90th Percentile MCS (PSDA2)',...
    'Mean                MCS (PSDA2)');
L.Location ='SouthWest';
L.Box      ='off';

function Ex2_Callback(hObject, eventdata, handles)

delete(findall(handles.FIGpsda,'type','legend'))
delete(findall(handles.ax1,'type','line')); drawnow
handles.ax1.NextPlot        = 'add';
handles.ax1.ColorOrderIndex = 1;

% loads and plots results from independent calculations
% -------------------------------------------------------------------------
load PSDA2_Validation_Ex1_2 D_vector percentile_10_PC_example ...
percentile_50_PC_example percentile_90_PC_example PC_term_0_summed_array
plot(handles.ax1,D_vector,percentile_10_PC_example)
plot(handles.ax1,D_vector,percentile_50_PC_example)
plot(handles.ax1,D_vector,percentile_90_PC_example)
plot(handles.ax1,D_vector,PC_term_0_summed_array(:,1)')

% -------------------------------------------------------------------------
% computes calculations directly in the platform to test psda_BT2007_pce
% -------------------------------------------------------------------------
load PSDA2_Validation_Ex1_2 Sa_vector Hazard_MC_samples PC_Coefficients_Sum_Hazards
param.sampling    = 'stable';
param.avgHazard   = 'N';
param.integration = 'pce';
realSa            = 100;
realD             = 100;
d                 = D_vector;
Ts_param          = [1.0 0.3 100];
ky_param          = [0.1 0.3 100];
im                = Sa_vector(:);
PCE               = Hazard_MC_samples';
Cz                = PC_Coefficients_Sum_Hazards;
hd                = psda_BT2007_pce(d, Ts_param, ky_param, im, PCE,param,Cz,realSa,realD);
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,d,prctile(hd,10),'.')
plot(handles.ax1,d,prctile(hd,50),'.')
plot(handles.ax1,d,prctile(hd,90),'.')
plot(handles.ax1,d,mean(hd,1),'.')
L=legend(...
    '10th Percentile PCE (Benchmark)', ...
    '50th Percentile PCE (Benchmark)', ...
    '90th Percentile PCE (Benchmark)', ...
    'Mean                PCE(Benchmark)', ...
    '10th Percentile PCE (PSDA2)', ...
    '50th Percentile PCE (PSDA2)', ...
    '90th Percentile PCE (PSDA2)',...
    'Mean                PCE(PSDA2)');
L.Location ='SouthWest';
L.Box      ='off';

function Problem1b_MCS_Callback(hObject, eventdata, handles)

delete(findall(handles.FIGpsda,'type','legend'))
delete(findall(handles.ax1,'type','line')); drawnow
handles.ax1.NextPlot        = 'add';
handles.ax1.ColorOrderIndex = 1;

% loads and plots results from independent calculations
% -------------------------------------------------------------------------
load Problem1 x yMCS Sa_vector Hazard_MC_samples xi_hazard xi_samples
plot(handles.ax1,x,yMCS)

% -------------------------------------------------------------------------
% computes calculations directly in the platform to test psda_BT2007_pce
% -------------------------------------------------------------------------
param.sampling    = 'custom';
param.zrnd        = xi_hazard;
param.xrnd        = xi_samples;
param.avgHazard   = 'Y';
param.integration = 'mcs0';
realSa            = 1000;
realD             = 1000;
D_vector          = logsp(1,100,20);
d                 = D_vector;
Ts_param          = [1.0 0.3 100];
ky_param          = [0.1 0.3 100];
im                = Sa_vector(:);
PCE               = Hazard_MC_samples(1:realSa,:)';
Cz                = [];
hd                = psda_BT2007_pce(d, Ts_param, ky_param, im, PCE,param,Cz,realSa,realD);
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,d,prctile(hd,10),'.')
plot(handles.ax1,d,prctile(hd,50),'.')
plot(handles.ax1,d,prctile(hd,90),'.')
plot(handles.ax1,d,mean(hd,1),'.')
L=legend(...
    '10th Percentile MCS (Benchmark)', ...
    '50th Percentile MCS (Benchmark)', ...
    '90th Percentile MCS (Benchmark)', ...
    'Mean                MCS (Benchmark)', ...
    '10th Percentile MCS (PSDA2)', ...
    '50th Percentile MCS (PSDA2)', ...
    '90th Percentile MCS (PSDA2)',...
    'Mean                MCS (PSDA2)');
L.Location ='SouthWest';
L.Box      ='off';

function Problem1b_PCE_Callback(hObject, eventdata, handles)

delete(findall(handles.FIGpsda,'type','legend'))
delete(findall(handles.ax1,'type','line')); drawnow
handles.ax1.NextPlot        = 'add';
handles.ax1.ColorOrderIndex = 1;

% loads and plots results from independent calculations
% -------------------------------------------------------------------------
load Problem1 x yPCE Sa_vector Hazard_MC_samples xi_hazard xi_samples PC_Coefficients_Sum_Hazards
plot(handles.ax1,x,yPCE)

% -------------------------------------------------------------------------
% computes calculations directly in the platform to test psda_BT2007_pce
% -------------------------------------------------------------------------
param.sampling    = 'custom';
param.zrnd        = xi_hazard;
param.xrnd        = xi_samples;
param.avgHazard   = 'Y';
param.integration = 'pce';
realSa            = 1000;
realD             = 1000;
D_vector          = logsp(1,100,20);
d                 = D_vector;
Ts_param          = [1.0 0.3 100];
ky_param          = [0.1 0.3 100];
im                = Sa_vector(:);
PCE               = Hazard_MC_samples(1:realSa,:)';
Cz                = PC_Coefficients_Sum_Hazards;
hd                = psda_BT2007_pce(d, Ts_param, ky_param, im, PCE,param,Cz,realSa,realD);
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,d,prctile(hd,10),'.')
plot(handles.ax1,d,prctile(hd,50),'.')
plot(handles.ax1,d,prctile(hd,90),'.')
plot(handles.ax1,d,mean(hd,1),'.')
L=legend(...
    '10th Percentile PCE (Benchmark)', ...
    '50th Percentile PCE (Benchmark)', ...
    '90th Percentile PCE (Benchmark)', ...
    'Mean                PCE (Benchmark)', ...
    '10th Percentile PCE (PSDA2)', ...
    '50th Percentile PCE (PSDA2)', ...
    '90th Percentile PCE (PSDA2)',...
    'Mean                PCE (PSDA2)');
L.Location ='SouthWest';
L.Box      ='off';
