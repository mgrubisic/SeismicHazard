function varargout = PCE(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PCE_OpeningFcn, ...
    'gui_OutputFcn',  @PCE_OutputFcn, ...
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


function PCE_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.Exit_button.CData  = double(imread('exit.jpg'))/255;
handles.ax2Limits.CData    = double(imread('Limits.jpg'))/255;
handles.sys        = varargin{1};
handles.model      = varargin{2};
handles.opt        = varargin{3};
handles.h          = varargin{4};
handles.HazOptions = varargin{5};
isPCE              = varargin{6};

% ------------ REMOVES REGULAR MODELS -------------------------------------
handles.isREGULAR = find(horzcat(handles.model.isregular)==1);
handles.isPCE     = find(horzcat(handles.model.isregular)==0);
handles.model = handles.model(isPCE); 
[~,B]=setdiff(handles.sys.BRANCH(:,2),isPCE);
if ~isempty(B)
    handles.sys.BRANCH(B,:)=[];
    handles.sys.WEIGHT(B,:)=[];
    handles.sys.WEIGHT(:,4)=handles.sys.WEIGHT(:,4)/sum(handles.sys.WEIGHT(:,4));
end
% -------------------------------------------------------------------------

handles.menu_branch.String = {handles.model.id};
handles.menu_source.String = {handles.model(1).source.label};
handles.site_menu.String   = handles.h.p(:,1);
handles.IM_select.String   = IM2str(handles.opt.IM);

% setup ax2
xlabel(handles.ax2,handles.IM_select.String{1},'fontsize',10)
ylabel(handles.ax2,'Mean Rate of Exceedance','fontsize',10)
set(handles.ax2,'xscale','log','yscale','log')

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = PCE_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSD>
varargout{1} = [];
close(handles.figure1)

function IM_select_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

function site_menu_Callback(hObject, eventdata, handles)

function DisplayControl_Callback(hObject, eventdata, handles)

function RunButton_Callback(hObject, eventdata, handles)

delete(findall(handles.ax2,'tag','PCE'));
delete(findall(handles.ax2,'tag','perPCE'));
delete(findall(handles.ax2,'tag','perMCS'));
drawnow

branch_ptr = handles.menu_branch.Value;
source_ptr = handles.menu_source.Value;
site_ptr   = handles.site_menu.Value;
im_ptr     = handles.IM_select.Value;
RandType   = handles.rand_pop.String{handles.rand_pop.Value};
im      = handles.opt.im(:,im_ptr);
IM      = handles.opt.IM(im_ptr);
site    = cell2mat(handles.h.p(site_ptr,2:4));
Vs30    = handles.h.Vs30(site_ptr);
opt     = handles.opt;
source  = handles.model(branch_ptr).source(source_ptr);

gmpetype = source.gmpe.type;
if ~strcmp(gmpetype,'pce')
    m=warndlg(sprintf('GMM %s not valid for Polynomial Chaos Expansion',func2str(source.gmpe.handle)));
    uiwait(m);
    return
end

Nsim    = str2double(handles.Nsim.String);
[handles.PCE,handles.MCS,time1,time2]   = runPCE(source,IM,site,Vs30,im,opt.ellipsoid,Nsim,RandType);
delete(findall(handles.ax2,'tag','PCE'));
delete(findall(handles.ax2,'tag','perPCE'));
delete(findall(handles.ax2,'tag','perMCS'));
plot(handles.ax2,im,handles.PCE,'-','color',[1 1 1]*0.7,'tag','PCE','HandleVisibility','off')
handles.ax2.Layer='top';
handles.timePCE.String = sprintf('PCE : %3.2f s',time1);
handles.timeMCS.String = sprintf('MCS : %3.2f s',time2);
drawpercentiles(handles)

c2 = uicontextmenu;
uimenu(c2,'Label','Undock','Callback' ,{@figure2clipboard_uimenu,handles.ax2});
set(handles.ax2,'uicontextmenu',c2);

guidata(hObject,handles)

function drawpercentiles(handles)

if ~isfield(handles,'PCE')
    return
end
im_ptr     = handles.IM_select.Value;
im         = handles.opt.im(:,im_ptr);
per        = str2double(handles.percentiles.String);
per(isnan(per))=[];
Nper       = length(per);

if Nper>0
    perPCE     = zeros(length(im),Nper);
    perMCS     = zeros(length(im),Nper);
    leg        = cell(1,Nper);
    for i=1:Nper
        perPCE(:,i) = prctile(handles.PCE,per(i),2);
        perMCS(:,i) = prctile(handles.MCS,per(i),2);
        leg{i}   = sprintf('Percentile %g',per(i));
    end
else
    perPCE     = nan(length(im),1);
    perMCS     = nan(length(im),1);
end

delete(findall(handles.ax2,'tag','perPCE'));
delete(findall(handles.ax2,'tag','perMCS'));
handles.ax2.ColorOrderIndex=1;

plot(handles.ax2,im,perPCE,'linewidth',1.5,'tag','perPCE');
handles.ax2.ColorOrderIndex=1;
switch handles.displayMCS.Value
    case 0, VIS = 'off';
    case 1, VIS = 'on';
end
plot(handles.ax2,im,perMCS,'o','tag','perMCS','markerfacecolor','w','visible',VIS,'HandleVisibility','off');

delete(findall(handles.figure1,'tag','legend'))
if Nper>0
    Leg = legend(leg,'tag','legend');
    Leg.Location = 'SouthWest';
    Leg.Box      = 'off';
end

function menu_branch_Callback(hObject, eventdata, handles)

function menu_source_Callback(hObject, eventdata, handles)

function DisplayOptions_Callback(hObject, eventdata, handles)

function percentiles_Callback(hObject, eventdata, handles)
per  = str2double(handles.percentiles.String);
per(isnan(per))=[];
handles.percentiles.String=num2cell(per);
drawpercentiles(handles)

function DisplayBranchhes_Callback(hObject, eventdata, handles)

ch = findall(handles.ax2,'tag','PCE');
switch hObject.Value
    case 1, set(ch,'Visible','on');
    case 0, set(ch,'Visible','off');
end

function displayMCS_Callback(hObject, eventdata, handles)

ch = findall(handles.ax2,'tag','perMCS');
switch hObject.Value
    case 1, set(ch,'Visible','on');
    case 0, set(ch,'Visible','off');
end

function popupmenu10_CreateFcn(hObject, eventdata, handles)

function Nsim_Callback(hObject, eventdata, handles)

% --- Executes on button press in Exit_button.
function Exit_button_Callback(hObject, eventdata, handles)
close(handles.figure1)

function pushbutton8_Callback(hObject, eventdata, handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

function rand_pop_Callback(hObject, eventdata, handles)

function ax2Limits_Callback(hObject, eventdata, handles)
handles.ax2=ax2Control(handles.ax2);
