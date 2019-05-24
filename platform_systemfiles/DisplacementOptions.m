function varargout = DisplacementOptions(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DisplacementOptions_OpeningFcn, ...
    'gui_OutputFcn',  @DisplacementOptions_OutputFcn, ...
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

function DisplacementOptions_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

if nargin==5
    handles.pan2_1.String=varargin{1};
    haz = varargin{2};
    
    switch haz.mode
        case 1
            handles.mode1.Value=1;
            handles.pan1_1.Enable='on';
            handles.pan1_2.Enable='on';
            handles.pan1_3.Enable='on';
            handles.pan1_4.Enable='on';
            handles.pan1_5.Enable='on';
            handles.mode2.Value=0;
            handles.pan2_1.Enable='off';
            handles.pan2_2.Enable='off';
            handles.pan2_3.Enable='off';
        case 2
            handles.mode1.Value=0;
            handles.pan1_1.Enable='off';
            handles.pan1_2.Enable='off';
            handles.pan1_3.Enable='off';
            handles.pan1_4.Enable='off';
            handles.pan1_5.Enable='off';
            handles.mode2.Value=1;
            handles.pan2_1.Enable='on';
            handles.pan2_2.Enable='on';
            handles.pan2_3.Enable='on';
    end
    
    handles.pan1_1.Value=haz.param(1);
    handles.pan1_2.Value=haz.param(2);
    handles.pan1_3.Value=haz.param(3);
    handles.pan1_4.String=sprintf('%g',haz.param(4));
    handles.pan1_5.Value=haz.param(5);
    
    handles.pan2_1.Value=haz.param(6);
    handles.pan2_2.Value=haz.param(7);
    handles.pan2_3.Value=haz.param(8);
    
end

switch handles.pan1.SelectedObject.String
    case 'Percentile'
        handles.pan1_4.Enable='on';
    otherwise
        handles.pan1_4.Enable='off';
end

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = DisplacementOptions_OutputFcn(hObject, eventdata, handles)

switch handles.mode1.Value
    case 1
        haz.mode  = 1;
        
    case 0
        haz.mode=2;
end
haz.param = [...
    handles.pan1_1.Value,...
    handles.pan1_2.Value,...
    handles.pan1_3.Value,...
    str2double(handles.pan1_4.String),...
    handles.pan1_5.Value,...
    handles.pan2_1.Value,...
    handles.pan2_2.Value,...
    handles.pan2_3.Value];

Nbranches = length(handles.pan2_1.String);
rnd = rand(Nbranches,1); rnd = rnd/sum(rnd);
haz.rnd = rnd;
varargout{1} = haz;
delete(handles.figure1)

function pan1_SelectionChangedFcn(hObject, eventdata, handles)

switch handles.pan1.SelectedObject.String
    case 'Percentile'
        handles.pan1_4.Enable='on';
    otherwise
        handles.pan1_4.Enable='off';
end

function pan1_1_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>

function pan1_2_Callback(hObject, eventdata, handles)

function pan1_3_Callback(hObject, eventdata, handles)

function pan1_4_Callback(hObject, eventdata, handles)

val = str2double(hObject.String);
val = min(max(val,0),100);
hObject.String=sprintf('%g',val);
guidata(hObject,handles)

function pan1_5_Callback(hObject, eventdata, handles)

function pan2_2_Callback(hObject, eventdata, handles)

switch hObject.Value
    case 0
    case 1,handles.pan2_3.Value=0;
end
guidata(hObject,handles)

function pan2_3_Callback(hObject, eventdata, handles)

switch hObject.Value
    case 0
    case 1,handles.pan2_2.Value=0;
end
guidata(hObject,handles)

function mode1_Callback(hObject, eventdata, handles)

switch hObject.Value
    case 0
        handles.pan1_1.Enable='off';
        handles.pan1_2.Enable='off';
        handles.pan1_3.Enable='off';
        handles.pan1_4.Enable='off';
        handles.pan1_5.Enable='off';
        handles.mode2.Value=1;
        handles.pan2_1.Enable='on';
        handles.pan2_2.Enable='on';
    case 1
        handles.pan1_1.Enable='on';
        handles.pan1_2.Enable='on';
        handles.pan1_3.Enable='on';
        handles.pan1_4.Enable='on';
        handles.pan1_5.Enable='on';
        handles.mode2.Value=0;
        handles.pan2_1.Enable='off';
        handles.pan2_2.Enable='off';
end

switch handles.pan1.SelectedObject.String
    case 'Percentile'
        handles.pan1_4.Enable='on';
    otherwise
        handles.pan1_4.Enable='off';
end

guidata(hObject,handles)

function mode2_Callback(hObject, eventdata, handles)
switch hObject.Value
    case 0
        handles.mode1.Value=1;
        handles.pan1_1.Enable='on';
        handles.pan1_2.Enable='on';
        handles.pan1_3.Enable='on';
        handles.pan1_4.Enable='on';
        handles.pan1_5.Enable='on';
        handles.pan2_1.Enable='off';
        handles.pan2_2.Enable='off';
        handles.pan2_3.Enable='off';
    case 1
        handles.mode1.Value=0;
        handles.pan1_1.Enable='off';
        handles.pan1_2.Enable='off';
        handles.pan1_3.Enable='off';
        handles.pan1_4.Enable='off';
        handles.pan1_5.Enable='off';
        handles.pan2_1.Enable='on';
        handles.pan2_2.Enable='on';
        handles.pan2_3.Enable='on';
end
guidata(hObject,handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);
