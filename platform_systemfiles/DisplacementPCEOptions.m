function varargout = DisplacementPCEOptions(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DisplacementPCEOptions_OpeningFcn, ...
    'gui_OutputFcn',  @DisplacementPCEOptions_OutputFcn, ...
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

function DisplacementPCEOptions_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

if nargin==5
    param = varargin{1};
    handles.pan1_1.Value=param(1);
    handles.pan1_2.Value=param(2);
    handles.pan1_3.String=sprintf('%g',param(3));
    handles.pan1_4.Value=param(4);
    
    handles.pop1.String=varargin{2};
    handles.pop1.Value=param(5);
    
end

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = DisplacementPCEOptions_OutputFcn(hObject, eventdata, handles)

param = [...
    handles.pan1_1.Value,...
    handles.pan1_2.Value,...
    str2double(handles.pan1_3.String),...
    handles.pan1_4.Value,...
    handles.pop1.Value];
varargout{1} = param;
delete(handles.figure1)

function pan1_SelectionChangedFcn(hObject, eventdata, handles)

function pan1_1_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>

function pan1_2_Callback(hObject, eventdata, handles)

function pan1_3_Callback(hObject, eventdata, handles)
val = str2double(hObject.String);
val = min(max(val,0),100);
hObject.String=sprintf('%g',val);
guidata(hObject,handles)

function pan1_4_Callback(hObject, eventdata, handles)

guidata(hObject,handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

function pop1_Callback(hObject, eventdata, handles)
