function varargout = setMRebins(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @setMRebins_OpeningFcn, ...
    'gui_OutputFcn',  @setMRebins_OutputFcn, ...
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

function setMRebins_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>


fillTable(handles,'R')
fillTable(handles,'M')
fillTable(handles,'E')

guidata(hObject, handles);
% uiwait(handles.figure1);
function varargout = setMRebins_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSD>
varargout{1} = [];

function handles=fillTable(handles,T)

switch T
    case 'R'
        rmin = str2double(handles.Rmin.String);
        rmax = str2double(handles.Rmax.String);
        dr   = str2double(handles.dr.String);
        Rbin =[(rmin:dr:rmax-dr)',(rmin:dr:rmax-dr)'+dr];
        NR   = size(Rbin,1);
        txt  = cell(NR,2);
        for i=1:NR
            txt{i,1}=sprintf('%7g   %3g',Rbin(i,1),Rbin(i,2));
            txt{i,2}=sprintf('%10g',1/2*(Rbin(i,1)+Rbin(i,2)));
        end
        handles.tableR.Data = txt;
    case 'M'
        mmin = str2double(handles.Mmin.String);
        mmax = str2double(handles.Mmax.String);
        dm   = str2double(handles.dm.String);
        Mbin =[(mmin:dm:mmax-dm)',(mmin:dm:mmax-dm)'+dm];
        NM   = size(Mbin,1);
        txt  = cell(NM,2);
        for i=1:NM
            txt{i,1}=sprintf('%7.1f   %3.1f',Mbin(i,1),Mbin(i,2));
            txt{i,2}=sprintf('%10g',1/2*(Mbin(i,1)+Mbin(i,2)));
        end
        handles.tableM.Data = txt;
    case 'E'
        Emin = str2double(handles.Emin.String);
        Emax = str2double(handles.Emax.String);
        de   = str2double(handles.de.String);
        Ebin =[[-inf,Emin];(Emin:de:Emax-de)',(Emin:de:Emax-de)'+de;[Emax,inf]];
        NE   = size(Ebin,1);
        txt  = cell(NE,2);
        for i=1:NE
            txt{i,1}=sprintf('%7.1f  %3.1f',Ebin(i,1),Ebin(i,2));
            txt{i,2}=sprintf('%10g',1/2*(Ebin(i,1)+Ebin(i,2)));
        end
        handles.tableE.Data = txt;
end

function Rmin_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
fillTable(handles,'R')

function Rmax_Callback(hObject, eventdata, handles)
fillTable(handles,'R')

function dr_Callback(hObject, eventdata, handles)
fillTable(handles,'R')

function Mmin_Callback(hObject, eventdata, handles)
fillTable(handles,'M')

function Mmax_Callback(hObject, eventdata, handles)
fillTable(handles,'M')

function dm_Callback(hObject, eventdata, handles)
fillTable(handles,'M')

function Emin_Callback(hObject, eventdata, handles)
fillTable(handles,'E')

function Emax_Callback(hObject, eventdata, handles)
fillTable(handles,'E')

function de_Callback(hObject, eventdata, handles)
fillTable(handles,'E')
