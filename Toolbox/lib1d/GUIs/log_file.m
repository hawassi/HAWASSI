% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

function varargout = log_file(varargin)
% LOG_FILE MATLAB code for log_file.fig
%      LOG_FILE, by itself, creates a new LOG_FILE or raises the existing
%      singleton*.
%
%      H = LOG_FILE returns the handle to a new LOG_FILE or the handle to
%      the existing singleton*.
%
%      LOG_FILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOG_FILE.M with the given input arguments.
%
%      LOG_FILE('Property','Value',...) creates a new LOG_FILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before log_file_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to log_file_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help log_file

% Last Modified by GUIDE v2.5 07-Jan-2015 15:23:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @log_file_OpeningFcn, ...
                   'gui_OutputFcn',  @log_file_OutputFcn, ...
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


% --- Executes just before log_file is made visible.
function log_file_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to log_file (see VARARGIN)

% Choose default command line output for log_file
global IDdemoLog
handles.output = hObject;
handles.input=varargin;
% Update handles structure
guidata(hObject, handles);
set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')

%set(hObject,'units','normalized','Position',[0.1 0.3 0.175 0.52])
IDdemoLog=0;
if IDdemoLog==0
set(handles.figure1,'Name',['HAWASSI AB1: Log-file']);
else
set(handles.figure1,'Name',['HAWASSI AB1 demo: Log-file']);   
end
[pathstr,~,~] = fileparts(mfilename('fullpath'));
%javaFrame = get(hObject,'JavaFrame');
if ~isdeployed
      cd (pathstr); cd ..; pathstr=pwd;
 %     javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
else
 %     javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'))
end

if ~isempty(handles.input)
log_string=handles.input{1}.string;
set(handles.listbox1,'string',{log_string{:}})
end
% UIWAIT makes log_file wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = log_file_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if ishandle(handles.input{2})
set(handles.input{2},'value',0);
end
delete(hObject);
