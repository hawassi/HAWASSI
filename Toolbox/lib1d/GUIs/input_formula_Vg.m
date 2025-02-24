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

function varargout = input_formula_Vg(varargin)
% INPUT_FORMULA_VG MATLAB code for input_formula_Vg.fig
%      INPUT_FORMULA_VG, by itself, creates a new INPUT_FORMULA_VG or raises the existing
%      singleton*.
%
%      H = INPUT_FORMULA_VG returns the handle to a new INPUT_FORMULA_VG or the handle to
%      the existing singleton*.
%
%      INPUT_FORMULA_VG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUT_FORMULA_VG.M with the given input arguments.
%
%      INPUT_FORMULA_VG('Property','Value',...) creates a new INPUT_FORMULA_VG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before input_formula_Vg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to input_formula_Vg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help input_formula_Vg

% Last Modified by GUIDE v2.5 30-Jul-2015 11:09:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @input_formula_Vg_OpeningFcn, ...
                   'gui_OutputFcn',  @input_formula_Vg_OutputFcn, ...
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


% --- Executes just before input_formula_Vg is made visible.
function input_formula_Vg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to input_formula_Vg (see VARARGIN)

% Choose default command line output for input_formula_Vg
handles.output = hObject;
handles.input=varargin;
set(handles.figure1,'Name',['HAWASSI AB1: Input group velocity']);

% Update handles structure
guidata(hObject, handles);

% set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
% set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')

if ~isempty(handles.input)
if ~isempty(handles.input{1})
set(handles.edit1,'string',handles.input{1}.eq)
end
end

% UIWAIT makes input_formula_Vg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = input_formula_Vg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2}=get(handles.edit1,'string');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
if isfield(handles,'input')
set(handles.input{1}.hf,'value',0);
param=get(handles.edit1,'string');
set(handles.input{1}.hf,'userdata',param);
end

delete(hObject);
