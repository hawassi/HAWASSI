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

function varargout = About_AB(varargin)
% ABOUT_AB MATLAB code for About_AB.fig
%      ABOUT_AB, by itself, creates a new ABOUT_AB or raises the existing
%      singleton*.
%
%      H = ABOUT_AB returns the handle to a new ABOUT_AB or the handle to
%      the existing singleton*.
%
%      ABOUT_AB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUT_AB.M with the given input arguments.
%
%      ABOUT_AB('Property','Value',...) creates a new ABOUT_AB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before About_AB_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to About_AB_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help About_AB

% Last Modified by GUIDE v2.5 03-Aug-2015 12:03:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @About_AB_OpeningFcn, ...
                   'gui_OutputFcn',  @About_AB_OutputFcn, ...
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


% --- Executes just before About_AB is made visible.
function About_AB_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to About_AB (see VARARGIN)

% Choose default command line output for About_AB
global IDdemoAbout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


[pathstr,~,~] = fileparts(mfilename('fullpath'));
%javaFrame = get(hObject,'JavaFrame');
try
    axes(handles.axes2);     % Show the LabMath.
    if ~isdeployed
        cd (pathstr); cd ..; pathstr=pwd;
        imshow([pathstr,'\Misc\logo\hawassi_black.jpg']); %
        %  javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
    else
        imshow('hawassi_black.jpg');
        %   javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'))
    end
catch
end

IDdemoAbout=0;
if IDdemoAbout==0
set(handles.figure1,'Name',['HAWASSI AB1: About']);
set(handles.text3,'string','v.1.1 (02-12-2015)')
else
set(handles.figure1,'Name',['HAWASSI AB1 demo: About']);
set(handles.text3,'string',' v.1.1 demo')
end
log_string{1}='HAWASSI stands for';
log_string{2}='Hamiltonian Wave-Ship-Structure Interaction';
log_string{4}='HAWASSI–AB is a spatial-spectral implementation' ;
log_string{5}='           of the Analytic Boussinesq Model (AB)';
%log_string{6}='Presently the code is for simulation of wave-structure interactions,';
%log_string{7}='coupled wave-ship interaction is foreseen in future releases.';



set(handles.listbox2,'string',{log_string{:}})
% UIWAIT makes About_AB wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = About_AB_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
