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

function varargout = input_friction_interval(varargin)
% INPUT_FRICTION_INTERVAL MATLAB code for input_friction_interval.fig
%      INPUT_FRICTION_INTERVAL, by itself, creates a new INPUT_FRICTION_INTERVAL or raises the existing
%      singleton*.
%
%      H = INPUT_FRICTION_INTERVAL returns the handle to a new INPUT_FRICTION_INTERVAL or the handle to
%      the existing singleton*.
%
%      INPUT_FRICTION_INTERVAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUT_FRICTION_INTERVAL.M with the given input arguments.
%
%      INPUT_FRICTION_INTERVAL('Property','Value',...) creates a new INPUT_FRICTION_INTERVAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before input_friction_interval_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to input_friction_interval_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help input_friction_interval

% Last Modified by GUIDE v2.5 14-Aug-2015 10:40:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @input_friction_interval_OpeningFcn, ...
    'gui_OutputFcn',  @input_friction_interval_OutputFcn, ...
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


% --- Executes just before input_friction_interval is made visible.
function input_friction_interval_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to input_friction_interval (see VARARGIN)

% Choose default command line output for input_friction_interval
handles.output = hObject;
handles.input=varargin;
set(handles.figure1,'Name',['HAWASSI AB1: Input Cf.']);


% Update handles structure
guidata(hObject, handles);

if ~isempty(handles.input)
    if ~isempty(handles.input{1}.str)
            set(handles.listbox3,'Userdata',handles.input{1}.str.interval)
            set(handles.listbox4,'Userdata',handles.input{1}.str.cf)
            set(handles.listbox3,'string',handles.input{1}.str.interval)
            set(handles.listbox4,'string',handles.input{1}.str.cf)
        else
            set(handles.listbox3,'Userdata',[])
            set(handles.listbox4,'Userdata',[])
            set(handles.listbox3,'string','')
            set(handles.listbox4,'string','')
    end
end

% UIWAIT makes input_friction_interval wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = input_friction_interval_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
str.interval=get(handles.listbox3,'Userdata');
str.cf=get(handles.listbox4,'Userdata');

varargout{2}=str;


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



function edit_Interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Interval as text
%        str2double(get(hObject,'String')) returns contents of edit_Interval as a double




% --- Executes during object creation, after setting all properties.
function edit_Interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Cf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Cf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Cf as text
%        str2double(get(hObject,'String')) returns contents of edit_Cf as a double


% --- Executes during object creation, after setting all properties.
function edit_Cf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Cf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_inputdata.
function pushbutton_inputdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_inputdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Clear.
function pushbutton_Clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [statusbarObj]=JavaFrame_handling
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
jFrame    = get(handle(gcf),'JavaFrame');
try
    jRootPane = jFrame.fHG1Client.getWindow;    %<2014b
catch
    jRootPane = jFrame.fHG2Client.getWindow;    %>=2014b
end
statusbarObj = com.mathworks.mwswing.MJStatusBar;
jRootPane.setStatusBar(statusbarObj);


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
param=str2num(get(hObject,'String'));
set(handles.edit5,'userdata',param);
[statusbarObj]=JavaFrame_handling();

if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit5);
else %if correct
    statusbarObj.setText('');
    if  param(1)>param(2)
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit5);
    else %if correct
        statusbarObj.setText('');
    end
end





% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
param=str2num(get(hObject,'String'));
set(handles.edit6,'userdata',param);
[statusbarObj]=JavaFrame_handling();

if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit6);
else %if correct
    statusbarObj.setText('');
    
    if  param< 0
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit6);
    else %if correct
        statusbarObj.setText('');
    end
end




% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[statusbarObj]=JavaFrame_handling();

Interv   =get(handles.edit5,'userdata');
Cfric       =get(handles.edit6,'userdata');

if isempty(Interv)
    statusbarObj.setText('Input spatial interval');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit5);
    return;
end

if length(Interv)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit5);
    return;
else %if correct
    if  Interv(1)>Interv(2)
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit5);
        return;
    end
end

if isempty(Cfric)
    statusbarObj.setText('Input a friction coefficient');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit6);
    return;
end

if length(Cfric)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit6);
    return;
else
    if  Cfric< 0
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit6);
        return;
    end
end

intervSaved=get(handles.listbox3,'userdata');
CfricSaved=get(handles.listbox4,'userdata');
Idbreak=0;
if isempty(intervSaved)
  iter=1;
else
  Nprev=length(intervSaved);
  iter=Nprev+1;
  for i=1:Nprev
  IntervStrNow{i}=intervSaved{i};
  CfricStrNow{i}=CfricSaved{i};
  IntervSaved=str2num(intervSaved{i});
  if Interv(1)>=IntervSaved(1)&&Interv(2)<=IntervSaved(2)
      Idbreak=1;
      break;
  elseif Interv(1)<IntervSaved(1)&&...
          (Interv(2)>=IntervSaved(1)&&Interv(2)<=IntervSaved(2))
      Idbreak=1;
      break;
  elseif Interv(2)>IntervSaved(2)&&...
          (Interv(1)>=IntervSaved(1)&&Interv(1)<=IntervSaved(2))
      Idbreak=1;
      break;    
  end
  end
end
if Idbreak==1
    statusbarObj.setText('The input interval is overlapping with the existing.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit5);
 return;
end
 IntervStrNow{iter}=[num2str(Interv(1)),';',num2str(Interv(2))];
 CfricStrNow{iter}=num2str(Cfric);

set(handles.listbox3,'string',IntervStrNow);
set(handles.listbox3,'userdata',IntervStrNow);
set(handles.listbox4,'string',CfricStrNow);
set(handles.listbox4,'userdata',CfricStrNow);

set(handles.edit5,'userdata',[]);
set(handles.edit5,'string','');
set(handles.edit6,'userdata',[]);
set(handles.edit6,'string','');


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox3,'string','');
set(handles.listbox3,'userdata',[]);
set(handles.listbox4,'string','');
set(handles.listbox4,'userdata',[]);

set(handles.edit5,'userdata',[]);
set(handles.edit5,'string','');
set(handles.edit6,'userdata',[]);
set(handles.edit6,'string','');

[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isfield(handles,'input')
set(handles.input{1}.hf,'value',0);
str.interval=get(handles.listbox3,'Userdata');
str.cf=get(handles.listbox4,'Userdata');
set(handles.input{1}.hf,'userdata',str);
end

delete(hObject);
