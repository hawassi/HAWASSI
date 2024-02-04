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

function varargout = input_bdyassim1d(varargin)
% INPUT_BDYASSIM1D MATLAB code for input_bdyassim1d.fig
%      INPUT_BDYASSIM1D, by itself, creates a new INPUT_BDYASSIM1D or raises the existing
%      singleton*.
%
%      H = INPUT_BDYASSIM1D returns the handle to a new INPUT_BDYASSIM1D or the handle to
%      the existing singleton*.
%
%      INPUT_BDYASSIM1D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUT_BDYASSIM1D.M with the given input arguments.
%
%      INPUT_BDYASSIM1D('Property','Value',...) creates a new INPUT_BDYASSIM1D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before input_bdyassim1d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to input_bdyassim1d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help input_bdyassim1d

% Last Modified by GUIDE v2.5 21-Mar-2019 02:48:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @input_bdyassim1d_OpeningFcn, ...
                   'gui_OutputFcn',  @input_bdyassim1d_OutputFcn, ...
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


% --- Executes just before input_bdyassim1d is made visible.
function input_bdyassim1d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to input_bdyassim1d (see VARARGIN)

% Choose default command line output for input_bdyassim1d
handles.output = hObject;
handles.input=varargin;
set(handles.figure1,'Name',['HAWASSI AB1: Input Bdy. Assim']);

% Update handles structure
guidata(hObject, handles);
if ~isempty(handles.input)
    if ~isempty(handles.input{1}.userdata)
            set(handles.edit_xinterv,'Userdata',handles.input{1}.userdata.xinterv)
            set(handles.edit_xinterv,'string',[num2str(handles.input{1}.userdata.xinterv(1)),';',num2str(handles.input{1}.userdata.xinterv(2))])
            set(handles.edit_smoothfact,'Userdata',handles.input{1}.userdata.smoothfact)
            set(handles.edit_smoothfact,'string',num2str(handles.input{1}.userdata.smoothfact))
            set(handles.edit_tinterv,'Userdata',handles.input{1}.userdata.tinterv)
            set(handles.edit_tinterv,'string',[num2str(handles.input{1}.userdata.tinterv(1)),';',num2str(handles.input{1}.userdata.tinterv(2))])
            set(handles.edit3,'Userdata',handles.input{1}.userdata.tstep)
            set(handles.edit3,'string',num2str(handles.input{1}.userdata.tstep))
            set(handles.pushbutton_load_eta,'Userdata',handles.input{1}.userdata.eta)
            set(handles.pushbutton_load_eta,'string','loaded');
            set(handles.popupmenu_direction,'value',handles.input{1}.userdata.direction)
            if handles.input{1}.userdata.checkbox_velocity==1
            set(handles.checkbox_potential,'value',1)
            set(handles.pushbutton2,'enable','on');
            set(handles.pushbutton2,'userdata',handles.input{1}.userdata.velocity);
            set(handles.pushbutton2,'string','loaded');
            else
            set(handles.checkbox_potential,'value',0)
            set(handles.pushbutton2,'enable','off');
            set(handles.pushbutton2,'string','load');
            set(handles.pushbutton2,'userdata',[]);
            end
            if handles.input{1}.userdata.checkbox_nonlin_adj==1
            set(handles.checkbox_nonlin_adj,'value',1);
            set(handles.text_nonlinadj_length,'enable','on');
            set(handles.edit_nonlinadj_length,'enable','on');
            set(handles.edit_nonlinadj_length,'userdata',handles.input{1}.userdata.nonlinadj_length);
            set(handles.edit_nonlinadj_length,'string',num2str(handles.input{1}.userdata.nonlinadj_length));
            set(handles.text7,'enable','on');
            set(handles.text9,'enable','on');
            else
            set(handles.checkbox_nonlin_adj,'value',0);
            set(handles.text_nonlinadj_length,'enable','off');
            set(handles.edit_nonlinadj_length,'enable','of');
            set(handles.edit_nonlinadj_length,'userdata',[]);
            set(handles.edit_nonlinadj_length,'string','');
            set(handles.text7,'enable','off');
            set(handles.text9,'enable','off');    
            end
            set(handles.pushbutton_OK,'userdata',handles.input{1}.userdata)
            
    else
            set(handles.edit_xinterv,'Userdata',[])
            set(handles.edit_xinterv,'string','')
            set(handles.edit_smoothfact,'Userdata',[])
            set(handles.edit_smoothfact,'string','')
            set(handles.edit_tinterv,'Userdata',[])
            set(handles.edit_tinterv,'string','')
            set(handles.edit3,'Userdata',[])
            set(handles.edit3,'string','')
            set(handles.pushbutton_load_eta,'Userdata',[])
            set(handles.pushbutton_load_eta,'string','load');
            set(handles.checkbox_potential,'value',0)
            set(handles.pushbutton2,'enable','off');
            set(handles.pushbutton2,'string','load');
            set(handles.pushbutton2,'userdata',[]);  
            set(handles.checkbox_nonlin_adj,'value',0);
            set(handles.text_nonlinadj_length,'enable','off');
            set(handles.edit_nonlinadj_length,'enable','of');
            set(handles.edit_nonlinadj_length,'userdata',[]);
            set(handles.edit_nonlinadj_length,'string','');
            set(handles.text7,'enable','off');
            set(handles.text9,'enable','off');    
    end
else
           set(handles.edit_xinterv,'Userdata',[])
            set(handles.edit_xinterv,'string','')
            set(handles.edit_smoothfact,'Userdata',[])
            set(handles.edit_smoothfact,'string','')
            set(handles.edit_tinterv,'Userdata',[])
            set(handles.edit_tinterv,'string','')
            set(handles.edit3,'Userdata',[])
            set(handles.edit3,'string','')
            set(handles.pushbutton_load_eta,'Userdata',[])
            set(handles.pushbutton_load_eta,'string','load');
            set(handles.checkbox_potential,'value',0)
            set(handles.pushbutton2,'enable','off');
            set(handles.pushbutton2,'string','load');
            set(handles.pushbutton2,'userdata',[]);
            set(handles.checkbox_nonlin_adj,'value',0);
            set(handles.text_nonlinadj_length,'enable','off');
            set(handles.edit_nonlinadj_length,'enable','of');
            set(handles.edit_nonlinadj_length,'userdata',[]);
            set(handles.edit_nonlinadj_length,'string','');
            set(handles.text7,'enable','off');
            set(handles.text9,'enable','off'); 
end

% UIWAIT makes input_bdyassim1d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = input_bdyassim1d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

varargout{2}=get(handles.pushbutton_OK,'userdata');


function edit_xinterv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xinterv as text
%        str2double(get(hObject,'String')) returns contents of edit_xinterv as a double
param=str2num(get(hObject,'String'));
set(handles.edit_xinterv,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xinterv);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_xinterv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tinterv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tinterv as text
%        str2double(get(hObject,'String')) returns contents of edit_tinterv as a double
param=str2num(get(hObject,'String'));
set(handles.edit_tinterv,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_tinterv);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_tinterv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
param=str2num(get(hObject,'String'));
set(handles.edit3,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit3);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_load_eta.
function pushbutton_load_eta_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');
[file_name,directory]=uigetfile([handles.input{1}.projdir,'\','*.txt; *.dat; *.mat; *.asc'],'Load assimilation data');
if directory~=0
    temp=load([directory,file_name]);
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    FlagDat=1;
    if ~isstruct(my_data)
%         x=my_data(1,2:end);t=my_data(2:end,1);
%         if length(my_data(4,:))~=Nx
%             FlagDat=0;
%         end
%         if length(my_data(4:end,:))~=Ny*Nt
%             FlagDat=0;
%         end
    else
        FlagDat=0;
    end
    
    if FlagDat==1
        x=my_data(1,2:end);t=my_data(2:end,1);
        set(hObject,'userdata',my_data);
        set(handles.edit_xinterv,'userdata',[x(1) x(end)],...
            'string',[num2str(x(1)),';',num2str(x(end))]);
        smmothfact=round((x(end)-x(1))*0.1/(x(2)-x(1)));
        set(handles.edit_smoothfact,'userdata',smmothfact,...
            'string',num2str(smmothfact));
        set(handles.edit_tinterv,'userdata',[t(1) t(end)],...
            'string',[num2str(t(1)),';',num2str(t(end))]);
        dt=t(2)-t(1);
        set(handles.edit3,'userdata',dt,...
            'string',num2str(dt));
        statusbarObj.setText('data loaded');
        set(hObject,'string','loaded');
    else
        set(hObject,'string','load');
        statusbarObj.setText('data is not loaded');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(hObject);
    end
else
    set(hObject,'string','load');
    statusbarObj.setText('data is not loaded');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');
[file_name,directory]=uigetfile([handles.input{1}.projdir,'\','*.txt; *.dat; *.mat; *.asc'],'Load assimilation data');
if directory~=0
    temp=load([directory,file_name]);
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    FlagDat=1;
    if ~isstruct(my_data)
%         x=my_data(1,2:end);t=my_data(2:end,1);
%         if length(my_data(4,:))~=Nx
%             FlagDat=0;
%         end
%         if length(my_data(4:end,:))~=Ny*Nt
%             FlagDat=0;
%         end
    else
        FlagDat=0;
    end
    
    if FlagDat==1
        set(hObject,'userdata',my_data);
        statusbarObj.setText('data loaded');
        set(hObject,'string','loaded');
    else
        set(hObject,'string','load');
        statusbarObj.setText('data is not loaded');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(hObject);
    end
else
    set(hObject,'string','load');
    statusbarObj.setText('data is not loaded');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
end


% --- Executes on button press in checkbox_potential.
function checkbox_potential_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_potential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_potential
Id=get(handles.checkbox_potential,'value');
if Id==1
            set(handles.pushbutton2,'enable','on');
else
            set(handles.pushbutton2,'enable','off');
end

% --- Executes on button press in checkbox_nonlin_adj.
function checkbox_nonlin_adj_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nonlin_adj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nonlin_adj
 Id=get(handles.checkbox_nonlin_adj,'value');
 if Id==1
 set(handles.text_nonlinadj_length,'enable','on');
 set(handles.edit_nonlinadj_length,'enable','on');
 set(handles.text7,'enable','on');
 set(handles.text9,'enable','on');
 else
 set(handles.text_nonlinadj_length,'enable','off');
 set(handles.edit_nonlinadj_length,'enable','off');
 set(handles.text7,'enable','off');
 set(handles.text9,'enable','off');    
 end

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata.xinterv=get(handles.edit_xinterv,'userdata');
userdata.smoothfact=get(handles.edit_smoothfact,'userdata');
userdata.tinterv=get(handles.edit_tinterv,'userdata');
userdata.tstep=get(handles.edit3,'userdata');
userdata.eta=get(handles.pushbutton_load_eta,'userdata');
userdata.direction=get(handles.popupmenu_direction,'value');
userdata.checkbox_velocity=get(handles.checkbox_potential,'value');
userdata.velocity=get(handles.pushbutton2,'userdata');
userdata.checkbox_nonlin_adj=get(handles.checkbox_nonlin_adj,'value');
userdata.nonlinadj_length=get(handles.edit_nonlinadj_length,'userdata');

[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');



if isempty(userdata.eta)
  statusbarObj.setText('Specify elevation data!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.pushbutton_load_eta);
       return;  
end

if userdata.checkbox_velocity==1
   if isempty(userdata.velocity)
       statusbarObj.setText('Specify velocity data!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       set(handles.pushbutton2,'string','load')
       uicontrol(handles.pushbutton2);
       return;
   end 
end

if isempty(userdata.xinterv)|| length(userdata.xinterv)~=2
  statusbarObj.setText('Specify x interval!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.edit_xinterv);
       return;  
end

if isempty(userdata.smoothfact)|| length(userdata.smoothfact)~=1
       statusbarObj.setText('Specify smooth factor!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.edit_smoothfact);
       return;  
end

if isempty(userdata.tinterv)|| length(userdata.tinterv)~=2
  statusbarObj.setText('Specify t interval!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.edit_tinterv);
       return;  
end

if isempty(userdata.tstep) || length(userdata.tstep)~=1
  statusbarObj.setText('Specify time step!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.edit3);
       return;  
end

if userdata.checkbox_nonlin_adj==1
    if isempty(userdata.nonlinadj_length)|| length(userdata.nonlinadj_length)~=1
       statusbarObj.setText('Specify length of nonlinear adjustment!')
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red);
       uicontrol(handles.edit_nonlinadj_length);
       return;
   end 
end
set(hObject,'userdata',userdata);
set(handles.input{1}.hf,'userdata',userdata);
figure1_CloseRequestFcn(handles.figure1, eventdata, handles)


function edit_nonlinadj_length_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nonlinadj_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nonlinadj_length as text
%        str2double(get(hObject,'String')) returns contents of edit_nonlinadj_length as a double
param=str2num(get(hObject,'String'));
set(hObject,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_nonlinadj_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nonlinadj_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_direction.
function popupmenu_direction_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_direction


% --- Executes during object creation, after setting all properties.
function popupmenu_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_smoothfact_Callback(hObject, eventdata, handles)
% hObject    handle to edit_smoothfact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_smoothfact as text
%        str2double(get(hObject,'String')) returns contents of edit_smoothfact as a double
param=str2num(get(hObject,'String'));
set(hObject,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_smoothfact_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smoothfact (see GCBO)
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
userdata=get(handles.pushbutton_OK,'Userdata');
set(handles.input{1}.hf,'userdata',userdata);
end


delete(hObject);
