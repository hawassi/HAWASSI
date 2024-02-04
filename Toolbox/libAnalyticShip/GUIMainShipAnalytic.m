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

function varargout = GUIMainShipAnalytic(varargin)
% GUIMAINSHIPANALYTIC MATLAB code for GUIMainShipAnalytic.fig
%      GUIMAINSHIPANALYTIC, by itself, creates a new GUIMAINSHIPANALYTIC or raises the existing
%      singleton*.
%
%      H = GUIMAINSHIPANALYTIC returns the handle to a new GUIMAINSHIPANALYTIC or the handle to
%      the existing singleton*.
%
%      GUIMAINSHIPANALYTIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIMAINSHIPANALYTIC.M with the given input arguments.
%
%      GUIMAINSHIPANALYTIC('Property','Value',...) creates a new GUIMAINSHIPANALYTIC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIMainShipAnalytic_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIMainShipAnalytic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIMainShipAnalytic

% Last Modified by GUIDE v2.5 28-Sep-2020 07:46:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIMainShipAnalytic_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIMainShipAnalytic_OutputFcn, ...
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


% --- Executes just before GUIMainShipAnalytic is made visible.
function GUIMainShipAnalytic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIMainShipAnalytic (see VARARGIN)

% Choose default command line output for GUIMainShipAnalytic

handles.output = hObject;
handles.input=varargin{1};

pathstr=handles.input.pathnow;
projectdirectory=handles.input.projectdirectory;
additemtolistbox(handles.popupmenu_ProjDIr,projectdirectory);

handles.output = hObject;

axes(handles.axes2);     % Show the LabMath.
%javaFrame = get(hObject,'JavaFrame');
if ~isdeployed
    imshow([pathstr,'\Toolbox\lib1d\Misc\logo\hawassi_black.jpg']); %
    %   javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
else
    imshow('\Toolbox\lib1d\Misc\logo\hawassi_black.jpg');
    %  javaFrame.setFigureIcon(javax.swing.ImageIcon('\Misc\logo\icon_16.png'));
end

initialization_gui(handles);
global FlagWarn
FlagWarn=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIMainShipAnalytic wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIMainShipAnalytic_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_calc.
function pushbutton_calc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global FlagWarn
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');
keep_guis; %% close all figures except gui

Curdir   =cellstr(get(handles.popupmenu_ProjDIr,'String'));
GUISAmain.proj.projdir=Curdir{get(handles.popupmenu_ProjDIr,'Value')};
GUISAmain.proj.name=get(handles.edit_projname,'String');
GUISAmain.proj.module ='Analytic - Wave Ship';
GUISAmain.proj.projhist=handles.input.projecthist;
GUISAmain.proj.workdir=[GUISAmain.proj.projdir,'\',GUISAmain.proj.name,'\'];

GUISAmain.ship.shape=get(handles.popupmenu_ship_shape,'value');
GUISAmain.ship.shapedat=get(handles.popupmenu_ship_shape,'userdata');
GUISAmain.ship.draft=get(handles.edit_draft,'Userdata');
GUISAmain.ship.width=get(handles.edit_width,'Userdata');

GUISAmain.ship.cog=get(handles.edit_cog,'Userdata');
GUISAmain.ship.Nx=get(handles.edit_ship_Nx,'Userdata');
GUISAmain.ship.visc_sway=get(handles.edit_visc_sway,'Userdata');
GUISAmain.ship.visc_roll=get(handles.edit_visc_roll,'Userdata');
GUISAmain.ship.GM=get(handles.edit_GM,'Userdata');
GUISAmain.ship.gyradius=get(handles.edit_gyradius,'userdata');
GUISAmain.ship.mooring.check=get(handles.checkbox_mooring,'value');
GUISAmain.ship.mooring.Tn=get(handles.edit_Tn,'userdata');

GUISAmain.wave.freqVar=get(handles.popupmenu_freq,'value');
GUISAmain.wave.freq=get(handles.edit_period,'Userdata');
GUISAmain.wave.amplVar=get(handles.popupmenu_Amplitude,'value');
GUISAmain.wave.ampl=get(handles.edit_amplitude,'Userdata');

GUISAmain.bottom.depth=get(handles.edit_depth,'Userdata');
GUISAmain.bottom.slope=get(handles.edit_slope,'Userdata');
GUISAmain.calc.Id=get(handles.popupmenu_calculation,'value');
GUISAmain.calc.evmodes=get(handles.edit_Nevmodes,'Userdata');
GUISAmain.data=get(handles.pushbutton_loaddata,'Userdata');
GUISAmain.dataname=get(handles.text24,'string');
 
GUISAmain.normVartype=get(handles.popupmenu_normvartype,'value');

if isempty(GUISAmain.proj.projdir) || strcmpi(GUISAmain.proj.projdir,'--')
    statusbarObj.setText(' Specify a project directory!')
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.popupmenu_ProjDIr);
    return;
end

if isempty(GUISAmain.proj.name)
    statusbarObj.setText('Specify a Project name!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_projname);
    return
end

if GUISAmain.ship.shape==1 || GUISAmain.ship.shape==3 
    if isempty(GUISAmain.ship.width) || length(GUISAmain.ship.width)~=1
        statusbarObj.setText('Specify a width of the ship');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_width);
        return;
    end
end

if GUISAmain.ship.shape==4
     if isempty(GUISAmain.ship.shapedat) || size((GUISAmain.ship.shapedat),2)~=2
        statusbarObj.setText('Specify ship profile data');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_width);
        return;
     end
end

if GUISAmain.ship.shape~=4
    if isempty(GUISAmain.ship.draft) || length(GUISAmain.ship.draft)~=1
        statusbarObj.setText('Specify a draft of the ship');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_draft);
        return;
    end
end
if isempty(GUISAmain.ship.cog) || length(GUISAmain.ship.cog)~=2 
    statusbarObj.setText('Specify cog of the ship (x;z)');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cog);
    return;
end

if isempty(GUISAmain.ship.Nx) || length(GUISAmain.ship.Nx)~=1 
    statusbarObj.setText('Specify a number of point in  ship hull');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_ship_Nx);
    return;
end

if GUISAmain.ship.Nx<7 
    statusbarObj.setText('Number point in ship must be larger than 7');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_ship_Nx);
    return;
end

if isempty(GUISAmain.wave.freq)
    statusbarObj.setText('Specify frequency in interval (w1:dw:w2) [rad/s]');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_period);
    return; 
end

if isempty(GUISAmain.wave.ampl)
    statusbarObj.setText('Specify an amplitude');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_amplitude);
    return;
end



if isempty(GUISAmain.bottom.depth) || length(GUISAmain.bottom.depth)~=1
    statusbarObj.setText('Specify a depth');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_depth);
    return; 
end

if isempty(GUISAmain.bottom.slope) || length(GUISAmain.bottom.slope)~=1
    statusbarObj.setText('Specify a slope');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_slope);
    return; 
end

if isempty(GUISAmain.calc.evmodes)
    statusbarObj.setText('Specify number of evanescent modes');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Nevmodes);
    return; 
end

if GUISAmain.calc.Id==3
    if isempty(GUISAmain.ship.visc_sway) || length(GUISAmain.ship.visc_sway)~=1
        statusbarObj.setText('Specify a viscous damping coef for sway');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_visc_sway);
        return;
    end
    if isempty(GUISAmain.ship.visc_roll) || length(GUISAmain.ship.visc_roll)~=1
        statusbarObj.setText('Specify a viscous damping coef for roll');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_visc_roll);
        return;
    end
     if isempty(GUISAmain.ship.GM) || length(GUISAmain.ship.GM)~=1
        statusbarObj.setText('Specify a metacenter heigh (GM)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_GM);
        return;
     end  
    
     if isempty(GUISAmain.ship.gyradius) || length( GUISAmain.ship.gyradius)~=1
        statusbarObj.setText('Specify gyradius');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_gyradius);
        return;
     end  
     
     if GUISAmain.ship.mooring.check==1 && (length(GUISAmain.ship.mooring.Tn)~=3)
         statusbarObj.setText('Specify designed natural period for mooring (Tn2,Tn3,Tn4)');
         statusbarTxt = statusbarObj.getComponent(0);
         statusbarTxt.setForeground(java.awt.Color.red);
         uicontrol(handles.edit_Tn);
         return;
     end
end

save_state(GUISAmain);

if GUISAmain.calc.Id==3
    try
        tempF=load([GUISAmain.proj.workdir,'\datCalForces.mat']);
        tempA=load([GUISAmain.proj.workdir,'\datCalAddedMass.mat']);
        tempD=load([GUISAmain.proj.workdir,'\datCalDampCoef.mat']);
        
        datdiffrad.Forces=tempF.datCalForces;
        datdiffrad.AddedMass=tempA.datCalAddedMass;
        datdiffrad.DampCoef=tempD.datCalDampCoef;
        GUISAmain.dataDiffRad=datdiffrad;
    catch
        statusbarObj.setText('Please do diff. & rad calculation first');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(hObject);
        return;
    end
else
    GUISAmain.dataDiffRad=[];
end


main(GUISAmain);

% --- Executes on selection change in popupmenu_ship_shape.
function popupmenu_ship_shape_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ship_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ship_shape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ship_shape
Id=get(hObject,'Value');
set(handles.edit_draft,'enable','on');
if Id==2
set(handles.edit_width,'enable','off','string','','userdata',[]);
elseif Id==4 
set(handles.edit_width,'enable','off','string','','userdata',[]); 
set(handles.edit_draft,'enable','off','string','','userdata',[]); 

[statusbarObj]=JavaFrame_handling();

CurDir=cellstr(get(handles.popupmenu_ProjDIr, 'string'));
pathnow=CurDir{get(handles.popupmenu_ProjDIr,'Value')};
[file_name,directory]=uigetfile([pathnow,'/','*.csv;*.txt; *.dat; *.mat; *.asc'],'Load Data file');
statusbarObj.setText('Uploading data...');

if directory~=0
    temp=load([directory,file_name]);
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    if isstruct(my_data)==0 && size(my_data,2)==2
        set(hObject,'Userdata',my_data);
        statusbarObj.setText([file_name,' has been loaded']);
    else
        statusbarObj.setText('wrong input file!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        set(hObject,'Value',1);
        set(handles.edit_width,'enable','on','string','','userdata',[]);
        set(handles.edit_draft,'enable','on','string','','userdata',[]);
    end
else
    statusbarObj.setText('');
    set(hObject,'Userdata',[]);
    set(hObject,'Value',1);
    set(handles.edit_width,'enable','on','string','','userdata',[]);
    set(handles.edit_draft,'enable','on','string','','userdata',[]);
end

else
set(handles.edit_width,'enable','on');    
end

% --- Executes during object creation, after setting all properties.
function popupmenu_ship_shape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ship_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_width_Callback(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_width as text
%        str2double(get(hObject,'String')) returns contents of edit_width as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
  
end

% --- Executes during object creation, after setting all properties.
function edit_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_draft_Callback(hObject, eventdata, handles)
% hObject    handle to edit_draft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_draft as text
%        str2double(get(hObject,'String')) returns contents of edit_draft as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_draft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_draft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_period_Callback(hObject, eventdata, handles)
% hObject    handle to edit_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_period as text
%        str2double(get(hObject,'String')) returns contents of edit_period as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if isempty(param)
    statusbarObj.setText('Specify freq interval (w1:dw:w2) [rad/s]');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amplitude as text
%        str2double(get(hObject,'String')) returns contents of edit_amplitude as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)<1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_depth as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_depth);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_slope_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slope as text
%        str2double(get(hObject,'String')) returns contents of edit_slope as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_slope);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cog_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cog as text
%        str2double(get(hObject,'String')) returns contents of edit_cog as a double
param=str2num(get(hObject,'String'));
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('specify cog (x;z)');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_cog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_calculation.
function popupmenu_calculation_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_calculation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_calculation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_calculation
Id=get(hObject,'value');
if Id==3
    set(handles.edit_visc_sway,'enable','on');
    set(handles.edit_visc_roll,'enable','on');
    set(handles.edit_GM,'enable','on');
    set(handles.checkbox_mooring,'enable','on');
    set(handles.edit_gyradius,'enable','on');
    if get(handles.checkbox_mooring,'value')==1
        set(handles.edit_Tn,'enable','on');
    else
        set(handles.edit_Tn,'enable','off');
    end
else
set(handles.edit_visc_sway,'enable','off');
set(handles.edit_visc_roll,'enable','off');
set(handles.edit_GM,'enable','off');    
set(handles.checkbox_mooring,'enable','off');
set(handles.edit_Tn,'enable','off');
 set(handles.edit_gyradius,'enable','off');
end
if Id==2
 set(handles.edit_amplitude,'enable','off');   
else
 set(handles.edit_amplitude,'enable','on');   
end
set(handles.pushbutton_loaddata,'userdata',[]);
set(handles.text24,'string','no data');

Idoutvar=get(handles.popupmenu_normvartype,'value');
if Id==2 && Idoutvar>3 || Id==3 && Idoutvar>1
set(handles.popupmenu_normvartype,'value',1)
end


% --- Executes during object creation, after setting all properties.
function popupmenu_calculation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_calculation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_visc_sway_Callback(hObject, eventdata, handles)
% hObject    handle to edit_visc_sway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_visc_sway as text
%        str2double(get(hObject,'String')) returns contents of edit_visc_sway as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1 
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_visc_sway_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_visc_sway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_visc_roll_Callback(hObject, eventdata, handles)
% hObject    handle to edit_visc_roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_visc_roll as text
%        str2double(get(hObject,'String')) returns contents of edit_visc_roll as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_visc_roll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_visc_roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_projname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_projname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_projname as text
%        str2double(get(hObject,'String')) returns contents of edit_projname as a double
Name=get(hObject,'String');
Curdir   =cellstr(get(handles.popupmenu_ProjDIr,'String'));
Path=Curdir{get(handles.popupmenu_ProjDIr,'Value')};
[statusbarObj]=JavaFrame_handling();
global FlagWarn
if exist([Path,'/Output/',Name],'dir')
    statusbarObj.setText('Warning: Project exists already, it will be overwritten!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.blue);
    uicontrol(hObject);
    FlagWarn=1;
else
    statusbarObj.setText('');
    FlagWarn=0;
end

% --- Executes during object creation, after setting all properties.
function edit_projname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_projname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_loaddata.
function pushbutton_loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[statusbarObj]=JavaFrame_handling();

CurDir=cellstr(get(handles.popupmenu_ProjDIr, 'string'));
pathnow=CurDir{get(handles.popupmenu_ProjDIr,'Value')};
[file_name,directory]=uigetfile([pathnow,'/','*.csv;*.txt; *.dat; *.mat; *.asc'],'Load Data file');
statusbarObj.setText('Uploading data...');

if directory~=0
    temp=load([directory,file_name]);
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    if isstruct(my_data)==0
        set(hObject,'Userdata',my_data); 
        statusbarObj.setText([file_name,' has been loaded']);
        set(handles.text24,'string', file_name)
    else
        statusbarObj.setText('wrong input file!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        set(handles.text24,'string', 'no data')
    end
else
    statusbarObj.setText('');
    set(handles.text24,'string', 'no data')
    set(hObject,'Userdata',[]); 
end



function edit_GM_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_GM as text
%        str2double(get(hObject,'String')) returns contents of edit_GM as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_GM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Loadproject_Callback(hObject, eventdata, handles)
% hObject    handle to Loadproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata, handles)
HAWASSI_AB_startpage;

% --- Executes on button press in pushbutton_browse_dir.
function pushbutton_browse_dir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurDir=cellstr(get(handles.popupmenu_ProjDIr, 'string'));
pathNow=CurDir{get(handles.popupmenu_ProjDIr,'Value')};
curdir=uigetdir(pathNow,'browse a directory');

if curdir==0
    set(handles.popupmenu_ProjDIr,'value',1);
else
    additemtolistbox(handles.popupmenu_ProjDIr,curdir);
end

% --- Executes on selection change in listbox_dir.
function listbox_dir_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_dir


% --- Executes during object creation, after setting all properties.
function listbox_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_ProjDIr.
function popupmenu_ProjDIr_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ProjDIr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ProjDIr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ProjDIr
contents = cellstr(get(hObject,'String'));
Dir_id=get(hObject,'Value');

newstring{1}=contents{Dir_id};
j=1;
for i=2:length(contents)
    if j==Dir_id,  j=j+1;end;
    newstring{i}=contents{j};
    j=j+1;
end

set(hObject, 'string', newstring);
set(hObject, 'value', 1);

% --- Executes during object creation, after setting all properties.
function popupmenu_ProjDIr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ProjDIr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function h = additemtolistbox(h, newitem)
% ADDITEMTOLISTBOX - add a new items to the listbox
% H = ADDITEMTOLISTBOX(H, STRING)
% H listbox handle
% STRING a new item to display

oldstring = get(h, 'string');
if isempty(oldstring)
    newstring = newitem;
elseif ~iscell(oldstring)
    newstring = {newitem oldstring};
elseif  any(strcmp(oldstring,newitem)==1)
    newstring=oldstring;
else
    newstring = {newitem oldstring{:}};
end
set(h, 'string', newstring);


function initialization_gui(handles)


reset_handles(handles);
if handles.input.flagOpenProj~=0
input_handles(handles,handles.input.GUIinput);
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('Project file loaded.');
else
 set(handles.edit_projname,'string',handles.input.projectname)
end

set(handles.pushbutton_browse_dir,'tooltipString','Browse a project directory.')
set(handles.popupmenu_ProjDIr,'tooltipString','The directory list')
set(handles.edit_projname,'tooltipString',sprintf('Specify the Project name.\n The output will be saved in the directory with this name.'))
set(handles.edit_Nevmodes,'tooltipString','Specify number of evanescent modes.')
set(handles.edit_cog,'tooltipString','Specify center of gravity (x,y)');
set(handles.edit_ship_Nx,'tooltipString','Specify number of points in ship hull for discretization');
set(handles.edit_GM,'tooltipString','Specify ship metacenter height');  
set(handles.edit_Tn,'tooltipString',sprintf('Specify designed natural period due to mooring. \n Input (Tn2,Tn3,Tn4) for sway, heave and roll,\n specify 0 if a motion not applicable, ie: Tn2,0,Tn4'));
set(handles.edit_visc_sway,'tooltipString',sprintf('Specify drag coef for sway viscous damping.\n Typical value in (5,7) for barge B=60,T=12m above depth in (50,20) m'));
set(handles.edit_visc_roll,'tooltipString',sprintf('Specify drag coef for roll viscous damping.\n Typical value in (0.1,0.2)'));
set(handles.edit_period,'tooltipString',sprintf('Specify normalized wave input in interval (val1:delta val: valN)'));
set(handles.edit_depth,'tooltipString','Specify still water depth at center of the ship')
set(handles.edit_slope,'tooltipString','Specify a slope. 0 for flat bottom');
set(handles.popupmenu_normvartype,'tooltipString','Choose a normalized output variable. See the manual for detail.');

function [statusbarObj]=JavaFrame_handling()
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
try
    jFrame    = get(handle(gcf),'JavaFrame');
    jRootPane = jFrame.fHG1Client.getWindow;    %<2014b
catch
    jFrame    = get(handle(gcf),'JavaFrame');
    jRootPane = jFrame.fHG2Client.getWindow;    %>=2014b
end
statusbarObj = com.mathworks.mwswing.MJStatusBar;
jRootPane.setStatusBar(statusbarObj);

function keep_guis
fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
for fh = fig_h
    uih = findobj( fh, 'Type', 'uicontrol' );
    if isempty( uih )
        delete( fh );
    end
end



function edit_ship_Nx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ship_Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ship_Nx as text
%        str2double(get(hObject,'String')) returns contents of edit_ship_Nx as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_ship_Nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ship_Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nevmodes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Nevmodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Nevmodes as text
%        str2double(get(hObject,'String')) returns contents of edit_Nevmodes as a double
param=str2num(get(hObject,'String'));

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',param);
end

% --- Executes during object creation, after setting all properties.
function edit_Nevmodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nevmodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_state(state)
if ~isdir(state.proj.workdir)
    fc=fix(clock);
    state.proj.createddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
    state.proj.createdtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];
    mkdir(state.proj.workdir);
else
    if ~isfield(state.proj,'createddate')
        fc=fix(clock);
        state.proj.createddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
        state.proj.createdtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];
    end
end
fc=fix(clock);
state.proj.modifieddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
state.proj.modifiedtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];

state.proj.note='';
GUIinput=state;
save([state.proj.workdir,'\abprojSA_',state.proj.name,'.mat'], 'GUIinput')

if isempty(state.proj.projhist)
    state.proj.projhist = {state.proj.name, strcat(state.proj.createddate, ' - '...
        , state.proj.createdtime), strcat(state.proj.modifieddate, ' - ', ...
        state.proj.modifiedtime),state.proj.workdir, state.proj.note};
else
    if ~isempty(state.proj.projhist{1,1})
        flagNew=0;
        for ii=1:length(state.proj.projhist(:,1))
            if strcmpi(state.proj.projhist(ii,1),state.proj.name)==1 && strcmpi(state.proj.projhist(ii,4),state.proj.workdir)==1
                temp = state.proj.projhist{ii,2};
                state.proj.projhist(ii,:) =[];
                state.proj.projhist= [{state.proj.name, temp, strcat(state.proj.modifieddate, ' - ', ...
                    state.proj.modifiedtime), state.proj.workdir, state.proj.note};state.proj.projhist];
                clear temp
                flagNew=1;
                break;
            end
        end
        if flagNew==0
            state.proj.projhist = [{state.proj.name, strcat(state.proj.createddate, ' - '...
                , state.proj.createdtime), strcat(state.proj.modifieddate, ' - ', ...
                state.proj.modifiedtime),state.proj.workdir, state.proj.note};state.proj.projhist];
        end
    else
        state.proj.projhist = {state.proj.name, strcat(state.proj.createddate, ' - '...
            , state.proj.createdtime), strcat(state.proj.modifieddate, ' - ', ...
            state.proj.modifiedtime), state.proj.workdir, state.proj.note};
    end
end

if length(state.proj.projhist(:,1))>10
    state.proj.projhist = state.proj.projhist(1:10,:);
end
project = state.proj.projhist;

userfold=getenv('LOCALAPPDATA');
projhistfile=[userfold,'\projhistAB.mat'];

save(projhistfile,'project');


% --- Executes on button press in checkbox_Amplitude.
function checkbox_Amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Amplitude



function edit_kA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kA as text
%        str2double(get(hObject,'String')) returns contents of edit_kA as a double


% --- Executes during object creation, after setting all properties.
function edit_kA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_kA.
function checkbox_kA_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_kA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_kA

function reset_handles(handles)
set(handles.edit_projname,'string','');
set(handles.popupmenu_ship_shape,'value',1);
set(handles.edit_width,'Userdata',[],'string','');
set(handles.edit_draft,'Userdata',[],'string','');
set(handles.edit_cog,'Userdata',[],'string','');
set(handles.edit_ship_Nx,'Userdata',[],'string','');
set(handles.edit_visc_sway,'Userdata',[],'string','','enable','off');
set(handles.edit_visc_roll,'Userdata',[],'string','','enable','off');
set(handles.edit_GM,'Userdata',[],'string','','enable','off');
set(handles.edit_period,'Userdata',[],'string','');
set(handles.popupmenu_freq,'value',1);
set(handles.popupmenu_Amplitude,'value',1);
set(handles.edit_amplitude,'Userdata',[],'string','','enable','on');

set(handles.edit_depth,'Userdata',[],'string','');
set(handles.edit_slope,'Userdata',[],'string','');
set(handles.popupmenu_calculation,'value',1);
set(handles.edit_Nevmodes,'Userdata',[],'string','');
set(handles.checkbox_mooring,'value',0,'enable','off');
set(handles.edit_Tn,'enable','off','String','','Userdata',[]);
set(handles.edit_gyradius,'string','','userdata',[],'enable','off');
set(handles.popupmenu_normvartype,'value',1);


function input_handles(handles,input)
set(handles.popupmenu_ship_shape,'value',input.ship.shape);
set(handles.edit_projname,'string',input.proj.name);
if input.ship.shape==2
set(handles.edit_width,'Userdata',[],'string','','enable','off');   
set(handles.edit_draft,'Userdata',input.ship.draft,'string',num2str(input.ship.draft),'enable','on');
elseif input.ship.shape==4
set(handles.edit_width,'Userdata',[],'string','','enable','off');    
set(handles.edit_draft,'Userdata',[],'string','','enable','off');
set(handles.popupmenu_ship_shape,'userdata',input.ship.shapedat);
else
set(handles.edit_width,'Userdata',input.ship.width,'string',num2str(input.ship.width),'enable','on');
set(handles.edit_draft,'Userdata',input.ship.draft,'string',num2str(input.ship.draft),'enable','on');
end

set(handles.edit_cog,'Userdata',input.ship.cog,'string',[num2str(input.ship.cog(1)),',',num2str(input.ship.cog(1))]);
set(handles.edit_ship_Nx,'Userdata',input.ship.Nx,'string',num2str(input.ship.Nx));
set(handles.popupmenu_freq,'value',input.wave.freqVar);
set(handles.edit_period,'Userdata',input.wave.freq,'string',[num2str(input.wave.freq(1)),':',num2str(input.wave.freq(2)-input.wave.freq(1)),':',num2str(input.wave.freq(end))]);
set(handles.popupmenu_Amplitude,'value',input.wave.amplVar);
set(handles.edit_amplitude,'Userdata',input.wave.ampl,'string',num2str(input.wave.ampl));
set(handles.pushbutton_loaddata,'Userdata',input.data);
set(handles.text24,'string',input.dataname);


set(handles.edit_depth,'Userdata',input.bottom.depth,'string',num2str(input.bottom.depth));
set(handles.edit_slope,'Userdata',input.bottom.slope,'string',num2str(input.bottom.slope));
set(handles.popupmenu_calculation,'value',input.calc.Id);
set(handles.edit_Nevmodes,'Userdata',input.calc.evmodes,'string',num2str(input.calc.evmodes));
if input.calc.Id==3
set(handles.edit_visc_sway,'Userdata',input.ship.visc_sway,'string',num2str(input.ship.visc_sway),'enable','on');
set(handles.edit_visc_roll,'Userdata',input.ship.visc_roll,'string',num2str(input.ship.visc_roll),'enable','on');
set(handles.edit_GM,'Userdata',input.ship.GM,'string',num2str(input.ship.GM),'enable','on');
if isfield(input.ship,'mooring')
    set(handles.checkbox_mooring,'value',input.ship.mooring.check,'enable','on');
    if input.ship.mooring.check==1
        set(handles.edit_Tn,'enable','on','String',...
            [num2str(input.ship.mooring.Tn(1)),',',num2str(input.ship.mooring.Tn(2)),',',num2str(input.ship.mooring.Tn(3))]...
             ,'Userdata',input.ship.mooring.Tn);
    else
        set(handles.edit_Tn,'enable','off','String','','Userdata',[]);
    end
end
else
set(handles.edit_visc_sway,'Userdata',input.ship.visc_sway,'string',num2str(input.ship.visc_sway),'enable','off');
set(handles.edit_visc_roll,'Userdata',input.ship.visc_roll,'string',num2str(input.ship.visc_roll),'enable','off');
set(handles.edit_GM,'Userdata',input.ship.GM,'string',num2str(input.ship.GM),'enable','off');  
if isfield(input.ship,'mooring')
set(handles.checkbox_mooring,'value',input.ship.mooring.check,'enable','off');
set(handles.edit_Tn,'enable','off','String','','Userdata',[]);
end
end

if input.calc.Id==2
 set(handles.edit_amplitude,'enable','off');   
else
 set(handles.edit_amplitude,'enable','on');   
end

if isfield(input,'normVartype')
set(handles.popupmenu_normvartype,'value',input.normVartype);
end
if isfield(input.ship,'gyradius')
set(handles.edit_gyradius,'string',num2str(input.ship.gyradius),'userdata',input.ship.gyradius);
end


% --- Executes on button press in checkbox_mooring.
function checkbox_mooring_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mooring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mooring
Id=get(hObject,'Value');
if Id==1
  set(handles.edit_Tn,'enable','on');  
else
  set(handles.edit_Tn,'enable','off');  
end


function edit_Tn_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tn as text
%        str2double(get(hObject,'String')) returns contents of edit_Tn as a double
str=str2num(get(hObject,'String'));
 [statusbarObj]=JavaFrame_handling();

if isempty(str) || length(str)~=3
    statusbarObj.setText('Specify adjusted natural period [Tn1,Tn3,Tn5] for sway heave and roll. Specify 0 if not applicable.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
end
set(hObject,'Userdata',str);

% --- Executes during object creation, after setting all properties.
function edit_Tn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tn (see GCBO)
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
delete(hObject);
delete(findall(0, 'type', 'figure'));


% --- Executes on selection change in popupmenu_freq.
function popupmenu_freq_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_freq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_freq


% --- Executes during object creation, after setting all properties.
function popupmenu_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_Amplitude.
function popupmenu_Amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Amplitude contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Amplitude


% --- Executes during object creation, after setting all properties.
function popupmenu_Amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_output_var.
function popupmenu_output_var_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_output_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_output_var contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_output_var


% --- Executes during object creation, after setting all properties.
function popupmenu_output_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_output_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_normvartype.
function popupmenu_normvartype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_normvartype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_normvartype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_normvartype
Id=get(hObject,'value');
idcal=get(handles.popupmenu_calculation,'value');
if (idcal==2 && Id>3) || ( idcal==3 && Id>1)
set(hObject,'value',1);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_normvartype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_normvartype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gyradius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gyradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gyradius as text
%        str2double(get(hObject,'String')) returns contents of edit_gyradius as a double
str=str2num(get(hObject,'String'));
 [statusbarObj]=JavaFrame_handling();

if isempty(str) || length(str)~=1
    statusbarObj.setText('Specify gyradius');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hObject);
else %if correct
    statusbarObj.setText('');    
    set(hObject,'Userdata',str);
end

% --- Executes during object creation, after setting all properties.
function edit_gyradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gyradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Manual_Callback(hObject, eventdata, handles)
% hObject    handle to Manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('\Toolbox\libAnalyticShip\Misc\Manual_AnalyticShip.pdf');