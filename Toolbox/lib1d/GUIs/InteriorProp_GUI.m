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

function varargout = InteriorProp_GUI(varargin)
% INTERIORPROP_GUI MATLAB code for InteriorProp_GUI.fig
%      INTERIORPROP_GUI, by itself, creates a new INTERIORPROP_GUI or raises the existing
%      singleton*.
%
%      H = INTERIORPROP_GUI returns the handle to a new INTERIORPROP_GUI or the handle to
%      the existing singleton*.
%
%      INTERIORPROP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERIORPROP_GUI.M with the given input arguments.
%
%      INTERIORPROP_GUI('Property','Value',...) creates a new INTERIORPROP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InteriorProp_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InteriorProp_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InteriorProp_GUI

% Last Modified by GUIDE v2.5 14-Jul-2015 13:08:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InteriorProp_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @InteriorProp_GUI_OutputFcn, ...
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


% --- Executes just before InteriorProp_GUI is made visible.
function InteriorProp_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InteriorProp_GUI (see VARARGIN)

% Choose default command line output for InteriorProp_GUI
reset_handles(handles,0);
handles.output = hObject;
handles.input=varargin;
% Update handles structure
guidata(hObject, handles);
%set(0,'ScreenPixelsPerInch',96)

pathstr=handles.input{1}.Proj.path;


axes(handles.logo1);
%javaFrame = get(hObject,'JavaFrame');
if ~isdeployed
    imshow([pathstr,'\Toolbox\lib1d\Misc\logo\hawassi_black.jpg']); %
  %  javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
else
    try
    imshow('hawassi_black.jpg');
  %  javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'))
    catch
    end
end

set(handles.figure1,'Name',['HAWASSI AB1: Internal Flow']);

%set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
%set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')


set(handles.edit_z1z2,'enable','off');
set(handles.edit_saveId,'visible','off');
set(handles.text44,'visible','off');


GIF_set=[1/100;inf];

set(handles.checkbox_PP2_GIF,'value',0);
set(handles.edit_PP2_GIF,'Enable','off');
set(handles.edit_PP2_GIF,'Userdata',GIF_set);
set(handles.edit_PP2_GIF,'String','1/100;inf');
%set(handles.checkbox_P,'value',1);



if isfield(handles.input{1},'dteta')
    projdir=handles.input{1}.Proj.projdir;
    additemtolistbox(handles.popupmenu_CurDir,projdir);
    set(handles.popupmenu3,'value',2);
    x=handles.input{1}.x;
    T=handles.input{1}.time;
    bathy=handles.input{1}.bathy;
    minbathy=min(handles.input{1}.bathy);

    if any(bathy>0)
    [MaxPlot, IndxMax]=max(max(handles.input{1}.eta)-bathy');
     maxeta=(1.2*MaxPlot+bathy(IndxMax));
    else
    maxeta=max(max(handles.input{1}.eta));
    end
    
    dzSuggest=(maxeta-minbathy)./50;
    PP_Z=[minbathy;maxeta];
    Fbdy=handles.input{1}.Fbdy;
    
    PP_X=[x(1)+Fbdy(1);x(end)-Fbdy(2)];
    set(handles.edit_X,'Userdata',PP_X);
    set(handles.edit_X,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
    set(handles.edit_xcoarse,'Userdata',1);
    set(handles.edit_xcoarse,'String','1');
    PP_T=[T(1);T(end)];
    set(handles.edit_T,'Userdata',PP_T);
    set(handles.edit_T,'String',[num2str(T(1)),';',num2str(T(end))]);
    set(handles.edit_tcoarse,'Userdata',1);
    set(handles.edit_tcoarse,'String','1');
    
    set(handles.edit_z1z2,'Userdata',PP_Z);
    set(handles.edit_z1z2,'String',[num2str(PP_Z(1),2),';',num2str(PP_Z(2),2)]);
    set(handles.edit_dz,'Userdata',dzSuggest);
    set(handles.edit_dz,'String',num2str(dzSuggest));
    
    set(findall(handles.uipanel1, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off')
else
    additemtolistbox(handles.popupmenu_CurDir,handles.input{1}.Proj.projdir);
    set(findall(handles.uipanel1, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off');
    set(handles.radiobutton1,'value',0);
    set(handles.radiobutton2,'value',0);
    set(handles.radiobutton3,'value',0);
    set(handles.radiobutton1,'enable','on');
    set(handles.radiobutton2,'enable','on');
    set(handles.radiobutton3,'enable','on');
    set(handles.edit_saveId_PP,'Enable','off');
    set(handles.checkbox_savedata_PP,'value',0);
    set(handles.popupmenu_saveFig,'Enable','off');
    set(handles.checkbox_savefigure,'value',0);
end

set(handles.edit_PP1_line_X,'enable','off');
set(handles.edit_PP1_line_Z,'enable','off');
set(handles.edit_PP1_line_T,'enable','off');
set(handles.edit_PP1_XT_Z,'enable','off');
set(handles.edit_PP1_ZT_X,'enable','off');
set(handles.edit_PP1_XZ_time,'enable','off');

%%%%%%%%%%%%%%%%ToolTips%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.popupmenu_CurDir,'tooltipString','The directory list')
set(handles.pushbutton_browse_Dir,'tooltipString','Specify the directory')
set(handles.popupmenu3,'tooltipString','Data status: Other: data need to be provided, please load.')
set(handles.edit_T,'tooltipString','Specify the time interval. [t_start;t_end]')
set(handles.edit_tcoarse,'tooltipString','Specify a coarse factor of the time. Input an integer [1,2,3, ..]. No coarse applied if input is 1. ')
set(handles.edit_X,'tooltipString','Specify the horizontal interval. [x_start;x_end]')
set(handles.edit_xcoarse,'tooltipString','Specify a coarse factor of the horiontal domain. Input an integer [1,2,3, ..]. No coarse applied if input is 1. ')
set(handles.edit_z1z2,'tooltipString','Specify the vertical interval. [z_bottom;z_top]')
set(handles.edit_dz,'tooltipString','Specify the vertical grid size. [dz]')
set(handles.checkbox_savedata,'tooltipString','Option to save the data of the interior properties.')
set(handles.edit_saveId,'tooltipString','Specify the ID of data. ID: 1,2,3,...')
set(handles.pushbutton_calculate_IP,'tooltipString','Start to calculate the interior properties.')

set(handles.popupmenu_Input,'tooltipString','Data status: Other: data need to be provided, please load.')
set(handles.radiobutton1,'tooltipString','Option to show the density plot and/or the line plot')
set(handles.checkbox_XZ,'tooltipString','Option to show the interior properties in the axes of the vertical spatial vs. the horizontal spatial ')
set(handles.checkbox_ZT,'tooltipString','Option to show the interior properties in the axes of the vertical spatial vs. the time')
set(handles.checkbox_XT,'tooltipString','Option to show the interior properties in the axes of the time vs. the horizontal spatial ')
set(handles.edit_PP1_XZ_time,'tooltipstring','Specify a single time. [t1]');
set(handles.edit_PP1_ZT_X,'tooltipstring','Specify a single horizontal position. [x1]');
set(handles.edit_PP1_XT_Z,'tooltipstring','Specify a single vertical position. [z1]');
set(handles.checkbox_PP1_line_X,'tooltipString','Option to show the interior properties in the horizontal spatial axes');
set(handles.checkbox_PP1_line_Z,'tooltipString','Option to show the interior properties in the vertical spatial axes');
set(handles.checkbox_PP1_line_T,'tooltipString','Option to show the interior properties in the time axes');
set(handles.edit_PP1_line_X,'tooltipString','Specify a single time and vertical position. [t1;z1]');
set(handles.edit_PP1_line_Z,'tooltipString','Specify a single time and horizontal position. [t1;x1]');
set(handles.edit_PP1_line_T,'tooltipString','Specify a single horizontal and vertical position. [x1;z1]');

set(handles.radiobutton2,'tooltipString','Option to show the the animation of the simulation')
set(handles.edit_PP2_T,'tooltipString','Specify the time interval. [t_start;t_end]');
set(handles.edit_PP2_X,'tooltipString','Specify the horizontal interval. [x_start;x_end]');
set(handles.edit_PP2_Z,'tooltipString','Specify the vertical interval. [z_bottom;z_top]');

set(handles.radiobutton3,'tooltipString','Option to show comparison between the simulation and the measurement');
set(handles.pushbutton_measdata,'tooltipString','Provide the measurement data');
set(handles.edit26,'tooltipString','This is the buoy position that loaded from the file');
set(handles.checkbox_PP3_densityplot,'tooltipString','Option to show the density plot.');
set(handles.checkbox_PP3_errorplot,'tooltipString','Option to show the error plot.');
set(handles.checkbox_PP3_lineplot,'tooltipString','Option to show the line plot.');
set(handles.checkbox39,'tooltipString','Option to show the multiple time traces at the vertical measurement positions')
set(handles.edit_dZmeas,'tooltipString','Specify a step of the measured data in the vertical direction. #step: 1,2,3, ... . If the input is 1, All data will be shown.')
set(handles.checkbox38,'tooltipString','Option to show a single time traces at a vertical measurement position');
set(handles.checkbox_timeMeanVariation,'tooltipString','Option to show the time mean variation plot.')

set(handles.checkbox25,'tooltipString','Option to adjust the x interval')
set(handles.edit_xlim,'tooltipString','Specify the horizontal interval. [x_min;x_max]')
set(handles.checkbox_setting_ylim,'tooltipString','Option to adjust the vertical interval')
set(handles.edit34,'tooltipString','Specify the vertical interval. [z_min;z_max]')
set(handles.checkbox22,'tooltipString','Option to adjust the colorbar interval')
set(handles.edit_colorbar,'tooltipString','Specify the the colorbar interval. [min;max]')
set(handles.checkbox_tcoarse,'tooltipString','Option to coarse the data in time')
set(handles.edit_PP_tcoarse,'tooltipString','Specify the coarsing factor. [1,2,...]')
set(handles.checkbox_PP2_GIF,'tooltipString','Option to set the output .gif file')
set(handles.edit_PP2_GIF,'tooltipString','Adjust the delay time and number of loop. Input format [delay_time;#loop], default: [1/100;inf]')

set(handles.checkbox_savedata_PP,'tooltipString','Option to save the data of the postprocessing.')
set(handles.edit_saveId_PP,'tooltipString','Specify the ID of data. ID: 1,2,3,...')

set(handles.pushbutton4,'tooltipString','Start to process.')
set(handles.pushbutton_closefig,'tooltipString','Closes all active figures.')
% UIWAIT makes InteriorProp_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InteriorProp_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_CurDir.
function popupmenu_CurDir_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_CurDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_CurDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_CurDir
contents = cellstr(get(hObject,'String'));
Dir_id=get(hObject,'Value');

newstring{1}=contents{Dir_id};
j=1;
for i=2:length(contents)
    if j==Dir_id,  j=j+1;end;
    newstring{i}=contents{j};
    j=j+1;
end

set(handles.popupmenu_CurDir, 'string', newstring);
set(handles.popupmenu_CurDir, 'value', 1);



% --- Executes during object creation, after setting all properties.
function popupmenu_CurDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_CurDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_browse_Dir.
function pushbutton_browse_Dir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse_Dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};
curdir=uigetdir(pathNow,'browse a directory');
if curdir==0
    set(handles.popupmenu_CurDir,'value',1);
else
    additemtolistbox(handles.popupmenu_CurDir,curdir);
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



% --- Executes on button press in pushbutton_load_input.
function pushbutton_load_input_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();

statusbarObj.setText('loading...');

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathNow,'\','*.mat'],'Load Data file');
if directory~=0
    %cd(directory);
    temp=load([directory,file_name]);
    
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    
    
    
    if isempty(my_data)|| ~isfield(my_data,'timeIP')
        set(handles.popupmenu_Input,'Value',1);
        statusbarObj.setText('No data, please load!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    else
        set(handles.popupmenu_Input,'Value',2)
        set(handles.pushbutton_load_input,'Enable','off')
        Tinterv=[my_data.timeIP(1);my_data.timeIP(end)];
        set(handles.edit_PP2_T,'Userdata',Tinterv);
        set(handles.edit_PP2_T,'String',[num2str(roundn(Tinterv(1),-2)),';',num2str(roundn(Tinterv(2),-2))]);
        Xinterv=[my_data.X(1);my_data.X(end)];
        set(handles.edit_PP2_X,'Userdata',Xinterv);
        set(handles.edit_PP2_X,'String',[num2str(roundn(Xinterv(1),-2)),';',num2str(roundn(Xinterv(2),-2))]);
        Zinterv=[my_data.Z(1);my_data.Z(end)];
        set(handles.edit_PP2_Z,'Userdata',Zinterv);
        set(handles.edit_PP2_Z,'String',[num2str(roundn(Zinterv(1),-2)),';',num2str(roundn(Zinterv(2),-2))]);
        
        set(handles.pushbutton_load_input,'Userdata',my_data);
        %cd(pathnow);
        statusbarObj.setText([file_name,' has been loaded']);
    end
else
    statusbarObj.setText('');
end


% --- Executes on selection change in popupmenu_Input.
function popupmenu_Input_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Input contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Input
contents = cellstr(get(hObject,'String'));
set(handles.popupmenu_Input,'value',1);
Cont=contents{1};

Data_Input_load =findobj('Tag','pushbutton_load_input');

if strcmp(Cont,'Other')
    set(Data_Input_load,'Enable','on');
else
    set(Data_Input_load,'Enable','off');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_P_dyn.
function checkbox_P_dyn_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_dyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_P_dyn


% --- Executes on button press in checkbox_P_NonlinDyn.
function checkbox_P_NonlinDyn_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_NonlinDyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_P_NonlinDyn


% --- Executes on button press in checkbox_ax.
function checkbox_ax_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ax


% --- Executes on button press in checkbox_az.
function checkbox_az_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_az


% --- Executes on button press in checkbox_A.
function checkbox_A_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_A


% --- Executes on button press in checkbox_P.
function checkbox_P_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_P


% --- Executes on button press in checkbox_Vx.
function checkbox_Vx_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Vx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Vx


% --- Executes on button press in checkbox_Vz.
function checkbox_Vz_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Vz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Vz


% --- Executes on button press in checkbox_V.
function checkbox_V_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_V


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');

IP_data=get(handles.pushbutton_load_input,'Userdata');
%assignin('base','IP_data',IP_data);

if isempty(IP_data)
    set(handles.popupmenu_Input,'Value',1);
    set(handles.pushbutton_load_input,'enable','on');
    statusbarObj.setText('No data, please load!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end


if ~isfield(IP_data,'Ptot')
    statusbarObj.setText('preparation..., please wait!');
    IP_data.Ptot=total_pressure(IP_data);
    statusbarObj.setText('');
    set(handles.pushbutton_load_input,'Userdata',IP_data)
else
    statusbarObj.setText('');
end


PP1=get(handles.radiobutton1,'Value');
PP2=get(handles.radiobutton2,'Value');
PP3=get(handles.radiobutton3,'Value');

if PP1+PP2+PP3==0
    statusbarObj.setText('Choose a posprocessing mode!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if PP1==1, PP=1;
elseif PP2==1, PP=2;
elseif PP3==1, PP=3;
end

if PP3==1
    if isempty(get(handles.pushbutton_measdata,'userdata'))
        statusbarObj.setText('Please load measurement data');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if get(handles.checkbox_PP3_lineplot,'value')+...
            get(handles.checkbox_PP3_densityplot,'value')+...
            get(handles.checkbox_PP3_errorplot,'value')+...
            get(handles.checkbox_timeMeanVariation,'value')==0
        
        statusbarObj.setText('Choose one plot (z-t axes plot, line plot or time mean plot)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if get(handles.checkbox_PP3_lineplot,'value')==1 && ...
            (get(handles.checkbox_PP3_densityplot,'value')+...
            get(handles.checkbox_PP3_errorplot,'value'))>0 ...
            &&  get(handles.checkbox_timeMeanVariation,'value')>0
        
        statusbarObj.setText('Choose one plot (z-t axes plot, line plot or time mean plot)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if get(handles.checkbox_PP3_lineplot,'value')==0 && ...
            (get(handles.checkbox_PP3_densityplot,'value')+...
            get(handles.checkbox_PP3_errorplot,'value'))>0 ...
            &&  get(handles.checkbox_timeMeanVariation,'value')==1
        
        statusbarObj.setText('Choose one plot (z-t axes plot, line plot or time mean plot)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    
    if get(handles.checkbox_PP3_lineplot,'value')==1 && ...
            (get(handles.checkbox_PP3_densityplot,'value')+...
            get(handles.checkbox_PP3_errorplot,'value'))>0 ...
            &&  get(handles.checkbox_timeMeanVariation,'value')==0
        
        statusbarObj.setText('Choose one plot (z-t axes plot, line plot or time mean plot)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if get(handles.checkbox_PP3_lineplot,'value')==1 && ...
            (get(handles.checkbox_PP3_densityplot,'value')+...
            get(handles.checkbox_PP3_errorplot,'value'))==0 ...
            &&  get(handles.checkbox_timeMeanVariation,'value')==1
        
        statusbarObj.setText('Choose one plot (z-t axes plot, line plot or time mean plot)');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    
    
    
    
    if get(handles.checkbox_PP3_lineplot,'value')==1
        if get(handles.checkbox38,'value')+get(handles.checkbox39,'value')==0
            statusbarObj.setText('Choose one type of Z interval');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            
        end
        
        if get(handles.checkbox38,'value')==1
            if isempty(get(handles.edit_PP3_Z,'userdata'))
                statusbarObj.setText('Please input z[m]!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                return;
            end
        elseif get(handles.checkbox39,'value')==1
            if isempty(get(handles.edit_dZmeas,'userdata'))
                statusbarObj.setText('Please input step of measurement data!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                return;
            end
            
        end
    end
end

if get(handles.checkbox_savedata_PP,'value')==1
    if  isempty(get(handles.edit_saveId_PP,'string'))
        statusbarObj.setText('Please input ID of the save data!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
end

Interior_PostProc(IP_data,PP,handles,statusbarObj);



% --- Executes on button press in checkbox_savedata_PP.
function checkbox_savedata_PP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savedata_PP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_savedata_PP
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_saveId_PP,'Enable','on');
else
    set(handles.edit_saveId_PP,'Enable','off');
end


function edit_saveId_PP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_saveId_PP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_saveId_PP as text
%        str2double(get(hObject,'String')) returns contents of edit_saveId_PP as a double


% --- Executes during object creation, after setting all properties.
function edit_saveId_PP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_saveId_PP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_closefig.
function pushbutton_closefig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_closefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keep_guis;

function keep_guis
fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
for fh = fig_h
    uih = findobj( fh, 'Type', 'uicontrol' );
    if isempty( uih )
        delete( fh );
    end
end


function edit_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_X as text
%        str2double(get(hObject,'String')) returns contents of edit_X as a double
param=str2num(get(hObject,'String'));
set(handles.edit_X,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_X);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_z1z2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_z1z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_z1z2 as text
%        str2double(get(hObject,'String')) returns contents of edit_z1z2 as a double
param=str2num(get(hObject,'String'));
set(handles.edit_z1z2,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_z1z2);
else
    statusbarObj.setText('');
end



% --- Executes during object creation, after setting all properties.
function edit_z1z2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_z1z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_savedata.
function checkbox_savedata_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_savedata
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_saveId,'visible','on');
    set(handles.text44,'visible','on');
else
    set(handles.edit_saveId,'visible','off');
    set(handles.text44,'visible','off');
end


function edit_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_T as text
%        str2double(get(hObject,'String')) returns contents of edit_T as a double
param=str2num(get(hObject,'String'));
set(handles.edit_T,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_T);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_saveId_Callback(hObject, eventdata, handles)
% hObject    handle to edit_saveId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_saveId as text
%        str2double(get(hObject,'String')) returns contents of edit_saveId as a double


% --- Executes during object creation, after setting all properties.
function edit_saveId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_saveId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();

statusbarObj.setText('loading...');

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathNow,'\','*.mat'],'Load Data file');
if directory~=0
    %cd(directory);
    temp=load([directory,file_name]);
    
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    
    
    
    
    if isempty(my_data)|| ~isfield(my_data,'dteta')
        set(handles.popupmenu3,'value',1);
        statusbarObj.setText('No data, please load!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    else
        set(handles.popupmenu3,'Value',2)
        set(handles.pushbutton5,'Enable','off')
        
        x=my_data.x;
        T=my_data.time;
        projdir=my_data.Proj.projdir;
        minbathy=min(my_data.bathy);
        maxeta=max(max(my_data.eta));
        dzSuggest=(maxeta-minbathy)./50;
        PP_Z=[minbathy;maxeta];
        
        additemtolistbox(handles.popupmenu_CurDir,projdir);
        
        if isfield(my_data,'Fbdy')
        Fbdy=my_data.Fbdy;
        PP_X=[x(1)+Fbdy(1);x(end)-Fbdy(2)];
        else
        PP_X=[x(1);x(end)]; 
        end
        set(handles.edit_X,'Userdata',PP_X);
        set(handles.edit_X,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
        set(handles.edit_xcoarse,'Userdata',1);
        set(handles.edit_xcoarse,'String','1');
        PP_T=[T(1);T(end)];
        set(handles.edit_T,'Userdata',PP_T);
        set(handles.edit_T,'String',[num2str(T(1)),';',num2str(T(end))]);
        set(handles.edit_tcoarse,'Userdata',1);
        set(handles.edit_tcoarse,'String','1');
        
        set(handles.edit_z1z2,'Userdata',PP_Z);
        set(handles.edit_z1z2,'String',[num2str(PP_Z(1),2),';',num2str(PP_Z(2),2)]);
        set(handles.edit_dz,'Userdata',dzSuggest);
        set(handles.edit_dz,'String',num2str(dzSuggest,2));
        
        
        set(findall(handles.uipanel1, '-property', 'enable'), 'enable', 'on');
        set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
        set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
        set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
        set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
        
        set(handles.pushbutton5,'Userdata',my_data);
        %cd(pathnow);
        statusbarObj.setText([file_name,' has been loaded']);
    end
else
    statusbarObj.setText('');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
contents = cellstr(get(hObject,'String'));
set(handles.popupmenu3,'value',1);
Cont=contents{1};

Data_Input_load =findobj('Tag','pushbutton5');

if strcmp(Cont,'Other')
    set(Data_Input_load,'Enable','on');
else
    set(Data_Input_load,'Enable','off');
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calculate_IP.
function pushbutton_calculate_IP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate_IP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');

if isfield(handles.input{1},'dteta')
    timeIP=handles.input{1}.time;
    x=handles.input{1}.x;
    k=handles.input{1}.k;
    nupeak=handles.input{1}.nupeak;
    Om=handles.input{1}.Om;
    bathy=handles.input{1}.bathy;
    eta=handles.input{1}.eta;
    u=handles.input{1}.u;
    dtphihat=handles.input{1}.dtphihat;
    dteta=handles.input{1}.dteta;
    savename=handles.input{1}.savename;
    IP.savename=savename;
    if isfield(handles.input{1},'H_min')
        H_min=handles.input{1}.H_min;
    else
        H_min=[];
    end
else
    
    IP_inputdata=get(handles.pushbutton5,'Userdata');
    %assignin('base','IP_data',IP_data);
    
    if isempty(IP_inputdata)
        set(handles.popupmenu3,'Value',1);
        set(handles.pushbutton5,'enable','on');
        statusbarObj.setText('No data, please load!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    else
        timeIP=IP_inputdata.time;
        x=IP_inputdata.x;
        k=IP_inputdata.k;
        nupeak=IP_inputdata.nupeak;
        Om=IP_inputdata.Om;
        bathy=IP_inputdata.bathy;
        eta=IP_inputdata.eta;
        u=IP_inputdata.u;
        dtphihat=IP_inputdata.dtphihat;
        dteta=IP_inputdata.dteta;
        savename=IP_inputdata.savename;
        IP.savename=savename;
        if isfield(IP_inputdata,'H_min')
            H_min=IP_inputdata.H_min;
        else
            H_min=[];
        end
        
    end
    
end
Xinterv=get(handles.edit_X,'Userdata');
xcoarse=get(handles.edit_xcoarse,'Userdata');
TintervSaved=get(handles.edit_T,'Userdata');
tcoarse=get(handles.edit_tcoarse,'Userdata');
zinterv=get(handles.edit_z1z2,'Userdata');
dz=get(handles.edit_dz,'Userdata');
save_check=get(handles.checkbox_savedata,'value');
SaveID=get(handles.edit_saveId,'string');




if isempty(TintervSaved) || length(TintervSaved)~=2
    statusbarObj.setText('Specify time interval!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end


if TintervSaved(1)>TintervSaved(2)
    statusbarObj.setText(['Error:Wrong input format! [t_start;t_end]']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if TintervSaved(1)<timeIP(1)
    statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(timeIP(1)),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if TintervSaved(2)>timeIP(end)
    statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(timeIP(end)),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if isempty(tcoarse) || length(tcoarse)~=1
    statusbarObj.setText('Specify t coarse.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if fix(tcoarse)~=tcoarse
    statusbarObj.setText('Error: tcoarse is not an integer.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if tcoarse<0
    statusbarObj.setText('Error: tcoarse is negative.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if tcoarse==0
    statusbarObj.setText('Error: tcoarse is zero.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if isempty(Xinterv) || length(Xinterv)~=2
    statusbarObj.setText('Specify x interval.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if Xinterv(1)>Xinterv(2)
    statusbarObj.setText('Error: wrong input format of x interval![x_left,x_right]');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if Xinterv(1)<x(1)
    statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(x(1)),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if Xinterv(2)>x(end)
    statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(x(end)),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end


if isempty(xcoarse) || length(xcoarse)~=1
    statusbarObj.setText('Specify x coarse.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if fix(xcoarse)~=xcoarse
    statusbarObj.setText('Error: xcoarse is not an integer.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if xcoarse<0
    statusbarObj.setText('Error: xcoarse is negative.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if xcoarse==0
    statusbarObj.setText('Error: xcoarse is zero.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end


if isempty(zinterv) || length(zinterv)~=2
    statusbarObj.setText('Specify vertical interval!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if  zinterv(1)>zinterv(2)
    statusbarObj.setText('Error: wrong input format of z interval![z_bottom,z_top]');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if zinterv(1)<min(bathy)
    statusbarObj.setText(['Error:Input value z< minimum value (z=',num2str(min(bathy)),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end  

maxZ=max(max(eta));
maxZ=maxZ+maxZ/4;
if zinterv(2)>maxZ
    statusbarObj.setText(['Error:Input value z> maximum value (z=',num2str(maxZ),')']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if isempty(dz) || length(dz)~=1
    statusbarObj.setText('Specify the vertical grid size (dz).');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if  dz<0
    statusbarObj.setText('Error: dz is negative ');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if  dz==0
    statusbarObj.setText('Error: dz is zero');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if save_check==1
if  isempty(SaveID)
    statusbarObj.setText('Specify the ID of data.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
end

indtI=closest(timeIP,TintervSaved(1));
indtF=closest(timeIP,TintervSaved(2));
dtStep=tcoarse;
IP.timeIP=timeIP(indtI:dtStep:indtF);
etaNow=eta(indtI:dtStep:indtF,:);
uNow  =u(indtI:dtStep:indtF,:);
dtphihatNow=dtphihat(indtI:dtStep:indtF,:);
dtetaNow=dteta(indtI:dtStep:indtF,:);
Indx1=closest(x,Xinterv(1));
Indx2=closest(x,Xinterv(2));
IP.meanEta=mean(etaNow(:,Indx1:Indx2));
Xinterv1 =[Xinterv(1);Xinterv(2);xcoarse];
zinterv1 =[zinterv(1);dz;zinterv(2)];
IDstop=0;
[Ptot,IP.dtPHI,IP.dxPHI,IP.dzPHI,IP.dt_dxPHI,IP.dt_dzPHI,IP.Z,IP.X,IDstop]=...
    InteriorProperties(x,Xinterv1,zinterv1,k,Om,nupeak,bathy,etaNow,uNow,...
    dtphihatNow,dtetaNow,IP.timeIP,H_min,IDstop);



[statusbarObj]=JavaFrame_handling();

savedata_Id=get(handles.checkbox_savedata,'value');

if IDstop==0;
if savedata_Id==1
    
    statusbarObj.setText('saving data..');
    
    
    CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
    projdir=CurDir{get(handles.popupmenu_CurDir,'Value')};
    
    global OS
    if strcmp(OS,'Windows')
        sf_savename = [projdir,'\',savename,'\'];
    else
        sf_savename = [projdir,'/',savename,'/'];
    end
    
    if ~isdir(sf_savename)
        mkdir(sf_savename);
    end
    
    save ('-v7.3',[sf_savename,savename,'_InteriorVars_',SaveID,'.mat'],'IP')
    
    statusbarObj.setText(['data saved as ',savename,'_InteriorVars_',SaveID,'.mat']);
    pause(1)
end

if ~isdeployed
    assignin('base','IP',IP);
end

IP.Ptot=Ptot;

set(handles.popupmenu_Input,'Value',2);
set(handles.pushbutton_load_input,'enable','off');
set(handles.pushbutton_load_input,'Userdata',IP);

Tinterv=[IP.timeIP(1);IP.timeIP(end)];
set(handles.edit_PP2_T,'Userdata',Tinterv);

set(handles.edit_PP2_T,'String',[num2str(roundn(Tinterv(1),-2)),';',num2str(roundn(Tinterv(2),-2))]);
Xinterv=[IP.X(1);IP.X(end)];
set(handles.edit_PP2_X,'Userdata',Xinterv);
set(handles.edit_PP2_X,'String',[num2str(roundn(Xinterv(1),-2)),';',num2str(roundn(Xinterv(2),-2))]);
Zinterv=[IP.Z(1);IP.Z(end)];
set(handles.edit_PP2_Z,'Userdata',Zinterv);
set(handles.edit_PP2_Z,'String',[num2str(roundn(Zinterv(1),-2)),';',num2str(roundn(Zinterv(2),-2))]);

set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'on')
set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton1,'value',0);
set(handles.radiobutton1,'enable','on');
set(handles.radiobutton2,'value',0);
set(handles.radiobutton2,'enable','on');
set(handles.radiobutton3,'value',0);
set(handles.radiobutton3,'enable','on');

statusbarObj.setText('ready for post-processing');
else
statusbarObj.setText('Terminated by user.');
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.blue);
end



function edit_tcoarse_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tcoarse as text
%        str2double(get(hObject,'String')) returns contents of edit_tcoarse as a double
param=str2num(get(hObject,'String'));
set(handles.edit_tcoarse,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_tcoarse);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_tcoarse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xcoarse_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xcoarse as text
%        str2double(get(hObject,'String')) returns contents of edit_xcoarse as a double
param=str2num(get(hObject,'String'));
set(handles.edit_xcoarse,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xcoarse);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_xcoarse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dz as text
%        str2double(get(hObject,'String')) returns contents of edit_dz as a double
param=str2num(get(hObject,'String'));
set(handles.edit_dz,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_dz);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
param=str2num(get(hObject,'String'));
set(handles.edit34,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit34);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_setting_ylim.
function checkbox_setting_ylim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_setting_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_setting_ylim
Id=get(hObject,'Value');
if Id==1
    set(handles.edit34,'enable','on');
else
    set(handles.edit34,'enable','off');
end

% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_xlim,'enable','on');
else
    set(handles.edit_xlim,'enable','off');
end


function edit_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_colorbar as text
%        str2double(get(hObject,'String')) returns contents of edit_colorbar as a double
param=str2num(get(hObject,'String'));
set(handles.edit_colorbar,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_colorbar);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_colorbar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_colorbar,'enable','On')
else
    set(handles.edit_colorbar,'enable','Off')
end

% --- Executes on button press in checkbox_PP2_GIF.
function checkbox_PP2_GIF_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP2_GIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP2_GIF
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP2_GIF,'enable','on');
else
    set(handles.edit_PP2_GIF,'enable','off');
end


function edit_PP2_GIF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP2_GIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP2_GIF as text
%        str2double(get(hObject,'String')) returns contents of edit_PP2_GIF as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP2_GIF,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP2_GIF);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP2_GIF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP2_GIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xlim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xlim as text
%        str2double(get(hObject,'String')) returns contents of edit_xlim as a double
param=str2num(get(hObject,'String'));
set(handles.edit_xlim,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xlim);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_xlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_tcoarse.
function checkbox_tcoarse_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tcoarse
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP_tcoarse,'enable','on');
else
    set(handles.edit_PP_tcoarse,'enable','off');
end


function edit_PP_tcorase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP_tcorase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP_tcorase as text
%        str2double(get(hObject,'String')) returns contents of edit_PP_tcorase as a double


% --- Executes during object creation, after setting all properties.
function edit_PP_tcorase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP_tcorase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
Id=get(hObject,'Value');

if Id==1
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'on')
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton2,'value',0);
    set(handles.radiobutton1,'value',0);
    set(handles.radiobutton2,'enable','on');
    set(handles.radiobutton1,'enable','on');
    if get(handles.checkbox_PP3_lineplot,'value')==1
        set(handles.checkbox38,'enable','on');
        set(handles.checkbox39,'enable','on');
        if get(handles.checkbox38,'value')==1
            set(handles.edit_PP3_Z,'enable','on');
        else
            set(handles.edit_PP3_Z,'enable','off');
        end
        if get(handles.checkbox39,'value')==1
            set(handles.edit_dZmeas,'enable','on');
        else
            set(handles.edit_dZmeas,'enable','off');
        end
        
    else
        set(handles.checkbox38,'enable','off');
        set(handles.checkbox39,'enable','off');
        set(handles.edit_PP3_Z,'enable','off');
        set(handles.edit_dZmeas,'enable','off');
    end
    
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'on')
    set(handles.checkbox_PP2_GIF,'enable','off')
    set(handles.edit_PP2_GIF,'enable','off')
    set(handles.checkbox_tcoarse,'enable','off')
    set(handles.edit_PP_tcoarse,'enable','off')
    if get(handles.checkbox22,'value')==0
        set(handles.edit_colorbar,'enable','off')
    end
    if get(handles.checkbox25,'value')==0
        set(handles.edit_xlim,'enable','off')
    end
    if get(handles.checkbox_setting_ylim,'value')==0
        set(handles.edit34,'enable','off')
    end
    if get(handles.checkbox_savedata_PP,'value')==0
        set(handles.edit_saveId_PP,'enable','off')
    end
    
    if get(handles.checkbox_tcoarse,'value')==0
        set(handles.edit_PP_tcoarse,'enable','off')
    end
    
     if get(handles.checkbox_savefigure,'value')==0
        set(handles.popupmenu_saveFig,'enable','off')
    end
    
else
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton3,'enable','on');
end
set(handles.edit_PP1_line_X,'enable','off')
set(handles.edit_dZmeas,'userdata',1);
set(handles.edit_dZmeas,'string','1');


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


% --- Executes on button press in pushbutton_measdata.
function pushbutton_measdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_measdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Input_gui = guidata(hObject);

[statusbarObj]=JavaFrame_handling();

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathNow,'\','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
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
    set(Input_gui.pushbutton_measdata,'Userdata',my_data);
    statusbarObj.setText([file_name,' has been loaded']);
    if isstruct(my_data)
        statusbarObj.setText('The data is not in correct format');
        set(Input_gui.pushbutton_measdata,'Userdata',[]);  
    else

        set(handles.edit26,'string',num2str(my_data(1,1)));
        set(handles.edit26,'enable','off');
    end
else
    statusbarObj.setText('');
end
set(handles.edit26,'enable','off');

% --- Executes on button press in checkbox_PP3_densityplot.
function checkbox_PP3_densityplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP3_densityplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP3_densityplot
Id=get(handles.checkbox_PP3_densityplot,'value');
if get(handles.checkbox_PP3_densityplot,'value')...
        +get(handles.checkbox_PP3_errorplot,'value')>0
    set(handles.checkbox22,'enable','on')
else
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end
if Id==1 
 set(handles.checkbox_PP3_lineplot,'value',0);
 set(handles.checkbox39,'enable','off');
 set(handles.checkbox38,'enable','off');
 set(handles.edit_dZmeas,'enable','off');
 set(handles.edit_PP3_Z,'enable','off');
 set(handles.checkbox_timeMeanVariation,'value',0)
end


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_line_X as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_line_X as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_PP3_lineplot.
function checkbox_PP3_lineplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP3_lineplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP3_lineplot
Id=get(hObject,'Value');
if Id==1
    set(handles.checkbox39,'value',	1);
    set(handles.checkbox38,'value',	0);
    set(handles.checkbox38,'enable','on');
    set(handles.checkbox39,'enable','on');
    if get(handles.checkbox38,'value')==1
        set(handles.edit_PP3_Z,'enable','on');
    else
        set(handles.edit_PP3_Z,'enable','off');
    end
    if get(handles.checkbox39,'value')==1
        set(handles.edit_dZmeas,'enable','on');
    else
        set(handles.edit_dZmeas,'enable','off');
    end
    
else
    set(handles.checkbox38,'enable','off');
    set(handles.checkbox39,'enable','off');
    set(handles.edit_PP3_Z,'enable','off');
    set(handles.edit_dZmeas,'enable','off');
    
end

if get(handles.checkbox_PP3_densityplot,'value')...
        +get(handles.checkbox_PP3_errorplot,'value')==0
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end

if Id==1 
 set(handles.checkbox_PP3_densityplot,'value',0);
 set(handles.checkbox_PP3_errorplot,'value',0);
 set(handles.checkbox_timeMeanVariation,'value',0)
end


% --- Executes on button press in checkbox38.
function checkbox38_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox38
Id=get(hObject,'Value');
if Id==1
    set(handles.checkbox39,'value',0);
    set(handles.edit_dZmeas,'enable','off');
    set(handles.edit_PP3_Z,'enable','on');
else
    set(handles.checkbox39,'value',1);
    set(handles.edit_dZmeas,'enable','on');
    set(handles.edit_PP3_Z,'enable','off');
end

% --- Executes on button press in checkbox39.
function checkbox39_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox39
Id=get(hObject,'Value');
if Id==1
    set(handles.checkbox38,'value',0);
    set(handles.edit_dZmeas,'enable','on');
    set(handles.edit_PP3_Z,'enable','off');
    
else
    set(handles.checkbox38,'value',1);
    set(handles.edit_dZmeas,'enable','off');
    set(handles.edit_PP3_Z,'enable','on');
end


function edit_dZmeas_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dZmeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dZmeas as text
%        str2double(get(hObject,'String')) returns contents of edit_dZmeas as a double
param=str2num(get(hObject,'String'));
set(handles.edit_dZmeas,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_dZmeas);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_dZmeas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dZmeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_PP3_errorplot.
function checkbox_PP3_errorplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP3_errorplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP3_errorplot
if get(handles.checkbox_PP3_densityplot,'value')...
        +get(handles.checkbox_PP3_errorplot,'value')>0
    set(handles.checkbox22,'enable','on')
else
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end
if get(handles.checkbox_PP3_errorplot,'value')==1 
 set(handles.checkbox_PP3_lineplot,'value',0);
 set(handles.checkbox39,'enable','off');
 set(handles.checkbox38,'enable','off');
 set(handles.edit_dZmeas,'enable','off');
 set(handles.edit_PP3_Z,'enable','off');
 set(handles.checkbox_timeMeanVariation,'value',0)
end


function edit_PP3_Z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP3_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP3_Z as text
%        str2double(get(hObject,'String')) returns contents of edit_PP3_Z as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP3_Z,'Userdata',param);

% --- Executes during object creation, after setting all properties.
function edit_PP3_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP3_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_timeMeanVariation.
function checkbox_timeMeanVariation_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_timeMeanVariation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_timeMeanVariation
Id=get(handles.checkbox_timeMeanVariation,'value');
if Id==1 
 set(handles.checkbox_PP3_lineplot,'value',0);
 set(handles.checkbox39,'enable','off');
 set(handles.checkbox38,'enable','off');
 set(handles.edit_dZmeas,'enable','off');
 set(handles.edit_PP3_Z,'enable','off');
 set(handles.checkbox_PP3_densityplot,'value',0)
 set(handles.checkbox_PP3_errorplot,'value',0)
end

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
Id=get(hObject,'Value');

if Id==1
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'on')
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton3,'value',0);
    set(handles.radiobutton1,'value',0);
    set(handles.radiobutton3,'enable','on');
    set(handles.radiobutton1,'enable','on');
    set(handles.edit_PP2_GIF,'enable','off')
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'on')
    set(handles.checkbox22,'enable','on')
    set(handles.checkbox25,'enable','off')
    if get(handles.checkbox22,'value')==0
        set(handles.edit_colorbar,'enable','off')
    else
        set(handles.edit_colorbar,'enable','on')
    end
    
    
    set(handles.edit_xlim,'enable','off')
    
    set(handles.checkbox_setting_ylim,'enable','off')
    set(handles.edit34,'enable','off')
  
    
    if get(handles.checkbox_PP2_GIF,'value')==0
        set(handles.edit_PP2_GIF,'enable','off')
    end
    if get(handles.checkbox_savedata_PP,'value')==0
        set(handles.edit_saveId_PP,'enable','off')
    end
    
    if get(handles.checkbox_tcoarse,'value')==0
        set(handles.edit_PP_tcoarse,'enable','off')
    end
     if get(handles.checkbox_savefigure,'value')==0
        set(handles.popupmenu_saveFig,'enable','off')
    end
    
else
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton2,'enable','on');
end




function edit_PP1_line_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_line_X as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_line_X as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_line_X,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_line_X);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_line_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PP2_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP2_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP2_X as text
%        str2double(get(hObject,'String')) returns contents of edit_PP2_X as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP2_X,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP2_X);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP2_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP2_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PP2_Z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP2_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP2_Z as text
%        str2double(get(hObject,'String')) returns contents of edit_PP2_Z as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP2_Z,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP2_Z);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_PP2_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP2_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PP1_XT_Z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_XT_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_XT_Z as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_XT_Z as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_XT_Z,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_XT_Z);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_XT_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_XT_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_XT.
function checkbox_XT_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_XT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_XT
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_XT_Z,'enable','on');
else
    set(handles.edit_PP1_XT_Z,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')>0
    set(handles.checkbox22,'enable','on')
else
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end


function edit_PP1_ZT_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_ZT_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_ZT_X as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_ZT_X as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_ZT_X,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_ZT_X);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_ZT_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_ZT_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ZT.
function checkbox_ZT_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ZT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ZT
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_ZT_X,'enable','on');
else
    set(handles.edit_PP1_ZT_X,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')>0
    set(handles.checkbox22,'enable','on')
else
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
Id=get(hObject,'Value');

if Id==1
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'on')
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton2,'value',0);
    set(handles.radiobutton3,'value',0);
    set(handles.radiobutton2,'enable','on');
    set(handles.radiobutton3,'enable','on');
    if get(handles.checkbox_XZ,'value')==0
        set(handles.edit_PP1_XZ_time,'enable','off');
    else
        set(handles.edit_PP1_XZ_time,'enable','on');
    end
    if get(handles.checkbox_ZT,'value')==0
        set(handles.edit_PP1_ZT_X,'enable','off');
    else
        set(handles.edit_PP1_ZT_X,'enable','on');
    end
    if get(handles.checkbox_XT,'value')==0
        set(handles.edit_PP1_XT_Z,'enable','off');
    else
        set(handles.edit_PP1_XT_Z,'enable','off');
    end
    if get(handles.checkbox_PP1_line_X,'value')==0
        set(handles.edit_PP1_line_X,'enable','off');
    else
        set(handles.edit_PP1_line_X,'enable','off');
    end
    if get(handles.checkbox_PP1_line_Z,'value')==0
        set(handles.edit_PP1_line_Z,'enable','off');
    else
        set(handles.edit_PP1_line_X,'enable','off');
    end
    if get(handles.checkbox_PP1_line_T,'value')==0
        set(handles.edit_PP1_line_T,'enable','off');
    else
        set(handles.edit_PP1_line_X,'enable','off');
    end
    
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'on')
    set(handles.checkbox_PP2_GIF,'enable','off')
    set(handles.edit_PP2_GIF,'enable','off')
    if get(handles.checkbox22,'value')==0
        set(handles.edit_colorbar,'enable','off')
    end
    if get(handles.checkbox25,'value')==0
        set(handles.edit_xlim,'enable','off')
    end
    if get(handles.checkbox_setting_ylim,'value')==0
        set(handles.edit34,'enable','off')
    end
    if get(handles.checkbox_savedata_PP,'value')==0
        set(handles.edit_saveId_PP,'enable','off')
    end
    
    set(handles.checkbox_tcoarse,'enable','off');
    set(handles.edit_PP_tcoarse,'enable','off')
   
     if get(handles.checkbox_savefigure,'value')==0
        set(handles.popupmenu_saveFig,'enable','off')
    end
    
else
    set(findall(handles.uipanel6, '-property', 'enable'), 'enable', 'off')
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
    set(handles.radiobutton1,'enable','on');
end


function edit_PP1_XZ_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_XZ_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_XZ_time as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_XZ_time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_XZ_time,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_XZ_time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_XZ_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_XZ_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_XZ.
function checkbox_XZ_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_XZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_XZ
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_XZ_time,'enable','on');
else
    set(handles.edit_PP1_XZ_time,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')>0
    set(handles.checkbox22,'enable','on')
else
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end

% --- Executes on button press in checkbox_PP1_line_X.
function checkbox_PP1_line_X_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_line_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_line_X
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_line_X,'enable','on');
else
    set(handles.edit_PP1_line_X,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')==0
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end


function edit_PP2_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP2_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP2_T as text
%        str2double(get(hObject,'String')) returns contents of edit_PP2_T as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP2_T,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP2_T);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_PP2_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP2_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PP1_line_Z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_line_Z as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_line_Z as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_line_Z,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_line_Z);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_line_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_PP1_line_Z.
function checkbox_PP1_line_Z_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_line_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_line_Z
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_line_Z,'enable','on');
else
    set(handles.edit_PP1_line_Z,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')==0
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end

% --- Executes on button press in checkbox_PP1_line_T.
function checkbox_PP1_line_T_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_line_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_line_T

Id=get(hObject,'Value');
if Id==1
    set(handles.edit_PP1_line_T,'enable','on');
else
    set(handles.edit_PP1_line_T,'enable','off');
end
if get(handles.checkbox_XZ,'value')+get(handles.checkbox_ZT,'value')...
        +get(handles.checkbox_XT,'value')==0
    set(handles.checkbox22,'enable','off')
    set(handles.edit_colorbar,'enable','off')
end


function edit_PP1_line_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP1_line_T as text
%        str2double(get(hObject,'String')) returns contents of edit_PP1_line_T as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP1_line_T,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP1_line_T);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP1_line_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP1_line_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reset_handles(handles,Flag)
if Flag==1
    [statusbarObj]=JavaFrame_handling();
    statusbarObj.setText('');
    
    JavProgressBar;
    jProgressBar.setVisible(0);
end

set(handles.edit_X,'Userdata',[]);
set(handles.edit_X,'String','');

set(handles.edit_z1z2,'Userdata',[]);
set(handles.edit_z1z2,'String','');
set(handles.popupmenu_Input,'Value',1);
set(handles.pushbutton_load_input,'Userdata',[]);
set(handles.edit_PP1_XZ_time,'Userdata',[]);
set(handles.edit_PP1_ZT_X,'Userdata',[]);
set(handles.edit_PP1_XT_Z,'Userdata',[]);
set(handles.edit_PP1_XZ_time,'String','');
set(handles.edit_PP1_ZT_X,'String','');
set(handles.edit_PP1_XT_Z,'String','');
set(handles.edit_PP2_T,'Userdata',[]);
set(handles.edit_PP2_T,'String','');
set(handles.edit_PP2_X,'Userdata',[]);
set(handles.edit_PP2_X,'String','');
set(handles.edit_PP2_Z,'Userdata',[]);
set(handles.edit_PP2_Z,'String','');
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',0);
set(handles.checkbox_XZ,'value',0);
set(handles.checkbox_ZT,'value',0);
set(handles.checkbox_XT,'value',0);
set(handles.checkbox_PP1_line_X,'value',0);
set(handles.checkbox_PP1_line_Z,'value',0);
set(handles.checkbox_PP1_line_T,'value',0);
set(handles.checkbox_P,'value',0);
set(handles.checkbox_P_dyn,'value',0);
set(handles.checkbox_P_NonlinDyn,'value',0);
set(handles.checkbox_V,'value',0);
set(handles.checkbox_Vx,'value',0);
set(handles.checkbox_Vz,'value',0);

if Flag==1
    statusbarObj.setText('Clear');
end



function edit_PP_tcoarse_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PP_tcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PP_tcoarse as text
%        str2double(get(hObject,'String')) returns contents of edit_PP_tcoarse as a double
param=str2num(get(hObject,'String'));
set(handles.edit_PP_tcoarse,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_PP_tcoarse);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_PP_tcoarse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PP_tcoarse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_savefigure.
function checkbox_savefigure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savefigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_savefigure
Id=get(hObject,'Value');
if Id==1
set(handles.popupmenu_saveFig,'Enable','on')    
else
set(handles.popupmenu_saveFig,'Enable','off')    
end

% --- Executes on selection change in popupmenu_saveFig.
function popupmenu_saveFig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_saveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_saveFig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_saveFig


% --- Executes during object creation, after setting all properties.
function popupmenu_saveFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_saveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
