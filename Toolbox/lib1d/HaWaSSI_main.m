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
% !
% !   HAWASSI AB 1D MAIN
% !
% !--------------------------------------------------------------------------------------
% !

function varargout = HaWaSSI_main(varargin)
% HAWASSI_MAIN MATLAB code for HaWaSSI_main.fig
%      HAWASSI_MAIN, by itself, creates a new HAWASSI_MAIN or raises the existing
%      singleton*.
%
%      H = HAWASSI_MAIN returns the handle to a new HAWASSI_MAIN or the handle to
%      the existing singleton*.
%
%      HAWASSI_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HAWASSI_MAIN.M with the given input arguments.
%
%      HAWASSI_MAIN('Property','Value',...) creates a new HAWASSI_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HaWaSSI_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HaWaSSI_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HaWaSSI_main

% Last Modified by GUIDE v2.5 11-Dec-2019 11:38:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HaWaSSI_main_OpeningFcn, ...
    'gui_OutputFcn',  @HaWaSSI_main_OutputFcn, ...
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


% --- Executes just before HaWaSSI_main is made visible.
function HaWaSSI_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HaWaSSI_main (see VARARGIN)

% Choose default command line output for HaWaSSI_main
global IDdemo


handles.output = hObject;
handles.input=varargin{1};

          
reset_handles(handles,0);
% Update handles structure
guidata(hObject, handles);


%set(hObject,'units','normalized','Position',[0.29 0.3 0.37 0.65])
%set(0,'ScreenPixelsPerInch',96)
pathstr=handles.input.pathnow;
projectdirectory=handles.input.projectdirectory;
additemtolistbox(handles.popupmenu_CurDir,projectdirectory);


try
    axes(handles.Logo2);     % Show the LabMath.
    
    %javaFrame = get(hObject,'JavaFrame');
    
    if ~isdeployed
        imshow([pathstr,'\Toolbox\lib1d\Misc\logo\hawassi_black.jpg']); %
        %   javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
    else
        imshow('hawassi_black.jpg');
        %  javaFrame.setFigureIcon(javax.swing.ImageIcon('\Misc\logo\icon_16.png'));
    end
catch
end

IDdemo=handles.input.IdDemo;
if IDdemo==0
    set(handles.figure1,'Name',['HAWASSI AB1: Main']);
else
    set(handles.figure1,'Name',['HAWASSI AB1 demo: Main']);
end


% set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
% set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')
% set(handles.text60,'FontName','symbol');
% set(handles.text70,'FontName','symbol');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if IDdemo==1
    %%%%%%%%%%%%%%%%%%%%%%%%Demo-Version%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.Interior,'Enable','off');
    set(handles.checkbox_interior,'Enable','off');
    set(handles.checkbox_Breaking,'Enable','off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

initialization_gui(handles)


global PreProc FlagWarn
PreProc=0;FlagWarn=0;
% UIWAIT makes HaWaSSI_main wait for user response (see UIRESUME)
%uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = HaWaSSI_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = [];%handles.output;

function initialization_gui(handles)
%%%%%Set tips
set(handles.popupmenu_CurDir,'tooltipString','The directory list')
set(handles.pushbutton_browse_Dir,'tooltipString','Specify the directory')
set(handles.edit_savename,'tooltipString',sprintf('Specify the Project name.\n The output will be saved in the directory with this name.'))
set(handles.edit_note,'tooltipString','Type your note here.')
set(handles.WaveType,'tooltipString','Choose an input.')
set(handles.WT_Other_data,'tooltipString','Load an influx signal.')

set(handles.edit_Tp,'tooltipString','Peak Period');
set(handles.edit_Jsg,'tooltipString',sprintf('Peak enhancement factor.\n Default value 3.3.'));
%set(handles.edit_Jsg,'userdata',3.3);set(handles.edit_Jsg,'string','3.3');
set(handles.edit_Jsg,'backgroundcolor',[0.871 0.922 0.98]);
set(handles.checkbox_Filter,'tooltipString','Option to restrict the frequency band of the time signals')
set(handles.edit_LFHF,'tooltipString',sprintf('Provides the desired frequency band.\n Input format [Low Freq; High Freq]'))
set(handles.edit_time,'tooltipString',sprintf('Specify the time interval.\n Input format: [t_start; t_end]'));
set(handles.edit_timestep,'tooltipString','Specify the time step. Input format: dt');
set(handles.checkbox_ramp,'tooltipString','Option to apply a ramp function in the start and the end of the time signal');
set(handles.edit_ramp,'tooltipString',sprintf('Specify length of the ramp. \n Default value 4'));

set(handles.dynModel,'tooltipString',sprintf('Choose a dynamic model.\n HS:Hamiltonian System \n with an order of nonlinearity (1,2,3,4).'))
set(handles.checkbox_Breaking,'tooltipString','Option to activate the breaking model');
set(handles.Breaking_param,'tooltipString',sprintf('Input a value in [0.6;1]. \n This is a kinematic criterion [U/C]. '));
set(handles.checkbox_def_paramB,'tooltipString','Option to use the default value of the termination parameters.');
set(handles.edit_paramB_UF_T,'tooltipString',sprintf('Specify the termination parameters.  \n Default value [UF/UI;T*/Tp]=[0.2;0.5].'));
set(handles.edit_paramB_UF_T,'backgroundcolor',[0.871 0.922 0.98]);
set(handles.popupmenu_dispersion,'tooltipString','Choose a dispersion model.')
set(handles.togglebutton_InputDispersion,'tooltipString','Input a dispersion formula. Use k and d as variables.')
set(handles.togglebutton_InputGroupVel,'tooltipString','Input a group velocity formula. Use k and d as variables.')

set(handles.popupmenu_influx,'tooltipString','Choose an influx method');
set(handles.popupmenu_influx_direction,'tooltipString',sprintf('Choose the propagation direction.\n Uni+ :uni direction to the right \n Uni- :uni direction to the left \n Bi: bi-direction'));
set(handles.edit_xinterv,'tooltipString',sprintf('Specify the horizontal spatial interval. \n Input format: x_start; x_end'));
set(handles.edit_Xinflux,'tooltipString','Specify the influx position. Input format: x_influx');
set(handles.edit_cutfracK,'tooltipString','Specify the cut fraction wavenumber for anti aliasing. Input an integer value.');

set(handles.popupmenu_bathy,'tooltipString','Choose a bathymetry');
set(handles.pushbutton_load_file_bathy,'tooltipString','Load a function file of the bathymetry');
set(handles.edit_depth_F,'tooltipString','Specify the depth')
set(handles.checkbox_bottom_friction,'tooltipString','Option to apply the friction.');
set(handles.togglebutton_friction,'tooltipstring','Input friction coefficients');
set(handles.checkbox_wall,'tooltipString', 'Option to activate a wall')
set(handles.edit_Xwall,'tooltipString','Specify the wall position.')
set(handles.popupmenu_WallType,'tooltipString','Choose a type of wall. uniform or freq. dependent');
set(handles.edit_Refl_percent,'tooltipString','Specify the reflection coefficient in [0,1]');
set(handles.togglebutton_inputEq,'tooltipString','Input formula the reflection coefficient. use f as frequency [rad/s] variable. ]');


set(handles.edit_FBL,'tooltipString','Specify the length of Fourier boundary [L;R] [m]');
set(handles.edit_Nonlin_Adj,'tooltipString','Specify the length of nonlinear adjustment (*peak wavelength)');


set(handles.checkbox_interior,'tooltipString','For post-processing data in specified interval')
set(handles.edit_Interior_Time,'tooltipString',sprintf('Specify the time interval. \n Input format: t_start; t_end'));
set(handles.edit_dtIP,'tooltipString',sprintf('Specify the timestep. \n Input format: a*dt, a=1,2,..., dt is the time-step of the wave-input'));

set(handles.checkbox_partition,'tooltipString',sprintf('Default: Number of partition is obtained automatically \n according to the size of input data and physical memory'))
set(handles.edit_partition,'tooltipString',sprintf('Specify number of the partition of the output data.\n Input an integer value [1,2, ...].'));
set(handles.checkbox_combineFiles,'tooltipString','Option to combine all parts of the output file')

set(handles.pushbutton1,'tooltipString','Start the pre-processing')
set(handles.togglebutton_logfile,'tooltipString','Show the log-file')
set(handles.pushbutton_run,'tooltipString','Start the dynamic calculation')
set(handles.popupmenu_IVP,'tooltipString',sprintf('Option for an initial condition.'))
set(handles.edit_IVP_A,'tooltipString',sprintf('Specify the amplitude'));
set(handles.edit_IVP_lambda,'tooltipString',sprintf('Specify standard deviation'));
set(handles.edit_IVP_x0,'tooltipString',sprintf('Specify the center position'));
set(handles.pushbutton_load_file_IVP,'tooltipString','Load a function file of the initial condition')


if handles.input.flagOpenProj==0
    %%% | Initially set the input box invinsible.
    set(handles.edit_savename,'String',handles.input.projectname);
    set(handles.edit_note,'String',handles.input.usernote);
    
    set(handles.pushbutton_run,'enable','off');
    set(handles.togglebutton_logfile,'enable','off');
    set(handles.togglebutton_logfile,'value',0);
    
    set(handles.WaveType,'value',1);
    set(handles.WT_Other_data,'Visible','off');
    set(handles.pushbutton_load_file_IVP,'Visible','off');
    set(handles.edit_Jsg,'Enable','off');
    set(handles.edit_Hs,'enable','off');
    set(handles.edit_Tp,'enable','off');
    set(handles.checkbox_ramp,'value',0)
    set(handles.edit_ramp,'visible','off')
    set(handles.text100,'visible','off')
    set(handles.text101,'visible','off')
    set(handles.edit_ramp,'userdata',4)
    set(handles.edit_ramp,'string','4')
    set(handles.checkbox_Filter,'Enable','off');
    set(handles.edit_LFHF,'visible','off');
    set(handles.text71,'visible','off');
    set(handles.text61,'visible','off');
    set(handles.checkbox_ramp,'Enable','off')
    set(handles.edit_ramp,'visible','off');
    set(handles.text101,'visible','off');
    set(handles.text100,'visible','off');
    set(handles.checkbox_Filter,'value',0);
    
    
    set(handles.Breaking_param,'Visible','off');
    set(handles.checkbox_def_paramB,'Value',1);
    set(handles.checkbox_def_paramB,'Visible','off');
    set(handles.edit_paramB_UF_T,'Visible','off');
    set(handles.text72,'Visible','off');
    set(handles.Breaking_param_text,'Visible','off');
    
    
    
    set(handles.text11,'Visible','On');
    set(handles.togglebutton_InputDispersion ,'Visible','off');
    set(handles.togglebutton_InputGroupVel ,'Visible','off');
    
    set(handles.edit_LFHF ,'Visible','off');
    set(handles.text71 ,'Visible','off');set(handles.text61 ,'Visible','off');
    set(handles.pushbutton_load_file_bathy ,'Visible','off');
    set(handles.edit_depth_F,'Visible','on');
    set(handles.edit_slope,'Enable','off');
    set(handles.popupmenu_depth_reff,'enable','off');
    set(handles.text112,'enable','off');
    
    set(handles.edit_Dmidref,'enable','off');
    set(handles.text118,'enable','off');
    set(handles.text119,'enable','off');
    
    set(handles.checkbox17,'value',0);
    set(handles.togglebutton_bdy_assim,'enable','off')
    set(handles.togglebutton_bdy_assim,'userdata',[])
    
    set(handles.checkbox_wind,'value',0);
    set(handles.edit_coef_wind,'enable','off');
    set(handles.edit_wind_tinterv,'enable','off');
    set(handles.edit_wind_xinterv,'enable','off');
    
    set(handles.edit_Xwall ,'Visible','off');
    set(handles.edit_Refl_percent,'Visible','off');
    set(handles.text68,'Visible','off');
    set(handles.text83,'Visible','off');
    set(handles.text84,'Visible','off');
    set(handles.text86,'Visible','off');
    set(handles.popupmenu_WallType,'Visible','off');
    set(handles.togglebutton_inputEq,'Visible','off');
    
    set(handles.edit_Interior_Time,'enable','off')
    set(handles.checkbox_partition,'value',1);
    set(handles.edit_partition,'enable','off');
    set(handles.checkbox_combineFiles,'value',1);
    set(handles.checkbox_combineFiles,'enable','off');
    
    if handles.input.ShipFlag==0
        set(handles.checkbox_ship,'value',0,'enable','off');
        set(handles.pushbutton_ship_setupGui,'enable','off')
    else
        set(handles.checkbox_ship,'value',1,'enable','on');
        set(handles.pushbutton_ship_setupGui,'enable','on')
    end
    
    datIniship.data=[];
    datIniship.user_shapedata=[];
    datIniship.FlagCloseGui=0;
    set(handles.pushbutton_ship_setupGui,'Userdata',datIniship)
    
    
    set(handles.edit_paramB_UF_T,'string','0.2;0.5')
    set(handles.edit_paramB_UF_T,'userdata',[0.2;0.5]);
    
    set(handles.edit_slope,'backgroundcolor',[0.871 0.922 0.98]);
    
    set(handles.checkbox_bottom_friction,'value',0);
    set(handles.togglebutton_friction,'Visible','off')
    set(handles.togglebutton_friction,'userdata',[])
    
    set(handles.edit_dtIP,'Enable','off');
    
    set(handles.popupmenu_IVP,'value',1)
    set(handles.edit_IVP_A,'enable','off');
    set(handles.edit_IVP_lambda,'enable','off');
    set(handles.edit_IVP_x0,'enable','off');
    set(handles.pushbutton_load_file_IVP,'visible','off')
    set(handles.edit_Nonlin_Adj,'enable','off');
    set(handles.togglebutton_inputEq,'userdata','1-min(f.^2,0.5)');
    DispEq='sign(k).*sqrt(9.81*k.*tanh(d.*k))';
    set(handles.togglebutton_InputDispersion,'userdata',DispEq);
    UgEq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)';
    set(handles.togglebutton_InputGroupVel,'userdata',UgEq);
else
    %[statusbarObj]=JavaFrame_handling();
    % statusbarObj.setText('Loading a project. . .');
    
    reset_handles(handles,0);
    input_handles_gui(handles.input.GUIinput,handles)
    % statusbarObj.setText('');
end
%


% --- Executes on selection change in WaveType.
function WaveType_Callback(hObject, eventdata, handles)
% hObject    handle to WaveType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns WaveType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from WaveType
contents = cellstr(get(hObject,'String'));
WT=contents{get(hObject,'Value')};


Hs  = findobj('Tag','edit_Hs');            %%
Tp  = findobj('Tag','edit_Tp');            %%
Jsg = findobj('Tag','edit_Jsg');            %%
WT_O_data = findobj('Tag','WT_Other_data');            %%

[statusbarObj]=JavaFrame_handling();

if strcmp(WT,'None')
    statusbarObj.setText('');
    set(Hs,'enable','off'); set(Tp,'enable','off');set(Jsg,'enable','off');
    set(WT_O_data,'Visible','off');
    set(handles.checkbox_Filter,'Enable','off');
    set(handles.edit_LFHF,'visible','off');
    set(handles.text71,'visible','off');
    set(handles.text61,'visible','off');
    set(handles.checkbox_ramp,'Enable','off')
    set(handles.edit_ramp,'visible','off');
    set(handles.text101,'visible','off');
    set(handles.text100,'visible','off');
    set(handles.text_Hs,'visible','off');
    set(handles.text_AInflux,'visible','on');
    set(handles.checkbox_ramp,'value',0)
    
    
elseif strcmp(WT,'Jonswap')
    statusbarObj.setText('');
    set(Hs,'enable','on'); set(Tp,'enable','on');set(Jsg,'enable','on');
    set(WT_O_data,'Visible','off');
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1;
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    set(handles.checkbox_ramp,'value',1)
    set(handles.edit_ramp,'visible','on');
    set(handles.text101,'visible','on');
    set(handles.text100,'visible','on');
    
    
    set(handles.text_Hs,'visible','on');
    set(handles.text_AInflux,'visible','off');
    set(handles.edit_Hs,'tooltipString','Significant wave height');
    
    
    
elseif strcmp(WT,'Harmonic')
    statusbarObj.setText('');
    set(Hs,'enable','on'); set(Tp,'enable','on');set(Jsg,'enable','off');
    set(WT_O_data,'Visible','off');
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1;
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    set(handles.checkbox_ramp,'value',1)
    set(handles.edit_ramp,'visible','on');
    set(handles.text101,'visible','on');
    set(handles.text100,'visible','on');
    set(handles.text_Hs,'visible','off');
    set(handles.text_AInflux,'visible','on');
    set(handles.edit_Hs,'tooltipString','Amplitude');
    
elseif strcmp(WT,'User-defined')
    statusbarObj.setText('');
    set(WT_O_data,'Visible','on');
    set(Hs,'enable','off'); set(Tp,'enable','off');set(Jsg,'enable','off');
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1;
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    set(handles.checkbox_ramp,'value',1)
    set(handles.edit_ramp,'visible','on');
    set(handles.text101,'visible','on');
    set(handles.text100,'visible','on');
    
    set(handles.text_Hs,'visible','off');
    set(handles.text_AInflux,'visible','on');
end


% --- Executes during object creation, after setting all properties.
function WaveType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PreProc
clc;
PreProc=1;
if get(handles.togglebutton_logfile,'value')==1
    delete(log_file);
    set(handles.togglebutton_logfile,'value',0);
end

Initialization(handles);





% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PreProc
set(handles.togglebutton_logfile,'userdata',[]); %flag data for logfile
if PreProc==1
    PreProc=0;
    set(handles.togglebutton_logfile,'enable','off');
    set(handles.pushbutton1,'enable','off');
    set(handles.togglebutton_logfile,'value',0);
    set(handles.pushbutton_run,'enable','off');
    
    jFrame = get(handle(gcf),'JavaFrame');
    try
        jRootPane = jFrame.fHG1Client.getWindow;    %<2014b
    catch
        jRootPane = jFrame.fHG2Client.getWindow;    %>=2014b
    end
    statusbarObj = com.mathworks.mwswing.MJStatusBar;
    jRootPane.setStatusBar(statusbarObj);
    statusbarObj.setText('');
    
    preproc_data=get(handles.pushbutton1,'userdata');
    par     =preproc_data.par;
    influx  =preproc_data.influx;
    model   =preproc_data.model;
    Proj    =preproc_data.Proj;
    bath    =preproc_data.bath;
    options =preproc_data.options;
    input   =preproc_data.input;
    IVP     =preproc_data.IVP;
    Oprt    =preproc_data.Oprt;
    shippar =preproc_data.shippar;
    bdyassim=preproc_data.bdyassim;
    
    if shippar.check==1
        set(handles.pushbutton_ship_setupGui,'Enable','off')
    end
          
    ODEsolve;
    
else
    clc
    disp('Please do pre-processing')
end

set(handles.pushbutton_run,'enable','off');
set(handles.togglebutton_logfile,'enable','on');
set(handles.pushbutton1,'enable','on');
set(handles.togglebutton_logfile,'value',0);


% --- Executes on button press in WT_Other_data.
function WT_Other_data_Callback(hObject, eventdata, handles)
% hObject    handle to WT_Other_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Input_gui = guidata(hObject);

Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
pathnow=Curdir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load Data file');

[statusbarObj]=JavaFrame_handling();

if directory~=0
    statusbarObj.setText('loading data');
    temp=load([directory,file_name]);
    if ~isdeployed
        clc;
    end
    
    
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_influx=temp.(namevar{1});
    else
        my_influx=temp;
    end
    clearvars temp;
    
    if max(max(isnan(my_influx)))==0
        if isstruct(my_influx)==0
            t_init=my_influx(1,1);t_end=my_influx(end,1);
            dt=my_influx(2,1)-t_init;
            %         Xinflux=my_influx(1,2);
            %         set(handles.edit_Xinflux,'Userdata',Xinflux);
            %         set(handles.edit_Xinflux,'String',num2str(Xinflux,3));
            
            set(handles.edit_time,'Userdata',[t_init t_end]);
            set(handles.edit_time,'string',[num2str(roundn(t_init,-2)),';', num2str(roundn(t_end,-2))]);
            set(handles.edit_timestep,'userdata',dt);
            set(handles.edit_timestep,'string',num2str(roundn(dt,-2)))
            set(handles.edit_timestep,'string',num2str(roundn(dt,-2)))
            
            set(Input_gui.WT_Other_data,'Userdata',my_influx);
            
            statusbarObj.setText([file_name,' has been loaded']);
        else
            statusbarObj.setText('wrong input file!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        end
    else
        statusbarObj.setText('Data is not loaded. There is NaN value!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    end
else
    statusbarObj.setText('No data loaded');
end

% --- Executes on button press in pushbutton_load_file_IVP.
function pushbutton_load_file_IVP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_file_IVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
pathnow=Curdir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load initial condition');

[statusbarObj]=JavaFrame_handling();

if directory~=0
    
    statusbarObj.setText('loading data');
    temp=load([directory,file_name]);
    
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    assignin('base','my_data',my_data)
    if isstruct(my_data)==0
        if length(my_data(:,1))<=2
            statusbarObj.setText('Wrong format the initial condition data');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        elseif length(my_data(1,:))==2 || length(my_data(1,:))==3 || length(my_data(1,:))==4
            IVPdat.data=my_data;
            IVPdat.name=file_name;
            set(handles.pushbutton_load_file_IVP,'Userdata',IVPdat); %record function name in handle
            statusbarObj.setText([file_name,' has been loaded']);
        else
            statusbarObj.setText('Wrong format the initial condition data');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        end
    else
        statusbarObj.setText('Wrong the initial condition data file');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    end
else
    statusbarObj.setText('No data loaded');
end


% --- Executes on selection change in dynModel.
function dynModel_Callback(hObject, eventdata, handles)
% hObject    handle to dynModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dynModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dynModel
contents = cellstr(get(hObject,'String'));
dynmodel=contents{get(hObject,'Value')};

global IDdemo
if IDdemo==1
    ID=get(hObject,'Value');
    [statusbarObj]=JavaFrame_handling();
    if ID==3||ID==4
        set(handles.dynModel,'value',2);
        statusbarObj.setText('HS3 or HS4 is not available in the demo version');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
    else
        statusbarObj.setText('');
    end
end

if strcmp(dynmodel,'HS1')||strcmp(dynmodel,'SWE1')
    set(handles.edit_Nonlin_Adj,'enable','off')
else
    set(handles.edit_Nonlin_Adj,'enable','on')
end
% --- Executes during object creation, after setting all properties.
function dynModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dynModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Breaking.
function checkbox_Breaking_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Breaking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Breaking
CheckBreak=get(hObject,'Value');

checkbox_def_param=findobj('Tag','checkbox_def_paramB');
Break_param=findobj('Tag','Breaking_param');
Break_param2=findobj('Tag','edit_paramB_UF_T');

if CheckBreak==1
    set(Break_param,'Visible','on');
    set(checkbox_def_param,'Visible','on');
    set(checkbox_def_param,'Enable','on');
    set(handles.text72,'Visible','on');
    set(handles.Breaking_param_text,'Visible','on');
    
    if get(handles.checkbox_def_paramB,'Value')==1
        set(Break_param2,'Enable','off');
        set(Break_param2,'Visible','on');
    else
        set(Break_param2,'Enable','on');
        set(Break_param2,'Visible','on');
    end
    set(Break_param2,'Visible','on');
else
    set(Break_param,'Visible','off');
    set(checkbox_def_param,'Visible','off');
    set(handles.text72,'Visible','off');
    set(handles.Breaking_param_text,'Visible','off');
    if get(handles.checkbox_def_paramB,'Value')==1
        set(Break_param2,'Enable','off');
        set(Break_param2,'Visible','off');
    else
        set(Break_param2,'Enable','on');
        set(Break_param2,'Visible','on');
    end
    set(Break_param2,'Visible','off');
end

function Breaking_param_Callback(hObject, eventdata, handles)
% hObject    handle to Breaking_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Breaking_param as text
%        str2double(get(hObject,'String')) returns contents of Breaking_param as a double
param=str2double(get(hObject,'String'));
set(handles.Breaking_param,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function Breaking_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Breaking_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_def_paramB.
function checkbox_def_paramB_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_def_paramB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_def_paramB
CheckdefBreak=get(hObject,'Value');

Break_param2=findobj('Tag','edit_paramB_UF_T');
if CheckdefBreak==1
    set(Break_param2,'Enable','off');
else
    set(Break_param2,'Enable','on');
end
set(Break_param2,'Visible','on');
set(handles.Breaking_param_text,'visible','on');
set(handles.text72,'visible','on');



function edit_paramB_UF_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_paramB_UF_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_paramB_UF_T as text
%        str2double(get(hObject,'String')) returns contents of edit_paramB_UF_T as a double
param=str2num(get(hObject,'String'));
set(handles.edit_paramB_UF_T,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_paramB_UF_T);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_paramB_UF_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_paramB_UF_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu_dispersion.
function popupmenu_dispersion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_dispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_dispersion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_dispersion
contents = cellstr(get(hObject,'String'));
Disp=contents{get(hObject,'Value')};

global IDdemo
if IDdemo==1
    if get(hObject,'Value')~=1
        set(handles.popupmenu_dispersion,'Value',1);
        [statusbarObj]=JavaFrame_handling();
        statusbarObj.setText('Only exact dispersion is available in the demo version');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
    end
else
    Disp_O       =findobj('Tag','togglebutton_InputDispersion');
    Group_vel_O  =findobj('Tag','togglebutton_InputGroupVel');
    
    if strcmp(Disp,'User-defined')
        set(Disp_O,'Visible','on'); set(Group_vel_O,'Visible','on');
        if isempty(get(handles.togglebutton_InputDispersion,'userdata'))
            DispEq='sign(k).*sqrt(9.81*k.*tanh(d.*k))';
            UgEq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)';
            set(handles.togglebutton_InputDispersion,'userdata',DispEq);
            set(handles.togglebutton_InputGroupVel,'userdata',UgEq);
        end
    else
        set(Disp_O,'Visible','off');  set(Group_vel_O,'Visible','off');
    end
end
% --- Executes during object creation, after setting all properties.
function popupmenu_dispersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_dispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_loadfile_disp.
function pushbutton_loadfile_disp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadfile_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_loadfile_groupVel.
function pushbutton_loadfile_groupVel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadfile_groupVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in checkbox_Filter.
function checkbox_Filter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Filter
check_LFHF=get(hObject,'Value');


if check_LFHF==1
    set(handles.edit_LFHF ,'Visible','on');
    set(handles.edit_LFHF ,'Enable','on');
    set(handles.text71 ,'Visible','on');set(handles.text61 ,'Visible','on');
else
    set(handles.edit_LFHF ,'Visible','off');
    set(handles.edit_LFHF ,'Enable','off');
    set(handles.text71 ,'Visible','off');set(handles.text61 ,'Visible','off');
end

function edit_LFHF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LFHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LFHF as text
%        str2double(get(hObject,'String')) returns contents of edit_LFHF as a double
param=str2num(get(hObject,'String'));
set(handles.edit_LFHF,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_LFHF);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_LFHF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LFHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_influx_direction.
function popupmenu_influx_direction_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_influx_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_influx_direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_influx_direction
contents = cellstr(get(hObject,'String'));
contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu_influx_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_influx_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_influx.
function popupmenu_influx_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_influx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_influx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_influx
contents = cellstr(get(hObject,'String'));
contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu_influx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_influx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time as text
%        str2double(get(hObject,'String')) returns contents of edit_time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_time,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
    if param(1)>param(2)
        statusbarObj.setText('Error:t_start>t_end');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_time);
    else %if correct
        statusbarObj.setText('');
    end
end



% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_timestep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timestep as text
%        str2double(get(hObject,'String')) returns contents of edit_timestep as a double
param=str2num(get(hObject,'String'));
set(handles.edit_timestep,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_timestep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p as text
%        str2double(get(hObject,'String')) returns contents of edit_p as a double
param=str2num(get(hObject,'String'));
set(handles.edit_p,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
end

if param<0
    statusbarObj.setText('Error: p is negative.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
end




% --- Executes during object creation, after setting all properties.
function edit_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
    
    if param(1)>param(2)
        statusbarObj.setText('Error:x_start>x_end');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_time);
    else %if correct
        statusbarObj.setText('');
    end
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


% --- Executes on selection change in popupmenu_bathy.
function popupmenu_bathy_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bathy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bathy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bathy
contents = cellstr(get(hObject,'String'));
Bathy=contents{get(hObject,'Value')};

depth_F      =findobj('Tag','edit_depth_F');
depth_O      =findobj('Tag','pushbutton_load_file_bathy');
Slope     =findobj('Tag','edit_slope');

if strcmp(Bathy,'Flat')
    set(depth_F ,'tooltipString','Specify the depth');
    set(depth_F ,'Enable','on');
    set(depth_O ,'Visible','off');
    set(Slope ,'Enable','off');
    set(handles.popupmenu_depth_reff,'enable','off');
    set(handles.text112,'enable','off');
    set(handles.edit_Dmidref,'enable','off');
    set(handles.text118,'enable','off');
    set(handles.text119,'enable','off');
    
    
elseif strcmp(Bathy,'Slope')
    set(depth_F ,'tooltipString','Specify the depth. Input format: depth [at the left; at the right]');
    set(Slope,'tooltipString','Specify the slope parameters. Input format: gradient; position of the foot of slope');
    set(depth_F ,'Enable','on');
    set(depth_O ,'Visible','off');
    set(Slope ,'Enable','on');
    set(handles.popupmenu_depth_reff,'enable','on');
    set(handles.text112,'enable','on');
elseif strcmp(Bathy,'Slope (Run-Up)')
    set(depth_F ,'tooltipString','Specify a maximum depth and a minimum shoreline depth (ie: 2% Hs)');
    set(Slope,'tooltipString','Specify the slope parameters. Input format: gradient; shore position');
    set(depth_F ,'Enable','on');
    set(depth_O ,'Visible','off');
    set(Slope ,'Enable','on');
    set(handles.popupmenu_depth_reff,'enable','on','value',2);
    set(handles.text112,'enable','on');
else
    set(depth_F ,'Enable','off');
    set(depth_O,'Visible','on');
    set(Slope ,'Enable','off');
    set(depth_O,'tooltipString','Load bathymetry data')
    set(handles.popupmenu_depth_reff,'enable','on');
    set(handles.text112,'enable','on');
end
set(depth_O,'Enable','on');

% --- Executes during object creation, after setting all properties.
function popupmenu_bathy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bathy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depth_F_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth_F as text
%        str2double(get(hObject,'String')) returns contents of edit_depth_F as a double
param=str2num(get(hObject,'String'));
Id =get(handles.popupmenu_bathy,'Value');
[statusbarObj]=JavaFrame_handling();
if Id==2||Id==3
    if length(param)~=2
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
    else %if correct
        statusbarObj.setText('');
    end
elseif Id==1 || Id==4
    if length(param)~=1
        statusbarObj.setText('wrong input format');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
    else %if correct
        statusbarObj.setText('');
    end
end
set(handles.edit_depth_F,'Userdata',param);

% --- Executes during object creation, after setting all properties.
function edit_depth_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_load_file_bathy.
function pushbutton_load_file_bathy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_file_bathy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
pathnow=Curdir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load bathymetry');

[statusbarObj]=JavaFrame_handling();

if directory~=0
    
    statusbarObj.setText('loading data');
    temp=load([directory,file_name]);
    
    if isstruct(temp)
        namevar = fieldnames(temp);
        my_data=temp.(namevar{1});
    else
        my_data=temp;
    end
    clearvars temp;
    
    if isstruct(my_data)==0
        bath.filename=file_name;
        bath.data=my_data;
        
        if length(my_data(1,:))~=2 || length(my_data(:,1))<=2
            statusbarObj.setText('Wrong format the bathymetry data');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
        else
            set(handles.pushbutton_load_file_bathy,'Userdata',bath); %record function name in handle
            statusbarObj.setText([file_name,' has been loaded']);
            if any(my_data(:,2)>=0)
                statusbarObj.setText('Warning: there are positive values of the bathymetry. It indicates Run-up case.');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.blue);
                
                set(handles.togglebutton_friction,'Visible','on')
                set(handles.checkbox_bottom_friction,'value',1)
                set(handles.popupmenu_depth_reff,'value',2);
                set(handles.edit_depth_F,'enable','on');
                set(handles.edit_depth_F,'tooltipstring','Specify a minimum depth for shore line (ie: 2% Hs)');
            else
                set(handles.togglebutton_friction,'Visible','off')
                set(handles.checkbox_bottom_friction,'value',0)
                 set(handles.edit_depth_F,'enable','off');
            end
        end
        
    else
        statusbarObj.setText('Wrong data file');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    end
    
else
    statusbarObj.setText('No data loaded');
end



function edit_FBL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FBL as text
%        str2double(get(hObject,'String')) returns contents of edit_FBL as a double
param=str2num(get(hObject,'String'));
set(handles.edit_FBL,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('Specify length of FBL at the left and right bdy. [L;R]');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_FBL);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_FBL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_wall.
function checkbox_wall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_wall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_wall
check_wall= get(hObject,'Value');


if check_wall==1
    set(handles.edit_Xwall ,'Visible','on');
    set(handles.edit_Xwall,'tooltipString','Specify the wall position.')
    set(handles.popupmenu_WallType,'Visible','on');
    set(handles.popupmenu_WallType,'tooltipString','Choose a type of wall. uniform or freq. dependent');
    set(handles.popupmenu_wall_method,'Visible','on');
    set(handles.text_wall_method,'Visible','on');
    Id=get(handles.popupmenu_WallType,'value');
    set(handles.text86,'Visible','on');
    if Id==1
        set(handles.text84,'Visible','on');
        set(handles.edit_Refl_percent,'Visible','on');
        set(handles.edit_Refl_percent,'tooltipString','Specify the reflection coefficient in [0,1]');
        set(handles.togglebutton_inputEq,'Visible','off');
    else
        set(handles.edit_Refl_percent,'Visible','off');
        set(handles.text84,'Visible','off');
        set(handles.togglebutton_inputEq,'Visible','on');
        set(handles.togglebutton_inputEq,'tooltipString','Input formula the reflection coefficient. use f as frequency [rad/s] variable. ]');
    end
    set(handles.text68,'Visible','on');
    set(handles.text83,'Visible','on');
    
    if get(handles.dynModel,'value')>2
        set(handles.edit_Nonlin_Adj,'enable','on')
    end
else
    set(handles.popupmenu_wall_method,'Visible','off');
    set(handles.text_wall_method,'Visible','off');
    set(handles.edit_Xwall ,'Visible','off');
    set(handles.popupmenu_WallType,'Visible','off');
    set(handles.edit_Refl_percent,'Visible','off');
    set(handles.togglebutton_inputEq,'Visible','off');
    set(handles.text68,'Visible','off');
    set(handles.text83,'Visible','off');
    set(handles.text84,'Visible','off');
    set(handles.text86,'Visible','off');
end



function edit_note_Callback(hObject, eventdata, handles)
% hObject    handle to edit_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_note as text
%        str2double(get(hObject,'String')) returns contents of edit_note as a double
get(hObject,'String');

% --- Executes during object creation, after setting all properties.
function edit_note_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_note (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_savename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_savename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_savename as text
%        str2double(get(hObject,'String')) returns contents of edit_savename as a double
Name=get(hObject,'String');
Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
Path=Curdir{get(handles.popupmenu_CurDir,'Value')};
[statusbarObj]=JavaFrame_handling();
global FlagWarn
if exist([Path,'/Output/',Name],'dir')
    statusbarObj.setText('Warning: Project exists already, it will be overwritten!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.blue);
    uicontrol(handles.edit_savename);
    FlagWarn=1;
else
    statusbarObj.setText('');
    FlagWarn=0;
end

% --- Executes during object creation, after setting all properties.
function edit_savename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_savename (see GCBO)
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

GUIinput=state;
save([state.proj.workdir,'\abproj_',state.proj.name,'.mat'], 'GUIinput')

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



function GUImain=adjust_variables(GUImain)
%this function to adjust variables name to be same with variable in
%MainUser.m

if strcmp(GUImain.WT,'None')
    GUImain.WT=1;
elseif strcmp(GUImain.WT,'Harmonic')
    GUImain.WT=2;
elseif strcmp(GUImain.WT,'Jonswap')
    GUImain.WT=3;
elseif strcmp(GUImain.WT,'User-defined')
    GUImain.WT=4;
end


if strcmp(GUImain.evol,'HS1')
    GUImain.evol='HS';
    GUImain.nonlin=1;
elseif strcmp(GUImain.evol,'HS2')
    GUImain.evol='HS';
    GUImain.nonlin=2;
elseif strcmp(GUImain.evol,'HS3')
    GUImain.evol='HS';
    GUImain.nonlin=3;
elseif strcmp(GUImain.evol,'HS4')
    GUImain.evol='HS';
    GUImain.nonlin=4;
elseif strcmp(GUImain.evol,'SWE1')
    GUImain.evol='SWE';
    GUImain.nonlin=1;
elseif strcmp(GUImain.evol,'SWE2')
    GUImain.evol='SWE';
    GUImain.nonlin=2;
end

if GUImain.Break==0
    GUImain.Break='No';
else
    GUImain.Break='Yes';
end

if strcmp(GUImain.disp,'Exact')
    GUImain.disp='OmExact';
    GUImain.groupvel='UgExact';
    GUImain.detaC2='detaC2Exact';   % For runup
    GUImain.Cder='CderExact';   % For runup
    GUImain.dispEqInput='';
    GUImain.groupVelEqInput='';
elseif strcmp(GUImain.disp,'KdV')
    GUImain.disp='OmKdV';
    GUImain.groupvel='UgKdV';
    GUImain.detaC2='detaC2Exact';   % For runup
    GUImain.Cder='CderExact';   % For runup
    GUImain.dispEqInput='';
    GUImain.groupVelEqInput='';
elseif strcmp(GUImain.disp,'BBM')
    GUImain.disp='OmBBM';
    GUImain.groupvel='UgBBM';
    GUImain.detaC2='detaC2Exact';   % For runup
    GUImain.Cder='CderExact';   % For runup
    GUImain.dispEqInput='';
    GUImain.groupVelEqInput='';
elseif strcmp(GUImain.disp,'SWE')
    GUImain.disp='OmSWE';
    GUImain.groupvel='UgSWE';
    GUImain.detaC2='detaC2Exact';   % For runup
    GUImain.Cder='CderExact';   % For runup
    GUImain.dispEqInput='';
    GUImain.groupVelEqInput='';
else
    GUImain.disp='OmUser';
    GUImain.groupvel='UgUser';
    GUImain.detaC2='detaC2User';
    GUImain.dispEqInput=GUImain.disp_O_file;
    GUImain.groupVelEqInput=GUImain.groupVel_O_file;
end

if strcmp(GUImain.inftype,'Area')
    GUImain.inftype='GVArea';
elseif strcmp(GUImain.inftype,'Area-Short')
    GUImain.inftype='GVAreaShort';
end

if strcmp(GUImain.bathy,'Flat')
    GUImain.bathy='F';
elseif strcmp(GUImain.bathy,'Slope')
    GUImain.bathy='B';
elseif strcmp(GUImain.bathy,'Slope (Run-Up)')
    GUImain.bathy='BR';
    GUImain.evol='HS';
    %     if GUImain.nonlin>2
    %         GUImain.evol='HS';
    %     else
    %         GUImain.evol='HS2';
    %     end
else
   
    if max(any(GUImain.bathyO.data(:,2)>=0))==1
        GUImain.bathy='UR';
        GUImain.evol='HS';
        %         if GUImain.nonlin>2
        %             GUImain.evol='HS';
        %         else
        %             GUImain.evol='HS2';
        %         end
    else
        GUImain.bathy='U';
    end
end

if GUImain.Wall==0
    GUImain.Wall='No';
else
    GUImain.Wall='Yes';
end




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
delete(findall(0, 'type', 'figure'));

% clc;
% close all;


function keep_guis
fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
for fh = fig_h
    uih = findobj( fh, 'Type', 'uicontrol' );
    if isempty( uih )
        delete( fh );
    end
end

function Initialization(handles)

global FlagWarn
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');

set(handles.pushbutton_run,'enable','off');
set(handles.togglebutton_logfile,'enable','off');

keep_guis; %% close all figures except gui

GUImain.proj.path=handles.input.pathnow;
GUImain.proj.module =handles.input.module;
GUImain.proj.projhist=handles.input.projecthist;
Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
GUImain.proj.projdir=Curdir{get(handles.popupmenu_CurDir,'Value')};
GUImain.proj.name=get(handles.edit_savename,'String');
GUImain.proj.workdir=[GUImain.proj.projdir,'\',GUImain.proj.name,'\'];
GUImain.proj.note=get(handles.edit_note,'String');

GUImain.IVP.type=get(handles.popupmenu_IVP,'value');
IVPname         =cellstr(get(handles.popupmenu_IVP,'String'));
GUImain.IVP.typename=IVPname{get(handles.popupmenu_IVP,'value')};
GUImain.IVP.file=get(handles.pushbutton_load_file_IVP,'Userdata');
GUImain.IVP.A=get(handles.edit_IVP_A,'Userdata');
GUImain.IVP.lambda=get(handles.edit_IVP_lambda,'Userdata');
GUImain.IVP.x0=get(handles.edit_IVP_x0,'Userdata');

wavetype=cellstr(get(handles.WaveType,'String'));
GUImain.WT=wavetype{get(handles.WaveType,'Value')};

if get(handles.WaveType,'Value')==2
    GUImain.Hs= get(handles.edit_Hs,'Userdata')*2; % harmonic input is amplitude
elseif get(handles.WaveType,'Value')==3
    GUImain.Hs= get(handles.edit_Hs,'Userdata');
end

GUImain.Tp= get(handles.edit_Tp,'Userdata');
GUImain.Jsg= get(handles.edit_Jsg,'Userdata');
GUImain.WT_O_data=get(handles.WT_Other_data,'Userdata');

GUImain.Filter= get(handles.checkbox_Filter,'Value');
GUImain.Filter_LFHF=get(handles.edit_LFHF,'Userdata');
GUImain.ramp.check=get(handles.checkbox_ramp,'Value');
GUImain.ramp.val=get(handles.edit_ramp,'userdata');


modeltype=cellstr(get(handles.dynModel,'String'));
GUImain.evol=modeltype{get(handles.dynModel,'Value')};
GUImain.Break= get(handles.checkbox_Breaking,'Value');
GUImain.Break_param= get(handles.Breaking_param,'Userdata');
GUImain.BreakDef= get(handles.checkbox_def_paramB,'Value');
GUImain.Break_paramDef= get(handles.edit_paramB_UF_T,'Userdata');

disptype=cellstr(get(handles.popupmenu_dispersion,'String'));
GUImain.disp=disptype{get(handles.popupmenu_dispersion,'Value')};
GUImain.disp_O_file=get(handles.togglebutton_InputDispersion,'Userdata');
GUImain.groupVel_O_file=get(handles.togglebutton_InputGroupVel,'Userdata');
Influxtype=cellstr(get(handles.popupmenu_influx,'String'));
GUImain.inftype=Influxtype{get(handles.popupmenu_influx,'Value')};
Influxdir=cellstr(get(handles.popupmenu_influx_direction,'String'));
GUImain.infdir=Influxdir{get(handles.popupmenu_influx_direction,'Value')};
GUImain.Assim.check=get(handles.checkbox17,'value');
GUImain.Assim.data=get(handles.togglebutton_bdy_assim,'userdata');

GUImain.tinterv=get(handles.edit_time,'Userdata');
GUImain.tstep=get(handles.edit_timestep,'Userdata');
GUImain.xinterv=get(handles.edit_xinterv,'Userdata');
GUImain.Xinflux=get(handles.edit_Xinflux,'Userdata');
GUImain.dx=get(handles.edit_p,'Userdata');
GUImain.cutfracK=get(handles.edit_cutfracK,'Userdata');
GUImain.FBL=get(handles.edit_FBL,'Userdata');
GUImain.NonlinAdj=get(handles.edit_Nonlin_Adj,'Userdata');
bathytype=cellstr(get(handles.popupmenu_bathy,'String'));
GUImain.bathy=bathytype{get(handles.popupmenu_bathy,'Value')};

GUImain.depth=get(handles.edit_depth_F,'Userdata');
GUImain.slope=get(handles.edit_slope,'Userdata');
GUImain.bathyO=get(handles.pushbutton_load_file_bathy,'Userdata');
GUImain.bathyInterp=get(handles.popupmenu_depth_reff,'Value');
GUImain.DmidRef=get(handles.edit_Dmidref,'userdata');

GUImain.bottomfriction_check=get(handles.checkbox_bottom_friction,'Value');
GUImain.frictiondata=get(handles.togglebutton_friction,'Userdata');

GUImain.wind.check=get(handles.checkbox_wind,'value');
GUImain.wind.coef=get(handles.edit_coef_wind,'userdata');
GUImain.wind.tinterv=get(handles.edit_wind_tinterv,'userdata');
GUImain.wind.xinterv=get(handles.edit_wind_xinterv,'userdata');

GUImain.Wall= get(handles.checkbox_wall,'Value');
GUImain.Wall_position= get(handles.edit_Xwall,'Userdata');
GUImain.Wall_Coef= get(handles.edit_Refl_percent,'Userdata');
GUImain.Wall_method=get(handles.popupmenu_wall_method,'Value');
GUImain.Wall_type= get(handles.popupmenu_WallType,'Value');
GUImain.Wall_reflEq=get(handles.togglebutton_inputEq,'Userdata');
GUImain.InteriorProp=get(handles.checkbox_interior,'Value');
GUImain.IP_time=get(handles.edit_Interior_Time,'Userdata');
GUImain.IP_dt=get(handles.edit_dtIP,'Userdata');
GUImain.partition=get(handles.edit_partition,'userdata');
GUImain.partitionId=get(handles.checkbox_partition,'value');
GUImain.partitionCombines=get(handles.checkbox_combineFiles,'Value');
GUImain.shippar=get(handles.pushbutton_ship_setupGui,'userdata');
GUImain.shippar.check=get(handles.checkbox_ship,'value');
%clc;
if isempty(GUImain.proj.projdir) || strcmpi(GUImain.proj.projdir,'--')
    statusbarObj.setText(' Specify a project directory!')
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.popupmenu_CurDir);
    return;
end

if isempty(GUImain.proj.name)
    statusbarObj.setText('Specify a Project name!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_savename);
    return
end


if GUImain.Break==1
    if strcmp(GUImain.evol,'HS1')||strcmp(GUImain.evol,'SWE1') ||strcmp(GUImain.evol,'SWE2')
        statusbarObj.setText('Use a nonlinear model (HS2,3,4)!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.dynModel);
        return;
    else
        statusbarObj.setText('');
    end
    
    if isempty(GUImain.Break_param)
        statusbarObj.setText('Specify an initiation breaking coefficient');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.dynModel);
        return;
    end
    
    if length(GUImain.Break_param)~=1
        statusbarObj.setText('Input only a value of the initiation breaking coefficient');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.Breaking_param);
        return;
    end
    
    if GUImain.Break_param<0
        statusbarObj.setText('Error: negative input of the breaking initiation parameter');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.Breaking_param);
        return;
    end
    
    if GUImain.Break_param==0
        statusbarObj.setText('Error: zero input of the breaking initiation parameter');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.Breaking_param);
        return;
    end
    
    if GUImain.BreakDef==0
        if isempty(GUImain.Break_paramDef)
            statusbarObj.setText('Specify the termination breaking parameters');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_paramB_UF_T);
            return;
        end
        
        if length(GUImain.Break_paramDef)~=2
            statusbarObj.setText('Erros: Wrong input of the termination breaking parameters');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_paramB_UF_T);
            return;
        end
        
        
        if any(GUImain.Break_paramDef<0)
            statusbarObj.setText('Error: negative input of the breaking termination parameter');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_paramB_UF_T);
            return;
        end
        
        if any(GUImain.Break_paramDef==0)
            statusbarObj.setText('Error: zero input of the breaking termination parameter');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_paramB_UF_T);
            return;
        end
        
        
    end
end

if strcmp(GUImain.disp,'User-defined')
    if isempty(GUImain.disp_O_file)
        statusbarObj.setText('Input a dispersion formula. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputDispersion);
        return;
    end
    if isempty(strfind(GUImain.disp_O_file,'k'))
        statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputDispersion);
        return;
    end
    if isempty(strfind(GUImain.disp_O_file,'d'))
        statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputDispersion);
        return;
    end
    if isempty(GUImain.groupVel_O_file)
        statusbarObj.setText('Input a group velocity formula. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputGroupVel);
        return;
    end
    if isempty(strfind(GUImain.groupVel_O_file,'k'))
        statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputDispersion);
        return;
    end
    if isempty(strfind(GUImain.groupVel_O_file,'d'))
        statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_InputDispersion);
        return;
    end
    
end


if GUImain.IVP.type==2||GUImain.IVP.type==3
    if isempty(GUImain.IVP.A)
        statusbarObj.setText('Specify the amplitude of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_A);
        return;
    end
    
    if length(GUImain.IVP.A)~=1
        statusbarObj.setText('Input only a value of the amplitude of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_A);
        return;
    end
    
    if isempty(GUImain.IVP.lambda)
        statusbarObj.setText('Specify the standard deviation of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_lambda);
        return;
    end
    if length(GUImain.IVP.lambda)~=1
        statusbarObj.setText('Input only a value of the standard deviation of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_lambda);
        return;
    end
    
    if isempty(GUImain.IVP.x0)
        statusbarObj.setText('Specify the center position');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_x0);
        return;
    end
    
    if length(GUImain.IVP.x0)~=1
        statusbarObj.setText('Input only a value of the center position');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_x0);
        return;
    end
    
end

if GUImain.IVP.type==4
    
    if isfield(GUImain.IVP.file,'data')
        if isempty(GUImain.IVP.file)
            statusbarObj.setText('Load an input data of the initial condition');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_IVP);
            return;
        end
    else %adjust from previous setting
        statusbarObj.setText('Load an input data of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_IVP);
        return;
    end
    
    if length(GUImain.IVP.file.data(1,:))==2 || length(GUImain.IVP.file.data(1,:))==3 || length(GUImain.IVP.file.data(1,:))==4
        if GUImain.shippar.check==1 && length(GUImain.IVP.file.data(1,:))~=4
            statusbarObj.setText('Wrong data of the initial condition');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_IVP);
            return;
        elseif GUImain.shippar.check==0 && length(GUImain.IVP.file.data(1,:))==4
            statusbarObj.setText('Wrong data of the initial condition');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_IVP);
            return;
        end
        
    else
        statusbarObj.setText('Wrong data of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_IVP);
        return;
    end
    
    if length(GUImain.IVP.file.data(:,1))<=2
        statusbarObj.setText('Wrong data of the initial condition');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_IVP);
        return;
    end
    
end

if strcmp(GUImain.WT,'Harmonic') || strcmp(GUImain.WT,'Jonswap')
    if isempty(GUImain.Hs)
        
        if strcmp(GUImain.WT,'Harmonic')
            statusbarObj.setText('Input Ampitude');
        else
            statusbarObj.setText('Input Hs');
        end
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Hs);
        return;
    end
    
    if length(GUImain.Hs)~=1
        if strcmp(GUImain.WT,'Harmonic')
            statusbarObj.setText('Input only a value of Ampitude');
        else
            statusbarObj.setText('Input only a value of Hs');
        end
        
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Hs);
        return;
    end
    
    if GUImain.Hs==0
        if strcmp(GUImain.WT,'Harmonic')
            statusbarObj.setText('Error: A is zero');
        else
            statusbarObj.setText('Error: Hs is zero');
        end
        
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Hs);
        return;
    end
    
    if GUImain.Hs<0
        if strcmp(GUImain.WT,'Harmonic')
            statusbarObj.setText('Error: A is negative');
        else
            statusbarObj.setText('Error: Hs is negative');
        end
        
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Hs);
        return;
    end
    
    if isempty(GUImain.Tp)
        statusbarObj.setText('Input Tp');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Tp);
        return;
    end
    
    if length(GUImain.Tp)~=1
        statusbarObj.setText('Input only a value of Tp');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Tp);
        return;
    end
    
    if GUImain.Tp<0
        statusbarObj.setText('Error: Tp is negative');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Tp);
        return;
    end
    
    if GUImain.Tp==0
        statusbarObj.setText('Error: Tp is zero');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Tp);
        return;
    end
end

if strcmp(GUImain.WT,'Jonswap')
    if isempty(GUImain.Jsg)
        statusbarObj.setText('Input gamma');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Jsg);
        return;
    end
    
    if length(GUImain.Jsg)~=1
        statusbarObj.setText('Input only a value of gamma');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Jsg);
        return;
    end
    
    if GUImain.Jsg<0
        statusbarObj.setText('Error: gamma is negative');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Jsg);
        return;
    end
    
    if GUImain.Jsg==0
        statusbarObj.setText('Error: gamma is zero');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Jsg);
        return;
    end
    
elseif strcmp(GUImain.WT,'User-defined')
    if isempty(GUImain.WT_O_data)
        statusbarObj.setText('Load an influx data');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
        uicontrol(handles.WT_Other_data);
        return;
    end
    
end

if GUImain.Filter==1
    if isempty(GUImain.Filter_LFHF)
        statusbarObj.setText('Input the freq. interval[LF;HF] [rad/s]!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_LFHF);
        return;
    end
    
    if length(GUImain.Filter_LFHF)~=2
        statusbarObj.setText('Error: wrong input format of the freq. interval.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_LFHF);
        return;
    end
    
    if GUImain.Filter_LFHF(1)>GUImain.Filter_LFHF(2)
        statusbarObj.setText('Error: wrong input format of the freq. interval.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_LFHF);
        return;
    end
end

if GUImain.ramp.check==1
    if isempty(GUImain.ramp.val)
        statusbarObj.setText('Input length of the ramp function.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_ramp);
        return;
    end
    
    if length(GUImain.ramp.val)~=1
        statusbarObj.setText('Input only a value of the length of the ramp .');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_ramp);
        return;
    end
    
    if  GUImain.ramp.val==0
        statusbarObj.setText('Error: the length of the ramp is zero.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_ramp);
        return;
    end
    
    if  GUImain.ramp.val<0
        statusbarObj.setText('Error: the length of the ramp is negative.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_ramp);
        return;
    end
    
end
if isempty(GUImain.tinterv)
    statusbarObj.setText('Input a time interval [s]!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
    return;
end
if  length(GUImain.tinterv)~=2
    statusbarObj.setText('Wrong input format of the time interval!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
    return;
end

if GUImain.tinterv(1)>GUImain.tinterv(2)
    statusbarObj.setText('Error:t_start>t_end');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
    return
end

if isempty(GUImain.tstep)
    statusbarObj.setText('Input a time step [s]!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_timestep);
    return;
end

if length(GUImain.tstep)~=1
    statusbarObj.setText('Input only a value of the time step!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_timestep);
    return;
end

if GUImain.tstep<0
    statusbarObj.setText('Error: time step is negative');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_timestep);
    return;
end

if GUImain.tstep==0
    statusbarObj.setText('Error: time step is zero');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_timestep);
    return;
end


if strcmp(GUImain.WT,'User-defined')
    if GUImain.tinterv(1)<floor(GUImain.WT_O_data(1,1));
        statusbarObj.setText(['Error:t_start< minimum time (t=',num2str(floor(GUImain.WT_O_data(1,1))),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_time);
        return;
    end
    
    dt_data=GUImain.WT_O_data(2,1)-GUImain.WT_O_data(1,1);
    if  GUImain.tstep<floor(dt_data);
        statusbarObj.setText('dt is too small, please increase!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
        uicontrol(handles.edit_timestep);
        return;
    end
end




if GUImain.Assim.check==1
    if isempty(GUImain.Assim.data)
        statusbarObj.setText('Specify assimilation data');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_bdy_assim);
        return;
    end
    if isfield(GUImain.Assim.data,'eta')
        if  isempty(GUImain.Assim.data.eta)
            statusbarObj.setText('Specify assimilation elevation data!')
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_bdy_assim);
            return;
        end
    end
    
    if isfield(GUImain.Assim.data,'velocity')
        if GUImain.Assim.data.checkbox_velocity==1
            if isempty(GUImain.Assim.data.velocity)
                statusbarObj.setText('Specify assimilation velocity data!')
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.togglebutton_bdy_assim);
                return;
            end
        end
    end
    
    if isempty(GUImain.Assim.data.xinterv)|| length(GUImain.Assim.data.xinterv)~=2
        statusbarObj.setText('Specify x interval in assimilation setup!')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_bdy_assim);
        return;
    end
    
    if isempty(GUImain.Assim.data.smoothfact)|| length(GUImain.Assim.data.smoothfact)~=1
        statusbarObj.setText('Specify smooth factor in assimilation setup!')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_bdy_assim);
        return;
    end
    
    if isempty(GUImain.Assim.data.tinterv)|| length(GUImain.Assim.data.tinterv)~=2
        statusbarObj.setText('Specify t interval in assimilation setup!')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_bdy_assim);
        return;
    end
    
    if isempty(GUImain.Assim.data.tstep) || length(GUImain.Assim.data.tstep)~=1
        statusbarObj.setText('Specify time step in assimilation setup!')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_bdy_assim);
        return;
    end
    
    if GUImain.Assim.data.checkbox_nonlin_adj==1
        if isempty(GUImain.Assim.data.nonlinadj_length)|| length(GUImain.Assim.data.nonlinadj_length)~=1
            statusbarObj.setText('Specify length of nonlinear adjustment in assimilation setup!')
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_bdy_assim);
            return;
        end
    end
    
end

if GUImain.IVP.type==1 && strcmp(GUImain.WT,'None') && GUImain.Assim.check==0
    statusbarObj.setText('There is no generation. Wave input and initial consdition are zero');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if GUImain.wind.check==1
    if length(GUImain.wind.coef)~=1
        statusbarObj.setText('Specify a wind term coef.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_coef_wind);
        return;
    end
    
    if length(GUImain.wind.tinterv)~=2
        statusbarObj.setText('Specify a time interval [tstart; tend] for wind forcing.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_wind_tinterv);
        return;
    end
    
    if length(GUImain.wind.xinterv)~=2
        statusbarObj.setText('Specify a x interval [Xstart; Xend] for wind forcing.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_wind_xinterv);
        return;
    end
    
end

if isempty(GUImain.xinterv)
    statusbarObj.setText('Input the spatial interval');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xinterv);
    return;
end

if length((GUImain.xinterv))~=2
    statusbarObj.setText('Wrong input format of the spatial interval');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xinterv);
    return;
end

if GUImain.xinterv(1)>GUImain.xinterv(2)
    statusbarObj.setText('Error:x_start>x_end');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_xinterv);
    return
end

if isempty(GUImain.Xinflux)
    statusbarObj.setText('Input the influx position');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Xinflux);
    return;
end

if  length((GUImain.Xinflux))~=1
    statusbarObj.setText('Input only a value of the influx position');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Xinflux);
    return;
end

if GUImain.Xinflux<GUImain.xinterv(1)||GUImain.Xinflux>GUImain.xinterv(2)
    statusbarObj.setText('Error: Influx position is out of the simulation domain');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Xinflux);
    return;
end

if ~isempty(GUImain.IVP.x0)
    if GUImain.IVP.x0<GUImain.xinterv(1) || GUImain.IVP.x0>GUImain.xinterv(2)
        statusbarObj.setText('Error: the center position of initial condition is out of the simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_IVP_x0);
        return;
    end
end

if GUImain.IVP.type==4
    if isfield(GUImain.IVP.file,'data')
        if (GUImain.IVP.file.data(1,1)<GUImain.xinterv(1) && ...
                GUImain.IVP.file.data(end,1)<GUImain.xinterv(1))
            statusbarObj.setText('Initial condition data is of the simulation domain');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_IVP);
            return;
        end
        
        if (GUImain.IVP.file.data(1,1)>GUImain.xinterv(2) && ...
                GUImain.IVP.file.data(end,1)>GUImain.xinterv(2))
            statusbarObj.setText('Initial condition data is out of the simulation domain');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_IVP);
            return;
        end
    end
end

if isempty(GUImain.dx)
    statusbarObj.setText('Input grid size dx');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_p);
    return;
end

if length(GUImain.dx)~=1
    statusbarObj.setText('Input only a value of p');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_p);
    return;
end

if GUImain.dx<=0
    statusbarObj.setText('Error: dx <= 0');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_p);
    return;
end

Nx=round((GUImain.xinterv(2)-GUImain.xinterv(1))*10000/...
    (GUImain.dx*10000))+1;

if any(factor(Nx)>13)
    statusbarObj.setText(...
        '>> Number of points in x axis has a largest prime factor>13. Please adjust dx!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_p);
    return;
end

if isempty(GUImain.cutfracK)
    statusbarObj.setText('Input a cutfraction wave number');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
    return;
end

if  length(GUImain.cutfracK)~=1
    statusbarObj.setText('Input only a value of the cutfraction wave number');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
    return;
end

if fix(GUImain.cutfracK)~=GUImain.cutfracK
    statusbarObj.setText('Error: the Cutfrac k is not an integer');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
    return;
end

if  GUImain.cutfracK<0
    statusbarObj.setText('Error: the Cutfrac k is negative');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
    return;
end

if  GUImain.cutfracK==0
    statusbarObj.setText('Error: the Cutfrac k is zero');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
    return;
end

if strcmp(GUImain.bathy,'Flat')
    if isempty(GUImain.depth)
        statusbarObj.setText('Input a depth!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if length(GUImain.depth)~=1
        statusbarObj.setText('Input only a value of the depth!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if GUImain.depth<0
        statusbarObj.setText('Error: depth is negative');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if GUImain.depth==0
        statusbarObj.setText('Error: depth is zero');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    set(handles.edit_depth_F ,'tooltipString','Input format: depth');
elseif strcmp(GUImain.bathy,'Slope')
    if isempty(GUImain.depth)
        statusbarObj.setText('Input depth [at the left; at the right]!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    if length(GUImain.depth)~=2
        statusbarObj.setText('Wrong input format of the depth!')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if  GUImain.depth(1)<0
        statusbarObj.setText('Error: depth is negative');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    if  GUImain.depth(2)<0
        statusbarObj.setText('Error: depth is negative');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if any(GUImain.depth==0)
        statusbarObj.setText('Error: depth is zero');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if isempty(GUImain.slope)||length(GUImain.slope)~=2
        statusbarObj.setText('Input [slope; slope start position]!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    
    if GUImain.slope(1)<0
        statusbarObj.setText('Error: slope is negative!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    if GUImain.slope(1)==0
        statusbarObj.setText('Error: slope is zero!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    if  GUImain.slope(2)<GUImain.xinterv(1)||GUImain.slope(2)>GUImain.xinterv(2)
        statusbarObj.setText('Error: slope start position out of simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    if  GUImain.bathyInterp==2
        if length(GUImain.DmidRef)~=1||isempty(GUImain.DmidRef)
            statusbarObj.setText('Specify the middle reference depth');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Dmidref);
            return;
        end
        if GUImain.DmidRef>=max(GUImain.depth(2),GUImain.depth(1))
            statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Dmidref);
            return;
        end
        if GUImain.DmidRef<=min(GUImain.depth(2),GUImain.depth(1))
            statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Dmidref);
            return;
        end
    end
    
    set(handles.edit_depth_F ,'tooltipString','Input format: depth [at the left; at the right]');
    set(handles.edit_slope,'tooltipString','Input format: slope; position of the foot of slope');
    
elseif strcmp(GUImain.bathy,'Slope (Run-Up)')
    
   
    if isempty(GUImain.depth)
        statusbarObj.setText('Input a maximum and minimum depth!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if length(GUImain.depth)~=2
        statusbarObj.setText('Wrong input format of the depth!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if GUImain.depth<0
        statusbarObj.setText('Error: depth is negative!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_depth_F);
        return;
    end
    
    if isempty(GUImain.slope)||length(GUImain.slope)~=2
        statusbarObj.setText('Input [slope; shore position]!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    if GUImain.slope(1)<0
        statusbarObj.setText('Error: slope is negative!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    if  GUImain.slope(2)<GUImain.xinterv(1)||GUImain.slope(2)>GUImain.xinterv(2)
        statusbarObj.setText('Error: shore position out of simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_slope);
        return;
    end
    
    
    if strcmp(GUImain.evol,'HS1')||strcmp(GUImain.evol,'HS4')|| ...
            strcmp(GUImain.evol,'SWE1')%||strcmp(GUImain.evol,'HS3')
        
        statusbarObj.setText('Use HS2 or HS3 dynamic model for the run-up case!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.dynModel);
        return;
    end
    
    
    if strcmp(GUImain.disp,'User-defined')
        statusbarObj.setText('User-defined dispersion is not applicable in the run-up case!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.popupmenu_dispersion);
        return;
    end
    
    if  GUImain.bathyInterp==2
        if length(GUImain.DmidRef)~=1||isempty(GUImain.DmidRef)
            statusbarObj.setText('Specify the middle reference depth');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Dmidref);
            return;
        end
        
        if GUImain.DmidRef>max(GUImain.depth)
            statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Dmidref);
            return;
        end
    end
elseif strcmp(GUImain.bathy,'User-defined')
    if isfield(GUImain.bathyO,'filename')   % Adjust from previous setting
        if isempty(GUImain.bathyO)
            statusbarObj.setText('Input a bathymetry file!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_load_file_bathy);
            return;
        end
    else
        statusbarObj.setText('Input a bathymetry file!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_bathy);
        return;
    end
    
    
    if length(GUImain.bathyO.data(1,:))~=2 || length(GUImain.bathyO.data(:,1))<=2
        statusbarObj.setText('Error: Wrong input bathymetry!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_bathy);
        return;
    end
    
    
    
    if (GUImain.bathyO.data(1,1)<GUImain.xinterv(1) && ...
            GUImain.bathyO.data(end,1)<GUImain.xinterv(1))
        statusbarObj.setText('Bathymetry data is out of the simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_bathy);
        return;
    end
    
    if (GUImain.bathyO.data(1,1)>GUImain.xinterv(2) && ...
            GUImain.bathyO.data(end,1)>GUImain.xinterv(2))
        statusbarObj.setText('Bathymetry data is out of the simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton_load_file_bathy);
        return;
    end
    
    
    if any(GUImain.bathyO.data(:,2)>=0)
        if strcmp(GUImain.evol,'HS1')||strcmp(GUImain.evol,'HS4')||   ...
                strcmp(GUImain.evol,'SWE1')%||strcmp(GUImain.evol,'HS3') ...
            statusbarObj.setText('Use HS2 or HS3 dynamic model for the run-up case!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.dynModel);
            return;
        end
        
        if strcmp(GUImain.disp,'User-defined')
            statusbarObj.setText('User-defined dispersion is not applicable in the run-up case!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.popupmenu_dispersion);
            return;
        end
        
        if  GUImain.bathyInterp==2
            if length(GUImain.DmidRef)~=1||isempty(GUImain.DmidRef)
                statusbarObj.setText('Specify the middle reference depth');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_Dmidref);
                return;
            end
            
            if GUImain.DmidRef>max(-GUImain.bathyO.data(:,2))
                statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_Dmidref);
                return;
            end
            
            
            if isempty(GUImain.depth)
                statusbarObj.setText('Input a minimum depth for shoreline!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_depth_F);
                return;
            end
            
            if length(GUImain.depth)~=1
                statusbarObj.setText('Wrong input format of the depth!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_depth_F);
                return;
            end
            
            if GUImain.depth<0
                statusbarObj.setText('Error: depth is negative!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_depth_F);
                return;
            end
            
    
        end
    else
        if  GUImain.bathyInterp==2
            if length(GUImain.DmidRef)~=1||isempty(GUImain.DmidRef)
                statusbarObj.setText('Specify the middle reference depth');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_Dmidref);
                return;
            end
            
            if GUImain.DmidRef>=max(-GUImain.bathyO.data(:,2))
                statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_Dmidref);
                return;
            end
            if GUImain.DmidRef<=min(-GUImain.bathyO.data(:,2))
                statusbarObj.setText('The middle reference depth must be in between max(D) and min(D)');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit_Dmidref);
                return;
            end
        end
    end
    
    
end

if GUImain.bottomfriction_check==1
    if isempty(GUImain.frictiondata)
        statusbarObj.setText('Specify the bottom friction coefficients!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.togglebutton_friction);
        return
    else
        if isempty(GUImain.frictiondata.cf)
            statusbarObj.setText('Specify the bottom friction coefficients!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_friction);
            return
        end
        if isempty(GUImain.frictiondata.interval)
            statusbarObj.setText('Specify the bottom friction intervals!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_friction);
            return
        end
    end
end


if GUImain.Wall==1
    if strcmp(GUImain.evol,'SWE1')||strcmp(GUImain.evol,'SWE2')
        statusbarObj.setText('Use only the HS dynamic model for embedded wall!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.dynModel);
        return;
    end
    
    if isempty(GUImain.Wall_position)
        statusbarObj.setText('Input a wall position!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Xwall);
        return
    end
    
    if length(GUImain.Wall_position)~=1
        statusbarObj.setText('Input only a value wall position!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Xwall);
        return
    end
    
    if GUImain.Wall_position<GUImain.xinterv(1)||GUImain.Wall_position>GUImain.xinterv(2)
        statusbarObj.setText('Error: wall position out of the simulation domain');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Xwall);
        return
    end
    
    if GUImain.Wall_type==1
        if isempty(GUImain.Wall_Coef)
            statusbarObj.setText('Specify the reflection coefficient in [0,1]!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Refl_percent);
            return
        end
        if length(GUImain.Wall_Coef)~=1
            statusbarObj.setText('Input only a value of the reflection coefficient!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Refl_percent);
            return
        end
        
        if GUImain.Wall_Coef<0
            statusbarObj.setText('Error: the reflection coefficient is negative!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Refl_percent);
            return
        end
        if GUImain.Wall_Coef>1
            statusbarObj.setText('Error: the reflection coefficient is larger than one!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Refl_percent);
            return
        end
    else
        
        if isempty(GUImain.Wall_reflEq)
            statusbarObj.setText('Input formula of the reflection coefficient dependent on  frequency!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_inputEq);
            return
        end
        
        if isempty(strfind(GUImain.Wall_reflEq,'f'))
            statusbarObj.setText('Wrong format of the reflection coefficient. Use f (Frequency) as a variable');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.togglebutton_inputEq);
            return
        end
    end
    
end

if isempty(GUImain.FBL)
    statusbarObj.setText('Input a length of the Fourier boundary!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_FBL);
    return;
end

if  length(GUImain.FBL)~=2
    statusbarObj.setText('Specify length of FBL at the left and right bdy. [L;R].');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_FBL);
    return;
end

if  any(GUImain.FBL)<0
    statusbarObj.setText('Error: the length of Fourier boundary is negative!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_FBL);
    return;
end

if strcmp(GUImain.evol,'HS1')||strcmp(GUImain.evol,'SWE1')
else
    if isempty(GUImain.NonlinAdj)
        statusbarObj.setText('Input a length of the Nonlinear Adjustment!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Nonlin_Adj);
        return;
    end
    
    if  length(GUImain.NonlinAdj)~=1
        statusbarObj.setText('Input a value of the length of the Nonlinear Adjustment!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Nonlin_Adj);
        return;
    end
    
    if GUImain.NonlinAdj<0
        statusbarObj.setText('the length of the nonlinear adjustment is negative!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Nonlin_Adj);
        return;
    end
    
    
end


if GUImain.shippar.check==1
    
    shipparam=GUImain.shippar.data;
    user_shapedata=GUImain.shippar.user_shapedata;
    
    Ns=length(shipparam(:,1));
    
    for ii=1:Ns
        if any(isempty(cell2mat(shipparam(ii,3))))||...
                any(isempty(cell2mat(shipparam(ii,4))))||...
                any(isempty(cell2mat(shipparam(ii,5))))||...
                any(isempty(cell2mat(shipparam(ii,6))))||...
                any(isempty(cell2mat(shipparam(ii,7))))
%               ||any(isnan(cell2mat(shipparam(ii,3:7))))
            statusbarObj.setText('Specify ship parameter');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_ship_setupGui);
            return;
        elseif any(strcmp(shipparam(:,1),'User-defined')) && isempty(user_shapedata)
            statusbarObj.setText('Specify ship shape data...');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_ship_setupGui);
            return;
        end
    end
    
    if strcmp(shipparam(ii,2),'Pitch') ||strcmp(shipparam(ii,2),'Free')
        try
            if isempty(cell2mat(shipparam(ii,8))) || isnan(cell2mat(shipparam(ii,8)))
                statusbarObj.setText('Specify z center rotation');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.pushbutton_ship_setupGui);
                return;
            end
        catch
            if isempty(cell2mat(shipparam(ii,8))) || isnan(str2num(cell2mat(shipparam(ii,8))))
                statusbarObj.setText('Specify z center rotation');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.pushbutton_ship_setupGui);
                return;
            end
            
        end
    end
    if strcmp(shipparam(ii,9),'yes')
        try
            if strcmp(shipparam(ii,2),'Surge')
                if isempty(cell2mat(shipparam(ii,10))) || isnan(cell2mat(shipparam(ii,10)))
                    statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                    statusbarTxt = statusbarObj.getComponent(0);
                    statusbarTxt.setForeground(java.awt.Color.red);
                    uicontrol(handles.pushbutton_ship_setupGui);
                    return;
                end
            elseif strcmp(shipparam(ii,2),'Heave')
                if isempty(cell2mat(shipparam(ii,11))) || isnan(cell2mat(shipparam(ii,11)))
                    statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                    statusbarTxt = statusbarObj.getComponent(0);
                    statusbarTxt.setForeground(java.awt.Color.red);
                    uicontrol(handles.pushbutton_ship_setupGui);
                    return;
                end
            elseif strcmp(shipparam(ii,2),'Pitch')
                if isempty(cell2mat(shipparam(ii,12))) || isnan(cell2mat(shipparam(ii,12)))
                    statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                    statusbarTxt = statusbarObj.getComponent(0);
                    statusbarTxt.setForeground(java.awt.Color.red);
                    uicontrol(handles.pushbutton_ship_setupGui);
                    return;
                end
            else
                if isempty(cell2mat(shipparam(ii,10))) || isnan(cell2mat(shipparam(ii,10)))...
                        ||isempty(cell2mat(shipparam(ii,11))) || isnan(cell2mat(shipparam(ii,11)))...
                        ||isempty(cell2mat(shipparam(ii,12))) || isnan(cell2mat(shipparam(ii,12)))
                    statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                    statusbarTxt = statusbarObj.getComponent(0);
                    statusbarTxt.setForeground(java.awt.Color.red);
                    uicontrol(handles.pushbutton_ship_setupGui);
                    return;
                end
            end
        catch
           if strcmp(shipparam(ii,2),'Surge') 
            if isempty(str2num(cell2mat(shipparam(ii,10)))) || isnan(str2num(cell2mat(shipparam(ii,10))))
                statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.pushbutton_ship_setupGui);
                return;
            end
           elseif strcmp(shipparam(ii,2),'Heave') 
            if isempty(str2num(cell2mat(shipparam(ii,11)))) || isnan(str2num(cell2mat(shipparam(ii,11))))
                statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.pushbutton_ship_setupGui);
                return;
            end
           elseif strcmp(shipparam(ii,2),'Pitch') 
            if isempty(str2num(cell2mat(shipparam(ii,12)))) || isnan(str2num(cell2mat(shipparam(ii,12))))
                statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.pushbutton_ship_setupGui);
                return;
            end
           else
               if isempty(str2num(cell2mat(shipparam(ii,10)))) || isnan(str2num(cell2mat(shipparam(ii,10))))...
                       ||isempty(str2num(cell2mat(shipparam(ii,11)))) || isnan(str2num(cell2mat(shipparam(ii,11))))...
                       ||isempty(str2num(cell2mat(shipparam(ii,12)))) || isnan(str2num(cell2mat(shipparam(ii,12))))
                   statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                   statusbarTxt = statusbarObj.getComponent(0);
                   statusbarTxt.setForeground(java.awt.Color.red);
                   uicontrol(handles.pushbutton_ship_setupGui);
                   return;
               end
           end
        end
    end
    
    if isfield(GUImain.shippar,'FlagCloseGui')
        if GUImain.shippar.FlagCloseGui==1
            statusbarObj.setText('Please close the ship-setup panel');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.pushbutton_ship_setupGui);
            return;
        end
    end
end


if GUImain.InteriorProp==1
    if isempty(GUImain.IP_time)
        statusbarObj.setText('Input a time interval for the internal properties!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Interior_Time);
        return;
    end
    
    if length(GUImain.IP_time)~=2
        statusbarObj.setText('Wrong input format of the time interval for the internal properties!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Interior_Time);
        return;
    end
    
    if  GUImain.IP_time(1)>GUImain.IP_time(2)
        statusbarObj.setText('Error: Wrong input format of the time interval of the internal properties!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Interior_Time);
        return;
    end
    
    if  GUImain.IP_time(1)<GUImain.tinterv(1)
        statusbarObj.setText('Error: t_start of Internal Properties < t_start of the main simulation!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Interior_Time);
        return;
    end
    
    
    if  GUImain.IP_time(2)>GUImain.tinterv(2)
        statusbarObj.setText('Error: t_end of Internal Properties > t_end of the main simulation!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_Interior_Time);
        return;
    end
    
    if isempty(GUImain.IP_dt)
        statusbarObj.setText('Input a time step for the internal properties (#dt)!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_dtIP);
        return;
    end
    
    if length(GUImain.IP_dt)~=1
        statusbarObj.setText('Input  only a value of the time step for the internal properties (#dt)!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_dtIP);
        return;
    end
    
    if  GUImain.IP_dt<0
        statusbarObj.setText('Error: time step for the internal properties is negative!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_dtIP);
        return;
    end
    
    if  GUImain.IP_dt==0
        statusbarObj.setText('Error: time step for the internal properties is zero!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_dtIP);
        return;
    end
    if  fix(GUImain.IP_dt)~=GUImain.IP_dt
        statusbarObj.setText('Input an integer value for the time step of the internal properties!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_dtIP);
        return;
    end
    
    
end

if GUImain.partitionId==0
    if isempty(GUImain.partition)
        statusbarObj.setText('Input a number of partition!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_partition);
        return;
    end
    
    if length(GUImain.partitionId)~=1
        statusbarObj.setText('Input only a value of the number of partition!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_partition);
        return;
    end
    
    if fix(GUImain.partition)~=GUImain.partition
        statusbarObj.setText('Error:#partition is not an integer');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_partition);
        return;
    end
    
    if  GUImain.partition<=0
        statusbarObj.setText('Error: #partition should be larger than 0');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_partition);
        return;
    end
end


if exist(GUImain.proj.workdir,'dir')
    if FlagWarn==0;
        statusbarObj.setText('Warning: Project exists already, it will be overwritten!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
        uicontrol(handles.edit_savename);
        FlagWarn=1;
        return
    else
        FlagWarn=0;
        statusbarObj.setText('');
    end
else
    statusbarObj.setText('');
    FlagWarn=0;
end


GUImain.FlagGui=1;

%uiresume(gcbf);
save_state(GUImain); %save all GUI input for future use.
GUImain=adjust_variables(GUImain); %adjust some variables to be same as in MainUser.m

statusbarObj.setText('Pre-processing...');

JavProgressBar;

set(jProgressBar,'Maximum',100, 'Value',0);
jProgressBar.setStringPainted( true );


MainUser;

statusbarObj.setText('Pre-processing done.');

set(handles.pushbutton_run,'enable','on');
set(handles.togglebutton_logfile,'enable','on');
set(handles.togglebutton_logfile,'userdata',[]);


if get(handles.WaveType,'Value')~=1
    numerical_setting_advises(handles,statusbarObj,jProgressBar);
end


%run([GUImain.path,'\Program\','MainUser.m']);



function numerical_setting_advises(handles,statusbarObj,jProgressBar)
preproc=get(handles.pushbutton1,'userdata');
Proj=preproc.Proj;
par=preproc.par;
influx=preproc.influx;
model=preproc.model;
Id_disp=get(handles.popupmenu_dispersion,'value');


if par.ode_Npartition>1
    set(handles.checkbox_combineFiles,'enable','on');
end

depth=influx.depth;
nFreqP=model.nonlinear;
%if nFreqP==1, nFreqP=nFreqP+1; end

Hs=influx.Hs;
% if Hs/2>=depth
%     jProgressBar.setVisible(0);
%     set(handles.pushbutton_run,'enable','off')
%     statusbarObj.setText('Error: Hs./2 >= depth');
%     statusbarTxt = statusbarObj.getComponent(0);
%     statusbarTxt.setForeground(java.awt.Color.red);
%     return;
% end

xleft=par.x(1);
xright=par.x(end);

pnow =round(log2(par.Nx));
cutfrac_now=par.cutfrac;

nu_peak=influx.nu_p;

%Om          = str2func(model.dispersion);
%omAdd       = model.OmFun;
Om          = str2func('OmExact');
omAdd       ='';
k_cut=invOm((nFreqP)*nu_peak,depth,Om,omAdd); %


cutfrac=nFreqP*2;  %nyquist k=1/2 kmax  as precaution use k=1/3 max
k_max=cutfrac*k_cut;

dx=2*pi/k_max;
Nx=(xright-xleft)./dx;
padvise=round(log2(Nx))+1;
dxadvise=(xright-xleft)*100000/((2^padvise-1))/100000;
cutfrac_advise=cutfrac;

if isnan(padvise)...  %For BBM
        || k_cut<0    %For KdV
    statusbarObj.setText('The dispersion is not reliable for the input wave frequencies');
    disp('-------------------------------Warning-------------------------------');
    disp('The dispersion is not reliable for the input wave frequencies');
    disp('---------------------------------------------------------------------');
    disp('%');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    set(handles.pushbutton_run,'enable','off');
    return;
end


kp_h=influx.kh;
if Id_disp==3
    if kp_h>1
        statusbarObj.setText('The dispersion is not reliable for the case of k*h>1');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        disp('-------------------------------Warning-------------------------------');
        disp('The dispersion is not reliable for the case of k*h>1');
        disp('---------------------------------------------------------------------');
        disp('%');
        
        set(handles.pushbutton_run,'enable','off');
        return;
    else
        statusbarObj.setText('');
    end
elseif Id_disp==4
    if kp_h>1.5
        statusbarObj.setText('The dispersion is not reliable for the case of k*h>3/2');
        
        disp('-------------------------------Warning-------------------------------');
        disp('The dispersion is not reliable for the case of k*h>3/2');
        disp('---------------------------------------------------------------------');
        disp('%');
        
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        set(handles.pushbutton_run,'enable','off');
        return;
    else
        statusbarObj.setText('');
    end
else
    statusbarObj.setText('');
end

diary off;

Fbdy_now=par.bf0;
lambda_peak=influx.lambda_p;

Fbdy_advise=floor(2*lambda_peak);

NonlinAdj_now=get(handles.edit_Nonlin_Adj,'userdata');
ka=influx.ka;
relWaveLength=influx.lambda_per_H;
NonlinAdj_advise=LengthNonlinAdj_advise(relWaveLength,ka);

if isempty(NonlinAdj_now)
    NonlinAdj_now=0;
end

if nu_peak>0.001
    
    if (pnow ~=padvise) || (cutfrac_now ~=cutfrac_advise) ...
            || (Fbdy_now(1) < Fbdy_advise) || (Fbdy_now(2) < Fbdy_advise) ...
            ||(NonlinAdj_now < NonlinAdj_advise)
        contents = cellstr(get(handles.dynModel,'String'));
        dynmodel=contents{get(handles.dynModel,'Value')};
        
        if strcmp(dynmodel,'HS1')
            statusbarObj.setText(['advise: dx<=',num2str(dxadvise),...
                '; cutfrac k=',num2str(cutfrac), ...
                '; Fourier Bdy.>=',num2str(Fbdy_advise)]);
        else
            statusbarObj.setText(['advise: dx<=',num2str(dxadvise),...
                '; cutfrac k=',num2str(cutfrac), ...
                '; Fourier Bdy.>=',num2str(Fbdy_advise), ...
                '; Nonlin. Adj.>=',num2str(NonlinAdj_advise),'.lambda_peak']);
        end
        
        
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
        
        % if pnow<padvise || cutfrac_now <cutfrac_advise
        % set(handles.pushbutton_run,'enable','off')
        % end
    end
else
    Influxdir=cellstr(get(handles.popupmenu_influx_direction,'String'));
    infdir=Influxdir{get(handles.popupmenu_influx_direction,'Value')};
    
    if strcmp(infdir,'Bi')==0
        statusbarObj.setText('Use only bi-directional propagation influx method.')
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        set(handles.pushbutton_run,'enable','off');
        set(handles.togglebutton_logfile,'enable','off');
        uicontrol(handles.popupmenu_influx_direction);
        return;
    else
        statusbarObj.setText('');
    end
end

function NonlinAdj_advise=LengthNonlinAdj_advise(relWaveLength,ka)
if relWaveLength<1/80 || ka<0.0075
    NonlinAdj_advise=2;
elseif  relWaveLength<1/40 || ka<0.015
    NonlinAdj_advise=3;
elseif  relWaveLength<1/20 || ka<0.03
    NonlinAdj_advise=4;
elseif  relWaveLength<1/10 || ka<0.06
    NonlinAdj_advise=5;
elseif  relWaveLength<1/2 || ka<0.12
    NonlinAdj_advise=7;
else
    NonlinAdj_advise=7;
end
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
%curdir=uigetdir(matlabroot);
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


% --- Executes on button press in checkbox_interior.
function checkbox_interior_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_interior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_interior
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_Interior_Time,'enable','on')
    set(handles.edit_dtIP,'enable','on')
else
    set(handles.edit_dtIP,'enable','off')
    set(handles.edit_Interior_Time,'enable','off')
end

% --- Executes on button press in togglebutton_logfile.
function togglebutton_logfile_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_logfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_logfile

preproc=get(handles.pushbutton1,'userdata');
par=preproc.par;
influx=preproc.influx;
model=preproc.model;
bath =preproc.bath;
Proj =preproc.Proj;
input=preproc.input;
IVP  =preproc.IVP;
options=preproc.options;
shippar=preproc.shippar;

Tog=get(hObject,'Value');

log.string{1}='*******************************************************************************';
log.string{2}=['HaWaSSI-AB simulation diary, date: ' ...
    num2str(par.fc(3)),'.' num2str(par.fc(2)),'.' num2str(par.fc(1)),...
    ' time:  ' num2str(par.fc(4)),'h.' num2str(par.fc(5)),'min.'];
log.string{3}='*******************************************************************************';
log.string{4}=['Project name : ',Proj.savename];
log.string{5}=['User Note     : ',Proj.UseNote];

log.string{7}='MODEL DESCRIPTION :';
log.string{8}=['Dynamic Model        : ', model.dyn(1:3)];
log.string{9}=['Dispersion Model     : ', model.dispersion(3:end)];
IdN=9;
if strcmp(model.breaking.check,'Yes')
    log.string{IdN+1}=['Breaking                  : '];
    log.string{IdN+2}=['   Initiation                : U/C= ', num2str(model.breaking.KBC)];
    log.string{IdN+3}=['   Termination          : uF/uI= ', num2str(model.breaking.TC),', T*/Tp= ', num2str(model.breaking.Tchar)];
    IdN=IdN+4;
else
    log.string{IdN+1}=['Breaking                  : No'];
    IdN=IdN+2;
end

log.string{IdN+1}='INFLUX DESCRIPTION : ';
if input.type==1
    log.string{IdN+2}= ['Signal type : ',influx.name];
end
if input.type~=1
    log.string{IdN+2}=['Depth (h)                                 : ', num2str(influx.depth),'[m]'];
    if strcmp(influx.name,'Harmonic')||strcmp(influx.name,'User-defined')
        log.string{IdN+3}= ['Signal type                              : ',influx.name];
    elseif strcmp(influx.name,'Jonswap')
        log.string{IdN+3}= ['Signal type                              : ',influx.name,' with gamma= ',num2str(influx.JS_gamma)];
    end
    if strcmp(influx.name,'Harmonic')
        log.string{IdN+4}= ['Amplitude (A)                          : ',num2str(roundn(influx.Hs/2,-3)),'[m]'];
    else
        log.string{IdN+4}= ['Significant wave Height (Hs) : ',num2str(roundn(influx.Hs,-3)),'[m]'];
    end
    log.string{IdN+5}= ['Peak period (Tp)                     : ',num2str(roundn(influx.Tp,-3)),'[s]'];
    log.string{IdN+6}= ['         Derived info : '];
    log.string{IdN+7}= ['         Peak frequency (nu)                     : ',num2str(roundn(influx.nu_p,-3)),'[rad/s]'];
    log.string{IdN+8}= ['         Peak wave-number (kp)               : ',num2str(roundn(influx.k_p,-3))];
    log.string{IdN+9}= ['         Peak wave-length (lambda)          : ',num2str(roundn(influx.lambda_p,-3)),'[m]'];
    log.string{IdN+10}= ['         Peak phase speed                        : ',num2str(roundn(influx.Cp_p,-3)),'[m/s]'];
    log.string{IdN+11}= ['         Peak group speed                        : ',num2str(roundn(influx.Vg_p,-3)),'[m/s]'];
    log.string{IdN+12}= ['         Steepness (kp*(Hs./2))                 : ',num2str(roundn(influx.ka,-3))];
    log.string{IdN+13}= ['         Relative wave-length (lambda/h)  : ',num2str(influx.lambda_per_H),' (',influx.category,')'];
    log.string{IdN+14}= ['                                            (kp*h)         : ', num2str(roundn(influx.kh,-3))];
    IdN=IdN+15;
else
    IdN=IdN+3;
end

log.string{IdN+1}=['INITIAL WAVE CONDITIONS : '];
log.string{IdN+2}=['Initial condition     : ',IVP.typename];
IdN=IdN+3;
if  IVP.type==2||IVP.type==3
    log.string{IdN+1}=['Amplitude              : ',num2str(IVP.A)];
    log.string{IdN+2}= ['Standard deviation : ',num2str(IVP.lambda)];
    log.string{IdN+3}= ['Center position    : ',num2str(IVP.x0)];
    IdN=IdN+4;
end

if IVP.type==4
    log.string{IdN}=['File name             : ',IVP.filename];
    IdN=IdN+1;
end


if influx.wind.check==1
    log.string{IdN+1}=['WIND INPUT   : YES'];
    log.string{IdN+2}=['     Coef    : ',num2str(influx.wind.coef)];
    log.string{IdN+3}=['     tinterv : [',num2str(influx.wind.tinterv(1)),';',num2str(influx.wind.tinterv(2)),'] [s]'];
    log.string{IdN+4}=['     xinterv : [',num2str(influx.wind.xinterv(1)),';',num2str(influx.wind.xinterv(2)),'] [m]'];
    IdN=IdN+5;
else
    log.string{IdN+1}=['WIND INPUT : NO'];
    IdN=IdN+2;
end

log.string{IdN+1}='NUMERICAL SETTINGS : ';
log.string{IdN+2}=['Spatial interval               : ','(', num2str(par.Xleft),';',num2str(par.Xright),') [m]'];
log.string{IdN+3}=['Length of Fourier bdy  [L;R] : [',num2str(par.bf0(1)),';',num2str(par.bf0(2)),'] [m]'];
log.string{IdN+4}=['Number of Nodes         : ',num2str(par.Nx)];
log.string{IdN+5}=['    Grid size (dx)           : ',num2str(par.dx),' [m]'];
log.string{IdN+6}=['    Cutfrac k                  : ',num2str(par.cutfrac)];
log.string{IdN+7}=['Time interval                 : ','(', num2str(par.t_init),';',num2str(par.t_end),') [s]'];
log.string{IdN+8}=['    Time step (dt)           : ',num2str(roundn(par.dt,-3)),' [s]'];
IdN=IdN+9;
if strcmp(bath.type,'F')
    log.string{IdN+1}=['Bathymetry     : ', bath.name];
    log.string{IdN+2}=['     Depth         : ', num2str(bath.depth),'[m]'];
    IdN=IdN+3;
elseif strcmp(bath.type,'B')
    log.string{IdN+1}=['Bathymetry     : ', bath.name];
    log.string{IdN+2}=['     Depth         : max= ', num2str(bath.par(1)),'[m],', ...
        ' min= ', num2str(bath.par(2)),'[m]'];
    log.string{IdN+3}=['     Slope         : gradient= ', num2str(bath.par(3)), ...
        ', start position : ', num2str(bath.par(4)),'[m]'];
    IdN=IdN+4;
elseif strcmp(bath.type,'BR')
    log.string{IdN+1}=['Bathymetry     : ', bath.name];
    log.string{IdN+2}=['     Depth         : Max= ', num2str(bath.par(1)),'[m]']; ...
        log.string{IdN+3}=['     Slope         : gradient= ', num2str(bath.par(2)), ', shore position : ', num2str(bath.par(3)),'[m]'];
    IdN=IdN+4;
else
    log.string{IdN+1}=['Bathymetry     : ', bath.name];
    log.string{IdN+2}=['     File-name    : ', bath.filename];
    IdN=IdN+3;
end

if bath.friction.check==1
    strF=bath.friction.data;
    NsF=length(strF.cf);
    for i=1:NsF-1
        strCf{i}=[strF.cf{i},' , ',];
        strInterv{i}=['(',strF.interval{i},'),'];
    end
    strCf{NsF}=[strF.cf{NsF}];
    strInterv{NsF}=['(',strF.interval{NsF},')'];
    
    log.string{IdN+1}=['Bottom friction: '];
    log.string{IdN+2}=['      coef         : ', strCf{1:NsF}];
    log.string{IdN+3}=['      interval     : ', strInterv{1:NsF},' [m]'];
    IdN=IdN+4;
else
    log.string{IdN+1}=['Bottom friction: No'];
    IdN=IdN+2;
end

if strcmp(bath.wall.check,'Yes')
    if bath.wall.type==1
        log.string{IdN+1}=['Wall: '];
        log.string{IdN+2}=['    Type       : Uniform'];
        log.string{IdN+3}=['    Position   : ',num2str(bath.wall.position), '[m]'];
        log.string{IdN+4}=['    Coef       : ',num2str(bath.wall.Coef)];
    else
        log.string{IdN+1}=['Wall:'];
        log.string{IdN+2}=['    Type       : Frequency dependent'];
        log.string{IdN+3}=['    Position   : ',num2str(bath.wall.position), '[m]'];
        log.string{IdN+4}=['    Function  : R(f)=', bath.wall.file_R];
    end
    IdN=IdN+5;
else
    log.string{IdN+1}=['Wall: No'];
    IdN=IdN+2;
end

if ~strcmp(model.influx.direction,'None')
    if strcmp(model.influx.type,'Point')
        log.string{IdN+1}=['Influx generation method: ', model.influx.direction,'; ', model.influx.type];
    else
        log.string{IdN+1}=['Influx generation method: ', model.influx.direction,'; ', model.influx.type(3:end)];
    end
    if ~isempty(bath.influx_AdjZone)
        log.string{IdN+2}=['    Nonlinear adjustment  : ', num2str(bath.influx_AdjZone), '*lambda_peak'];
    else
        log.string{IdN+2}=['    Nonlinear adjustment  : 0*lambda_peak'];
    end
    if input.ramp.check==1
        log.string{IdN+3}=['    Ramp                          : ', num2str(input.ramp.length), '*Tp'];
    else
        log.string{IdN+3}=['    Ramp                          : Off'];
    end
    
    if input.filter.check==1
        log.string{IdN+4}=['    Filter                           : On; Freq. interval: (',num2str(input.filter.LFreq),';',num2str(input.filter.HFreq),') rad/s'];
    else
        log.string{IdN+4}=['    Filter                           : Off'];
    end
    log.string{IdN+5}=['    Influx position             : ', num2str(influx.position),'[m]'];
    IdN=IdN+6;
else
    log.string{IdN+1}=['Influx generation method: ', model.influx.direction];
    IdN=IdN+2;
end


if shippar.check==1
    log.string{IdN+1}=['SHIP SETUP'];
    log.string{IdN+2}=['Nship: ',num2str(shippar.Nship)];
    IdN=IdN+2;iDD=1;
    for ii=1:shippar.Nship
        if strcmpi(shippar.data(ii,1),'Barge')
            shipshape='Barge';
        elseif strcmpi(shippar.data(ii,1),'Half-circle')
            shipshape='Half-circle';
        else
            shipshape='User-defined';
        end
        
        if strcmp(shippar.data(ii,2),'Fixed')
            shipMot='Fixed';
            AdMass=0;
        elseif strcmp(shippar.data(ii,2),'Heave')
            shipMot='Heave';
            AdMass=shippar.rad.Ma.z(ii);
        elseif strcmp(shippar.data(ii,2),'Surge')
            shipMot='Surge';
            AdMass=shippar.rad.Ma.x(ii);
        elseif strcmp(shippar.data(ii,2),'Pitch')
            shipMot='Pitch';
            AdMass=shippar.rad.Ma.theta(ii);
        else
            shipMot='Free';
            AdMassx=shippar.rad.Ma.x(ii);
            AdMassz=shippar.rad.Ma.z(ii);
            AdMasstheta=shippar.rad.Ma.theta(ii);
        end
        
        log.string{IdN+iDD}=['Ship #',num2str(ii),': Shape=',shipshape, ', Motion=',shipMot,...
            ', Length=',num2str(cell2mat(shippar.data(ii,3))),'[m] , Draft=',num2str(cell2mat(shippar.data(ii,4))),'[m]'];
        log.string{IdN+iDD+1}=['      Initial center ship position(Xc,Zc)= (',num2str(cell2mat(shippar.data(ii,5))),',',num2str(cell2mat(shippar.data(ii,6))),')[m]'];
        if strcmpi(shipMot,'Free')
            log.string{IdN+iDD+2}=['      Mass=', num2str(shippar.form.Mass(ii)), 'Moment Inertia=',num2str(shippar.form.MIner(ii)),...
                ', Added Mass (heave)= ', num2str(AdMassz),'(surge)= ', num2str(AdMassx),'(pitch)= ', num2str(AdMasstheta)];
        elseif strcmpi(shipMot,'Pitch')
            log.string{IdN+iDD+2}=['      Moment Inertia=', num2str(shippar.form.MIner(ii)),...
                ', Added Moment Inertia = ', num2str(AdMass)];
        else
            log.string{IdN+iDD+2}=['      Mass=', num2str(shippar.form.Mass(ii)),...
                ', Added Mass = ', num2str(AdMass)];
        end
        iDD=iDD+3;
    end
    IdN=IdN+shippar.Nship*3+1;
else
    log.string{IdN+1}=['SHIP:  NO'];
    IdN=IdN+2;
end





if options.interior.check==1
    log.string{IdN+1}='INTERIOR CALCULATION PREPARATION: ';
    log.string{IdN+2}=['Time interval        : (',num2str(options.interior.time(1)),';',num2str(options.interior.time(2)),') [s]'];
    log.string{IdN+3}=['Time step            : ',num2str(par.dt.*options.interior.time(3)),' [s]'];
    IdN=IdN+4;
end

log.string{IdN+1}='ODE SETTING:';
log.string{IdN+2}=['ODE solver           : ', par.odesolver];
log.string{IdN+3}=['Relative error tol.  : ', num2str(par.odetol)];
log.string{IdN+4}=['No of Partition       : ', num2str(par.ode_Npartition)];
log.string{IdN+5}='------------------------------------------------------------------------------';
IdN=IdN+6;
if ~isempty(get(handles.togglebutton_logfile,'userdata'))
    comptime=get(handles.togglebutton_logfile,'userdata');
    log.string{IdN+1}=['After simulation:'];
    log.string{IdN+2}=['CPUtime for ODEs: ',num2str(comptime.cpu,'%10.2e\n'),'[sec]'];
    log.string{IdN+3}=['Relative computation time: ',num2str(comptime.crel,2),' percent'];
    log.string{IdN+4}=['Data saved as : ',Proj.savename,'_simul.mat'];
    log.string{IdN+5}='------------------------------------------------------------------------------';
end

if Tog==1
    log_file(log,handles.togglebutton_logfile);
else
    delete(log_file);
end

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenProject_Callback(hObject, eventdata, handles) %%% New
% hObject    handle to OpenProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata, handles)
HAWASSI_AB_startpage;


function input_handles_gui(state,handles)
% cd(directory);
% cd(pathnow);

%     if isfield(state,'edit_CurDir') %adjusted for GUI_state old version
%         set(handles.edit_CurDir,'String',state.edit_CurDir);
%     end
%
%     if exist(state.path,'file')
%         additemtolistbox(handles.popupmenu_CurDir,state.path);
%     end
global IDdemo

if handles.input.flagOpenProj==2
    additemtolistbox(handles.popupmenu_CurDir,state.proj.projdir);
end

if isfield(state,'IVP')
    if state.IVP.type==5  %Previous GUI setting
        state.IVP.type=4;
    end
    if state.IVP.type==1 || state.IVP.type==4
        if state.IVP.type==1
            set(handles.popupmenu_IVP,'value',1);
        else
            set(handles.popupmenu_IVP,'value',4);
        end
        set(handles.edit_IVP_A,'userdata',[]);
        set(handles.edit_IVP_A,'string','');
        set(handles.edit_IVP_A,'enable','off');
        set(handles.edit_IVP_lambda,'userdata',[]);
        set(handles.edit_IVP_lambda,'string','');
        set(handles.edit_IVP_lambda,'enable','off');
        
        set(handles.edit_IVP_x0,'userdata',[]);
        set(handles.edit_IVP_x0,'string','');
        set(handles.edit_IVP_x0,'enable','off');
        
        set(handles.pushbutton_load_file_IVP,'visible','off')
    end
    if state.IVP.type==2 || state.IVP.type==3
        if state.IVP.type==2
            set(handles.popupmenu_IVP,'value',2);
        elseif state.IVP.type==3
            set(handles.popupmenu_IVP,'value',3);
        end
        set(handles.edit_IVP_A,'userdata',state.IVP.A);
        set(handles.edit_IVP_A,'string',num2str(state.IVP.A));
        set(handles.edit_IVP_A,'enable','on');
        set(handles.edit_IVP_lambda,'userdata',state.IVP.lambda);
        set(handles.edit_IVP_lambda,'string',num2str(state.IVP.lambda));
        set(handles.edit_IVP_lambda,'enable','on');
        
        set(handles.edit_IVP_x0,'userdata',state.IVP.x0);
        set(handles.edit_IVP_x0,'string',num2str(state.IVP.x0));
        set(handles.edit_IVP_x0,'enable','on');
        set(handles.pushbutton_load_file_IVP,'visible','off')
    end
    
    if state.IVP.type==4
        set(handles.pushbutton_load_file_IVP,'userdata',state.IVP.file);
        set(handles.pushbutton_load_file_IVP,'visible','on')
    end
    
else
    set(handles.popupmenu_IVP,'value',1);
    set(handles.edit_IVP_A,'userdata',[]);
    set(handles.edit_IVP_A,'string','');
    set(handles.edit_IVP_A,'enable','off');
    set(handles.edit_IVP_lambda,'userdata',[]);
    set(handles.edit_IVP_lambda,'string','');
    set(handles.edit_IVP_lambda,'enable','off');
    
    set(handles.edit_IVP_x0,'userdata',[]);
    set(handles.edit_IVP_x0,'string','');
    set(handles.edit_IVP_x0,'enable','off');
    
    set(handles.pushbutton_load_file_IVP,'visible','off')
end

if ~isfield(state,'IVP')&&strcmp(state.WT,'Initial Value Problem')
    set(handles.popupmenu_IVP,'value',4);
    set(handles.edit_IVP_A,'userdata',[]);
    set(handles.edit_IVP_A,'string','');
    set(handles.edit_IVP_A,'enable','off');
    set(handles.edit_IVP_lambda,'userdata',[]);
    set(handles.edit_IVP_lambda,'string','');
    set(handles.edit_IVP_lambda,'enable','off');
    
    set(handles.edit_IVP_x0,'userdata',[]);
    set(handles.edit_IVP_x0,'string','');
    set(handles.edit_IVP_x0,'enable','off');
    
    
    set(handles.pushbutton_load_file_IVP,'visible','on')
    if isfield(state,'WT_I_func')
        set(handles.pushbutton_load_file_IVP,'userdata',state.WT_I_func);
    end
end

if strcmp(state.WT,'Harmonic')
    stateWT=2;
elseif strcmp(state.WT,'Jonswap')
    stateWT=3;
elseif strcmp(state.WT,'User-defined')
    stateWT=4;
elseif strcmp(state.WT,'Other')  %adjustment for previous version GUImain_state
    stateWT=4;
else
    stateWT=1;
end



set(handles.WaveType,'Value',stateWT)

if stateWT==1
    set(handles.edit_Hs,'Enable','off');
    set(handles.edit_Tp,'Enable','off');
    set(handles.edit_Jsg,'Enable','Off');
    set(handles.checkbox_Filter,'Enable','off');
    set(handles.edit_LFHF,'visible','off');
    set(handles.text71,'visible','off');
    set(handles.text61,'visible','off');
    set(handles.checkbox_ramp,'Enable','off')
    set(handles.edit_ramp,'visible','off');
    set(handles.text101,'visible','off');
    set(handles.text100,'visible','off');
    set(handles.checkbox_Filter,'value',0);
    
elseif stateWT==2
    if isfield(state,'Hs')
        set(handles.edit_Hs,'Userdata',state.Hs/2);
        set(handles.edit_Tp,'Userdata',state.Tp);
        set(handles.edit_Hs,'String',num2str(state.Hs/2));
        set(handles.edit_Tp,'String',num2str(state.Tp));
    else
        set(handles.edit_Hs,'Userdata',state.WT_H_param(1)/2);
        set(handles.edit_Tp,'Userdata',state.WT_H_param(2));
        set(handles.edit_Hs,'String',num2str(state.WT_H_param(1)/2));
        set(handles.edit_Tp,'String',num2str(state.WT_H_param(2)));
    end
    set(handles.edit_Hs,'Enable','On');
    set(handles.edit_Tp,'Enable','On');
    set(handles.edit_Jsg,'Enable','Off');
    set(handles.WT_Other_data,'Visible','Off');
    
    set(handles.text_Hs,'Visible','Off');
    set(handles.text_AInflux,'Visible','On');
    
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    if get(handles.checkbox_ramp,'value')==1
        set(handles.edit_ramp,'visible','on');
        set(handles.text101,'visible','on');
        set(handles.text100,'visible','on');
    end
    
elseif stateWT==3
    if isfield(state,'Hs')
        set(handles.edit_Hs,'Userdata',state.Hs);
        set(handles.edit_Tp,'Userdata',state.Tp);
        set(handles.edit_Jsg,'Userdata',state.Jsg);
        set(handles.edit_Hs,'String',num2str(state.Hs));
        set(handles.edit_Tp,'String',num2str(state.Tp));
        set(handles.edit_Jsg,'String',num2str(state.Jsg));
    else
        set(handles.edit_Hs,'Userdata',state.WT_H_param(1));
        set(handles.edit_Tp,'Userdata',state.WT_H_param(2));
        set(handles.edit_Jsg,'Userdata',state.WT_H_param(3));
        set(handles.edit_Hs,'String',num2str(state.WT_H_param(1)));
        set(handles.edit_Tp,'String',num2str(state.WT_H_param(2)));
        set(handles.edit_Jsg,'String',num2str(state.WT_H_param(3)));
    end
    set(handles.edit_Hs,'Enable','On');
    set(handles.edit_Tp,'Enable','On');
    set(handles.edit_Jsg,'Enable','On');
    set(handles.WT_Other_data,'Visible','Off');
    
    set(handles.text_Hs,'Visible','On');
    set(handles.text_AInflux,'Visible','Off');
    
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1;
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    if get(handles.checkbox_ramp,'value')==1;
        set(handles.edit_ramp,'visible','on');
        set(handles.text101,'visible','on');
        set(handles.text100,'visible','on');
    end
    
elseif stateWT==4
    set(handles.edit_Hs,'Enable','Off');
    set(handles.edit_Tp,'Enable','Off');
    set(handles.edit_Jsg,'Enable','Off');
    set(handles.WT_Other_data,'Userdata',state.WT_O_data);
    set(handles.WT_Other_data,'Visible','On');
    
    set(handles.checkbox_Filter,'Enable','on');
    if get(handles.checkbox_Filter,'value')==1
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
    end
    set(handles.checkbox_ramp,'Enable','on')
    if get(handles.checkbox_ramp,'value')==1;
        set(handles.edit_ramp,'visible','on');
        set(handles.text101,'visible','on');
        set(handles.text100,'visible','on');
    end
    
end

if stateWT~=1
    set(handles.checkbox_Filter,'Value',state.Filter);
    if state.Filter==1
        set(handles.edit_LFHF,'visible','on');
        set(handles.text71,'visible','on');
        set(handles.text61,'visible','on');
        set(handles.edit_LFHF,'Userdata',state.Filter_LFHF);
        set(handles.edit_LFHF,'Enable','On');
        set(handles.edit_LFHF,'String',[num2str(state.Filter_LFHF(1)),';' ...
            num2str(state.Filter_LFHF(2))]);
    else
        set(handles.edit_LFHF,'Enable','Off');
        set(handles.edit_LFHF,'visible','off');
        set(handles.text71,'visible','off');
        set(handles.text61,'visible','off');
    end
    
    if isfield(state,'ramp')
        set(handles.checkbox_ramp,'value',state.ramp.check)
        if state.ramp.check==1
            set(handles.edit_ramp,'userdata',state.ramp.val) ;
            set(handles.edit_ramp,'string',num2str(state.ramp.val)) ;
            set(handles.edit_ramp,'visible','on') ;
            set(handles.text100,'visible','on')
            set(handles.text101,'visible','on')
        else
            set(handles.text100,'visible','off')
            set(handles.text101,'visible','off')
            set(handles.edit_ramp,'visible','off') ;
            set(handles.checkbox_ramp,'value',0)
        end
    else
        set(handles.checkbox_ramp,'value',1)
        set(handles.edit_ramp,'userdata',4) ;
        set(handles.edit_ramp,'string','4') ;
        set(handles.edit_ramp,'visible','on') ;
        set(handles.text100,'visible','on')
        set(handles.text101,'visible','on')
    end
end
if strcmp(state.evol,'HS1')
    statedynM=1;
elseif strcmp(state.evol,'HS2')
    statedynM=2;
elseif strcmp(state.evol,'HS3')
    statedynM=3;
elseif strcmp(state.evol,'HS4')
    statedynM=4;
elseif strcmp(state.evol,'SWE1')||strcmp(state.evol,'LSWE')
    statedynM=5;
elseif strcmp(state.evol,'SWE2')||strcmp(state.evol,'SWE')
    statedynM=6;
end

if IDdemo==1
    if statedynM>2, statedynM=2;end
end
set(handles.dynModel,'Value',statedynM)

if IDdemo==1
    if state.Break==1,state.Break=0;end
end

set(handles.checkbox_Breaking,'Value',state.Break);
set(handles.checkbox_def_paramB,'Value',state.BreakDef);

if state.Break==1
    set(handles.Breaking_param,'Userdata',state.Break_param);
    set(handles.Breaking_param,'Visible','On');
    set(handles.Breaking_param_text,'Visible','On');
    set(handles.text72,'Visible','On');
    set(handles.edit_paramB_UF_T,'Visible','On');
    
    set(handles.Breaking_param,'String',num2str(state.Break_param));
    if ~isempty(state.Break_paramDef)
        set(handles.edit_paramB_UF_T,'String',[num2str(state.Break_paramDef(1)),';',...
            num2str(state.Break_paramDef(2))]);
        set(handles.edit_paramB_UF_T,'userdata',state.Break_paramDef);
    else
        set(handles.edit_paramB_UF_T,'String','0.2;0.5');
        set(handles.edit_paramB_UF_T,'userdata',[0.2;0.5]);
    end
    set(handles.checkbox_def_paramB,'Visible','On');
    if state.BreakDef==0
        set(handles.edit_paramB_UF_T,'Enable','On');
    else
        set(handles.edit_paramB_UF_T,'Enable','Off');
    end
else
    set(handles.Breaking_param,'Visible','Off');
    set(handles.Breaking_param_text,'Visible','Off');
    set(handles.checkbox_def_paramB,'Visible','Off');
    set(handles.edit_paramB_UF_T,'Visible','Off');
    set(handles.text72,'Visible','Off');
    set(handles.edit_paramB_UF_T,'String','0.2;0.5');
    set(handles.edit_paramB_UF_T,'userdata',[0.2;0.5]);
end

set(handles.text11,'Visible','On')
if strcmp(state.disp,'Exact'),dispers=1;
elseif strcmp(state.disp,'SWE'),dispers=2;
elseif strcmp(state.disp,'KdV'),dispers=3;
elseif strcmp(state.disp,'BBM'),dispers=4;
else
    dispers=5;
    set(handles.togglebutton_InputDispersion,'Userdata',state.disp_O_file);
    set(handles.togglebutton_InputDispersion,'Visible','On');
    set(handles.togglebutton_InputGroupVel,'Userdata',state.groupVel_O_file);
    set(handles.togglebutton_InputGroupVel,'Visible','On');
end

if IDdemo==1
    if dispers>1, dispers=1;end
end
set(handles.popupmenu_dispersion,'Value',dispers)

if strcmp(state.inftype,'Area'),inftype=1;
elseif strcmp(state.inftype,'Area-Short'),inftype=2;
elseif strcmp(state.inftype,'Point'),inftype=3;
end
set(handles.popupmenu_influx,'Value',inftype)

if strcmp(state.infdir,'Uni+'),infdir=1;
elseif strcmp(state.infdir,'Uni-'),infdir=2;
elseif strcmp(state.infdir,'Bi'),infdir=3;
end

set(handles.popupmenu_influx_direction,'Value',infdir);

if isfield(state,'tstep')
    set(handles.edit_time,'String',[num2str(state.tinterv(1)),';',...
        num2str(state.tinterv(2))]);
    set(handles.edit_time,'Userdata',[state.tinterv(1) state.tinterv(2)]);
    set(handles.edit_timestep,'String',num2str(state.tstep));
    set(handles.edit_timestep,'Userdata',state.tstep);
else
    set(handles.edit_time,'String',[num2str(state.tinterv(1)),';',...
        num2str(state.tinterv(3))]);
    set(handles.edit_time,'Userdata',[state.tinterv(1) state.tinterv(3)]);
    set(handles.edit_timestep,'String',num2str(state.tinterv(2)));
    set(handles.edit_timestep,'Userdata',state.tinterv(2));
end
if isfield(state,'Assim')
    if state.Assim.check==1
        set(handles.checkbox17,'value',1);
        set(handles.togglebutton_bdy_assim,'enable','on',...
            'userdata',state.Assim.data);
    else
        set(handles.checkbox17,'value',0);
        set(handles.togglebutton_bdy_assim,'enable','off');
    end
else
    set(handles.checkbox17,'value',0);
    set(handles.togglebutton_bdy_assim,'enable','off');
end

if isfield(state,'wind')
    set(handles.checkbox_wind,'value',state.wind.check);
    set(handles.edit_coef_wind,'userdata',state.wind.coef);
    set(handles.edit_coef_wind,'string',num2str(state.wind.coef));
    if state.wind.check==1
        if isfield(state.wind,'tinterv')
            set(handles.edit_wind_tinterv,'String',[num2str(state.wind.tinterv(1)),';',...
                num2str(state.wind.tinterv(2))],'userdata',state.wind.tinterv);
        else
            set(handles.edit_wind_tinterv,'String',[num2str(state.tinterv(1)),';',...
                num2str(state.tinterv(3))],'userdata',[state.tinterv(1);state.tinterv(3)]);
        end
        
        if isfield(state.wind,'xinterv')
            set(handles.edit_wind_xinterv,'String',[num2str(state.wind.xinterv(1)),';',...
                num2str(state.wind.xinterv(2))],'userdata',state.wind.xinterv);
        else
            set(handles.edit_wind_xinterv,'String',[num2str(state.xinterv(1)),';',...
                num2str(state.xinterv(3))],'userdata',[state.xinterv(1);state.tinterv(3)]);
        end
        
        set(handles.edit_coef_wind,'enable','on');
        set(handles.edit_wind_tinterv,'enable','on');
        set(handles.edit_wind_xinterv,'enable','on');
    else
        set(handles.edit_coef_wind,'enabl','off');
        set(handles.edit_wind_tinterv,'enable','off');
        set(handles.edit_wind_xinterv,'enable','off');
    end
end


if isfield(state,'Xinflux') %adjusted for GUI_state old version
    set(handles.edit_xinterv,'String',[num2str(state.xinterv(1)),';',...
        num2str(state.xinterv(2))]);
    set(handles.edit_xinterv,'Userdata',state.xinterv);
    set(handles.edit_Xinflux,'String',num2str(state.Xinflux));
    set(handles.edit_Xinflux,'Userdata',state.Xinflux);
else
    set(handles.edit_xinterv,'String',[num2str(state.xinterv(1)),';',...
        num2str(state.xinterv(3))]);
    tempstate.xinterv(1)=state.xinterv(1);
    tempstate.xinterv(2)=state.xinterv(3);
    set(handles.edit_xinterv,'Userdata',tempstate.xinterv);
    set(handles.edit_Xinflux,'String',num2str(state.xinterv(2)));
    set(handles.edit_Xinflux,'Userdata',state.xinterv(2));
end


if isfield(state,'cutfracK')
    if ~isfield(state,'dx')
        dx=(state.xinterv(2)-state.xinterv(1))./((2^state.p)-1);
    else
        dx=state.dx;
    end
    set(handles.edit_p,'String',num2str(dx));
    set(handles.edit_p,'Userdata',dx);
    set(handles.edit_cutfracK,'String',num2str(state.cutfracK));
    set(handles.edit_cutfracK,'Userdata',state.cutfracK);
else
    if ~isfield(state,'dx')
        dx=(state.xinterv(2)-state.xinterv(1))./((2^state.p(1))-1);
    else
        dx=state.dx;
    end
    set(handles.edit_p,'String',num2str(dx));
    set(handles.edit_p,'Userdata',dx);
    set(handles.edit_cutfracK,'String',num2str(state.p(2)));
    set(handles.edit_cutfracK,'Userdata',state.p(2));
end




if isfield(state,'NonlinAdj')
    if length(state.FBL)==1
        set(handles.edit_FBL,'String',[num2str(state.FBL),';',num2str(state.FBL)]);
        set(handles.edit_FBL,'Userdata',[state.FBL;state.FBL]);
    else
        set(handles.edit_FBL,'String',[num2str(state.FBL(1)),';',num2str(state.FBL(2))]);
        set(handles.edit_FBL,'Userdata',state.FBL);
    end
    
    set(handles.edit_Nonlin_Adj,'String',num2str(state.NonlinAdj));
    set(handles.edit_Nonlin_Adj,'Userdata',state.NonlinAdj);
    if statedynM>1
        set(handles.edit_Nonlin_Adj,'enable','on');
    end
    
else
    set(handles.edit_FBL,'String',[num2str(state.FBL(1)),';',num2str(state.FBL(1))]);
    set(handles.edit_FBL,'Userdata',[state.FBL(1);state.FBL(1)]);
    
    set(handles.edit_Nonlin_Adj,'String',num2str(state.FBL(2)));
    set(handles.edit_Nonlin_Adj,'Userdata',state.FBL(2));
    if statedynM>1
        set(handles.edit_Nonlin_Adj,'enable','on');
    end
end

if statedynM==1 ||statedynM==5
    set(handles.edit_Nonlin_Adj,'enable','off');
    set(handles.edit_Nonlin_Adj,'Userdata',[]);
    set(handles.edit_Nonlin_Adj,'string','');
else
    set(handles.edit_Nonlin_Adj,'enable','on');
end

if strcmp(state.bathy,'Flat'),bathy=1;
    set(handles.edit_depth_F,'enable','on');
    if isfield(state,'depth')
        set(handles.edit_depth_F,'Userdata',state.depth);
        set(handles.edit_depth_F,'String',num2str(state.depth));
    else
        set(handles.edit_depth_F,'Userdata',state.bathyF);
        set(handles.edit_depth_F,'String',num2str(state.bathyF));
    end
    set(handles.edit_slope,'Enable','Off');
    set(handles.pushbutton_load_file_bathy,'Enable','Off');
    set(handles.popupmenu_depth_reff,'enable','off');
    set(handles.text112,'enable','off');
elseif strcmp(state.bathy,'Slope'),bathy=2;
    set(handles.edit_depth_F,'enable','on');
    if isfield(state,'depth')
        set(handles.edit_depth_F,'Userdata',state.depth);
        set(handles.edit_depth_F,'String',[num2str(state.depth(1)),';',num2str(state.depth(2))]);
        set(handles.edit_slope,'Userdata',state.slope);
        set(handles.edit_slope,'String',[num2str(state.slope(1)),';',num2str(state.slope(2))]);
    else
        set(handles.edit_depth_F,'Userdata',[state.bathyS(1) state.bathyS(2)]);
        set(handles.edit_depth_F,'String',[num2str(state.bathyS(1)),';',num2str(state.bathyS(2))]);
        set(handles.edit_slope,'Userdata',[state.bathyS(3) state.bathyS(4)]);
        set(handles.edit_slope,'String',[num2str(state.bathyS(3)),';',num2str(state.bathyS(4))]);
    end
    set(handles.edit_slope,'Enable','On');
    set(handles.pushbutton_load_file_bathy,'Visible','Off');
    set(handles.popupmenu_depth_reff,'enable','on');
    set(handles.text112,'enable','on');
    if isfield(state,'bathyInterp')
        set(handles.popupmenu_depth_reff,'value',state.bathyInterp);
    else
        set(handles.popupmenu_depth_reff,'value',1);
    end
    
elseif strcmp(state.bathy,'Run-up (Slope)') ||strcmp(state.bathy,'Slope (Run-Up)')
    bathy=3;
    set(handles.edit_depth_F,'enable','on');
    if isfield(state,'depth')
        set(handles.edit_depth_F,'Userdata',state.depth,'tooltipstring','Specify a maximum depth and a minimum depth for shoreline (ie. 2% Hs). Format (max D; min D)');
        if length(state.depth)==2
        set(handles.edit_depth_F,'String',[num2str(state.depth(1)),';',num2str(state.depth(2))]);
        else
        set(handles.edit_depth_F,'String',[num2str(state.depth(1))]);    
        end
        
        set(handles.edit_slope,'Userdata',state.slope);
        set(handles.edit_slope,'String',[num2str(state.slope(1)),';',num2str(state.slope(2))]);
    else
        set(handles.edit_depth_F,'Userdata',state.bathyRS(1));
        set(handles.edit_depth_F,'String',num2str(state.bathyRS(1)));
        set(handles.edit_slope,'Userdata',[state.bathyRS(2) state.bathyRS(3)]);
        set(handles.edit_slope,'String',[num2str(state.bathyRS(2)),';',num2str(state.bathyRS(3))]);
    end
    set(handles.edit_slope,'Enable','On');
    set(handles.pushbutton_load_file_bathy,'Visible','Off');
    
    set(handles.popupmenu_depth_reff,'enable','on');
    set(handles.text112,'enable','on');
    %         if isfield(state,'bathyInterp')
    %             set(handles.popupmenu_depth_reff,'value',state.bathyInterp);
    %         else
    %             set(handles.popupmenu_depth_reff,'value',1);
    %         end
    set(handles.popupmenu_depth_reff,'value',2);
    
elseif strcmp(state.bathy,'User-defined')||strcmp(state.bathy,'Other')  %adjustment for previous version GUImain_state
    bathy=4;
    set(handles.edit_depth_F,'enable','off');
    set(handles.edit_slope,'enable','off');
    set(handles.pushbutton_load_file_bathy,'Userdata',state.bathyO);
    set(handles.pushbutton_load_file_bathy,'Visible','On');
    set(handles.pushbutton_load_file_bathy,'Enable','On');
    
    set(handles.popupmenu_depth_reff,'enable','on');
    
    set(handles.text112,'enable','on');
    if isfield(state,'bathyInterp')
        set(handles.popupmenu_depth_reff,'value',state.bathyInterp);
    else
        set(handles.popupmenu_depth_reff,'value',1);
    end
    if any(state.bathyO.data(:,2)>=0)
        set(handles.popupmenu_depth_reff,'value',2);
        set(handles.edit_depth_F,'enable','on','userdata',state.depth,'String',num2str(state.depth),...
            'tooltipstring','a minimum depth for shoreline (ie. 2% Hs)');
     end
    
elseif strcmp(state.bathy,'Run-up (Other)') ||strcmp(state.bathy,'User-defined (Run-Up)')
    bathy=4;
    set(handles.pushbutton_load_file_bathy,'Userdata',state.bathyO);
    set(handles.edit_slope,'Enable','Off');
    %set(handles.edit_depth_F,'enable','off');
    set(handles.edit_depth_F,'enable','on','userdata',state.depth,'String',num2str(state.depth),...
        'tooltipstring','a minimum depth for shoreline (ie. 2% Hs)');
    set(handles.pushbutton_load_file_bathy,'Visible','On');
    set(handles.pushbutton_load_file_bathy,'Enable','On')
    
    set(handles.popupmenu_depth_reff,'enable','on');
    set(handles.text112,'enable','on');
    %         if isfield(state,'bathyInterp')
    %             set(handles.popupmenu_depth_reff,'value',state.bathyInterp);
    %         else
    %             set(handles.popupmenu_depth_reff,'value',1);
    %         end
    set(handles.popupmenu_depth_reff,'value',2);
end

if bathy>1
    if get(handles.popupmenu_depth_reff,'value')==2
        if isfield(state,'DmidRef')
            set(handles.edit_Dmidref,'enable','on','userdata',state.DmidRef,'string',num2str(state.DmidRef));
            set(handles.text118,'enable','on');
            set(handles.text119,'enable','on');
        else
            if bathy==2 || (bathy==4&& ~any(state.bathyO.data(:,2)>=0) )
                if bathy==2
                    Dmidref=(state.depth(1)+state.depth(2))/2;
                else
                    Dmidref=(min(-state.bathyO.data(:,2))+max(-state.bathyO.data(:,2)))/2;
                end
            else
                if bathy==3
                    Dmidref=state.depth/6;
                else
                    Dmidref=(max(-state.bathyO.data(:,2)))/6;
                end
            end
            set(handles.edit_Dmidref,'enable','on','userdata',Dmidref,'string',num2str(Dmidref));
            set(handles.text118,'enable','on');
            set(handles.text119,'enable','on');
        end
    else
        set(handles.edit_Dmidref,'enable','off','userdata',[],'string','');
        set(handles.text118,'enable','off');
        set(handles.text119,'enable','off');
    end
    
end

set(handles.edit_slope,'backgroundcolor',[0.871 0.922 0.98]);
set(handles.popupmenu_bathy,'Value',bathy);


set(handles.checkbox_bottom_friction,'Value',0);
set(handles.togglebutton_friction,'Visible','off');
set(handles.togglebutton_friction,'Userdata',[]);

if isfield(state,'bottomfriction_check')
    if state.bottomfriction_check==1
        set(handles.checkbox_bottom_friction,'Value',1);
        set(handles.togglebutton_friction,'Visible','on');
        if isfield(state,'frictiondata')
            set(handles.togglebutton_friction,'Userdata',state.frictiondata);
        else %from previous setting
            strF.interval{1}=[num2str(state.bottomfriction_interv(1)),';',...
                num2str(state.bottomfriction_interv(2))];
            strF.cf{1}      =num2str(state.bottomfriction_cf);
            set(handles.togglebutton_friction,'Userdata',strF);
        end
    end
end



set(handles.checkbox_wall,'Value',state.Wall);

if state.Wall==1
    if isfield(state,'Wall_position')
        if isfield(state,'Wall_percent') %previous setting
            valCoef=state.Wall_percent/100;
        else
            valCoef=state.Wall_Coef;
        end
        
        set(handles.edit_Xwall,'Userdata',state.Wall_position);
        set(handles.edit_Xwall,'String',num2str(state.Wall_position));
        set(handles.edit_Refl_percent,'Userdata',valCoef);
        set(handles.edit_Refl_percent,'String',num2str(valCoef));
    else
        set(handles.edit_Xwall,'Userdata',state.Wall_param(1));
        set(handles.edit_Xwall,'String',num2str(state.Wall_param(1)));
        set(handles.edit_Refl_percent,'Userdata',state.Wall_param(2));
        set(handles.edit_Refl_percent,'String',num2str(state.Wall_param(2)));
    end
    
    if isfield(state,'Wall_type')
        set(handles.popupmenu_WallType,'Value',state.Wall_type);
        if isfield(state,'Wall_reflEq')
            set(handles.togglebutton_inputEq,'Userdata',state.Wall_reflEq);
        else
            set(handles.togglebutton_inputEq,'Userdata','1-min(f.^2,0.5)');
        end
    else
        set(handles.popupmenu_WallType,'Value',1);
    end
    
    IdWallType=get(handles.popupmenu_WallType,'Value');
    
    
    if IDdemo==1
        if IdWallType>1
            IdWallType=1;
            set(handles.popupmenu_WallType,'Value',1);
        end
    end
    
    if IdWallType==1
        set(handles.text84,'Visible','on');
        set(handles.edit_Refl_percent,'Visible','on');
        set(handles.togglebutton_inputEq,'Visible','off');
        set(handles.popupmenu_WallType,'Visible','on');
        set(handles.text86,'Visible','on');
    else
        set(handles.edit_Refl_percent,'Visible','off');
        set(handles.text84,'Visible','off');
        set(handles.togglebutton_inputEq,'Visible','on');
        
        set(handles.popupmenu_WallType,'Visible','on');
        set(handles.text86,'Visible','on');
    end
    set(handles.edit_Xwall ,'Visible','on');
    set(handles.text68,'Visible','on');
    set(handles.text83,'Visible','on');
    set(handles.popupmenu_wall_method,'Visible','on');
    set(handles.text_wall_method,'Visible','on');
    
    if isfield(state,'Wall_method')
     set(handles.popupmenu_wall_method,'value',state.Wall_method)   
    end
else
    set(handles.popupmenu_wall_method,'Visible','off');
    set(handles.text_wall_method,'Visible','off');
    
    set(handles.edit_Xwall ,'Visible','off');
    set(handles.popupmenu_WallType,'Visible','off');
    set(handles.edit_Refl_percent,'Visible','off');
    set(handles.togglebutton_inputEq,'Visible','off');
    set(handles.text68,'Visible','off');
    set(handles.text83,'Visible','off');
    set(handles.text84,'Visible','off');
    set(handles.text86,'Visible','off');
end

    
if handles.input.ShipFlag==0
    set(handles.checkbox_ship,'value',0,'enable','off');
    set(handles.pushbutton_ship_setupGui,'enable','off')
else
    
    if isfield(state,'shippar')
        if state.shippar.check==1
            set(handles.checkbox_ship,'value',1);
            set(handles.pushbutton_ship_setupGui,'userdata',state.shippar,'enable','on');
        else
            set(handles.checkbox_ship,'value',0);
            datIniship.data=[];
            datIniship.user_shapedata=[];
            datIniship.FlagCloseGui=0;
            set(handles.pushbutton_ship_setupGui,'userdata',datIniship,'enable','off');
        end
    else
        set(handles.checkbox_ship,'value',0);
        datIniship.data=[];
        datIniship.user_shapedata=[];
        datIniship.FlagCloseGui=0;
        set(handles.pushbutton_ship_setupGui,'userdata',datIniship,'enable','off');
    end
end
if isfield(state,'InteriorProp')
    set(handles.checkbox_interior,'value',state.InteriorProp);
    if state.InteriorProp==1
        if isfield(state,'IP_dt')
            set(handles.edit_Interior_Time,'Userdata',state.IP_time);
            set(handles.edit_Interior_Time,'String',[num2str(state.IP_time(1)),';',...
                num2str(state.IP_time(2))]);
            set(handles.edit_dtIP,'Userdata',state.IP_dt);
            set(handles.edit_dtIP,'string',num2str(state.IP_dt));
            set(handles.edit_Interior_Time,'Enable','On');
            set(handles.edit_dtIP,'Enable','On');
        else
            set(handles.edit_Interior_Time,'Userdata',[state.IP_time(1) state.IP_time(2)]);
            set(handles.edit_Interior_Time,'String',[num2str(state.IP_time(1)),';',...
                num2str(state.IP_time(2))]);
            set(handles.edit_dtIP,'Userdata',state.IP_time(3));
            set(handles.edit_dtIP,'string',num2str(state.IP_time(3)));
            set(handles.edit_Interior_Time,'Enable','On');
            set(handles.edit_dtIP,'Enable','On');
        end
        
    end
else
    set(handles.checkbox_interior,'value',0);
    set(handles.edit_Interior_Time,'Userdata',[]);
    set(handles.edit_Interior_Time,'String','');
    set(handles.edit_dtIP,'Userdata',[]);
    set(handles.edit_dtIP,'string','');
    set(handles.edit_Interior_Time,'Enable','Off');
    set(handles.edit_dtIP,'Enable','Off');
end

if isfield(state,'partitionId')
    if state.partitionId==0
        set(handles.checkbox_partition,'value',0);
        set(handles.edit_partition,'userdata',state.partition);
        set(handles.edit_partition,'string',num2str(state.partition));
        set(handles.edit_partition,'Enable','on');
        set(handles.checkbox_combineFiles,'enable','on');
    else
        set(handles.checkbox_partition,'value',1);
        set(handles.edit_partition,'userdata',[]);
        set(handles.edit_partition,'string','');
        set(handles.edit_partition,'Enable','off');
        set(handles.checkbox_combineFiles,'enable','off');
    end
end

if isfield(state,'partitionCombines')
    set(handles.checkbox_combineFiles,'value',state.partitionCombines);
else
    set(handles.checkbox_combineFiles,'value',1);
end

set(handles.edit_savename,'String',state.proj.name);
set(handles.edit_note,'String',state.proj.note);

Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
projdir=Curdir{get(handles.popupmenu_CurDir,'Value')};

global FlagWarn
%    [statusbarObj]=JavaFrame_handling();

if exist([projdir,'/',state.proj.name],'dir')
    %         statusbarObj.setText('Warning: Project exists already, it will be overwritten!');
    %         statusbarTxt = statusbarObj.getComponent(0);
    %         statusbarTxt.setForeground(java.awt.Color.blue);
    uicontrol(handles.edit_savename);
    FlagWarn=1;
end

%%%%%%%%%%%%%%%%%%%%%Demo version%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IDdemo==1
    %         statusbarObj.setText('Some features are not available in the demo version');
    %         statusbarTxt = statusbarObj.getComponent(0);
    %         statusbarTxt.setForeground(java.awt.Color.blue);
end



% --------------------------------------------------------------------
function SaveProject_Callback(hObject, eventdata, handles)
% hObject    handle to SaveProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};
saveName=get(handles.edit_savename,'String');
ProjectDir = uiputfile([pathNow,'\abproj_',saveName,'.mat'],'Save project file');
if ProjectDir~=0
    save_handles(handles,ProjectDir);
end
% --------------------------------------------------------------------
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reset_handles(handles,1)


% --------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
quit;

% --------------------------------------------------------------------
function Modules_Callback(hObject, eventdata, handles)
% hObject    handle to Modules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Doc_Callback(hObject, eventdata, handles)
% hObject    handle to Doc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if ~isdeployed
 winopen('\Toolbox\lib1d\Misc\Manual_AB_v1_1.pdf');
 else
 winopen('Manual_AB_v1_1.pdf');    
 end

% --------------------------------------------------------------------
function PostProc_Callback(hObject, eventdata, handles)
% hObject    handle to PostProc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input.Proj.path=handles.input.pathnow;
input.Proj.projdir=handles.input.projectdirectory;
PostProc_GUI(input);

% --------------------------------------------------------------------
function influxing_Callback(hObject, eventdata, handles)
% hObject    handle to influxing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Interior_Callback(hObject, eventdata, handles)
% hObject    handle to Interior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input.Proj.path=handles.input.pathnow;
input.Proj.projdir=handles.input.projectdirectory;
InteriorProp_GUI(input);

% --------------------------------------------------------------------
function Calc_Callback(hObject, eventdata, handles)
% hObject    handle to Calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Calculator();

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



function edit_Interior_Time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Interior_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Interior_Time as text
%        str2double(get(hObject,'String')) returns contents of edit_Interior_Time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Interior_Time,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Interior_Time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Interior_Time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Interior_Time (see GCBO)
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

set(handles.WaveType,'Value',1)
set(handles.edit_Hs,'Userdata',[]);
set(handles.edit_Hs,'String','');
set(handles.edit_Tp,'Userdata',[]);
set(handles.edit_Tp,'String','');
set(handles.edit_Jsg,'Userdata',[]);
set(handles.edit_Jsg,'String','');
set(handles.WT_Other_data,'Userdata',[]);

set(handles.pushbutton_load_file_IVP,'Userdata',[]);
set(handles.checkbox_Filter,'Value',0);
set(handles.edit_LFHF,'Userdata',[]);
set(handles.edit_LFHF,'String','');
set(handles.edit_Xinflux,'Userdata',[]);
set(handles.edit_Xinflux,'string','');
set(handles.dynModel,'Value',1);
set(handles.checkbox_Breaking,'Value',0);
set(handles.checkbox_def_paramB,'Value',1);
set(handles.Breaking_param,'Userdata',[]);
set(handles.Breaking_param,'String','');
set(handles.edit_paramB_UF_T,'string','0.2;0.5')
set(handles.edit_paramB_UF_T,'userdata',[0.2;0.5]);
set(handles.togglebutton_InputDispersion,'Userdata',[]);
set(handles.togglebutton_InputGroupVel,'Userdata',[]);
set(handles.popupmenu_dispersion,'Value',1);
set(handles.popupmenu_influx,'Value',1);
set(handles.popupmenu_influx_direction,'Value',1);
set(handles.edit_time,'String','');
set(handles.edit_time,'Userdata',[]);
set(handles.edit_timestep,'String','');
set(handles.edit_timestep,'Userdata',[]);
set(handles.edit_xinterv,'String','');
set(handles.edit_xinterv,'Userdata',[]);
set(handles.edit_p,'String','');
set(handles.edit_p,'Userdata',[]);
set(handles.edit_cutfracK,'String','');
set(handles.edit_cutfracK,'Userdata',[]);
set(handles.edit_FBL,'String','');
set(handles.edit_FBL,'Userdata',[]);
set(handles.edit_Nonlin_Adj,'Userdata',[]);
set(handles.edit_Nonlin_Adj,'string','');
set(handles.edit_depth_F,'Userdata',[]);
set(handles.edit_depth_F,'String','');
set(handles.pushbutton_load_file_bathy,'Userdata',[]);
set(handles.edit_slope,'Userdata',[]);
set(handles.edit_slope,'String','');
set(handles.edit_Dmidref,'enable','off','userdata',[],'string','');
set(handles.text118,'enable','off');
set(handles.text119,'enable','off');
set(handles.checkbox_bottom_friction,'value',0);
set(handles.togglebutton_friction,'Visible','off')
set(handles.togglebutton_friction,'userdata',[])
set(handles.popupmenu_bathy,'Value',1);
set(handles.checkbox_wall,'Value',0);
set(handles.edit_Xwall,'Userdata',[]);
set(handles.edit_Xwall,'String','');
set(handles.edit_Refl_percent,'Userdata',[]);
set(handles.edit_Refl_percent,'String','');
set(handles.checkbox_interior,'value',0);
set(handles.edit_Interior_Time,'Userdata',[]);
set(handles.edit_Interior_Time,'String','');
set(handles.edit_Interior_Time,'enable','off');
set(handles.edit_dtIP,'Userdata',[]);
set(handles.edit_dtIP,'String','');
set(handles.edit_dtIP,'enable','off');
set(handles.edit_savename,'String','');
set(handles.edit_note,'String','');
set(handles.edit_partition,'userdata',[]);
set(handles.edit_partition,'string','');
set(handles.edit_partition,'enable','off');
set(handles.checkbox_partition,'value',1);
set(handles.checkbox_combineFiles,'value',1);
set(handles.checkbox_combineFiles,'enable','off');
set(handles.checkbox_ramp,'value',0);
set(handles.edit_ramp,'userdata',4);
set(handles.edit_ramp,'string','4');
set(handles.edit_ramp,'visible','off');
set(handles.text100,'visible','off');
set(handles.text101,'visible','off');

set(handles.checkbox17,'value',0);
set(handles.togglebutton_bdy_assim,'enable','off')
set(handles.togglebutton_bdy_assim,'userdata',[])

set(handles.checkbox_wind,'value',0);
set(handles.edit_coef_wind,'enable','off');
set(handles.edit_wind_tinterv,'enable','off');
set(handles.edit_wind_xinterv,'enable','off');

set(handles.checkbox_ship,'value',0);
set(handles.pushbutton_ship_setupGui,'userdata',[],'enable','off');

Hs  = findobj('Tag','edit_Hs');            %%
Tp  = findobj('Tag','edit_Tp');            %%
Jsg = findobj('Tag','edit_Jsg');            %%
WT_I_func = findobj('Tag','pushbutton_load_file_IVP');            %%
WT_O_data = findobj('Tag','WT_Other_data');            %%
set(Hs,'enable','off'); set(Tp,'enable','off');set(Jsg,'enable','off');
set(WT_I_func,'Visible','off');set(WT_O_data,'Visible','off');

Disp_O       =findobj('Tag','togglebutton_InputDispersion');
Group_vel_O  =findobj('Tag','togglebutton_InputGroupVel');
set(Disp_O,'Visible','off');  set(Group_vel_O,'Visible','off');

checkbox_def_param=findobj('Tag','checkbox_def_paramB');
Break_param=findobj('Tag','Breaking_param');
Break_param2=findobj('Tag','edit_paramB_UF_T');

set(Break_param,'Visible','off');
set(checkbox_def_param,'Visible','off');
set(handles.text72,'Visible','off');
set(handles.Breaking_param_text,'Visible','off');
if get(handles.checkbox_def_paramB,'Value')==1
    set(Break_param2,'Enable','off');
    set(Break_param2,'Visible','off');
else
    set(Break_param2,'Enable','on');
    set(Break_param2,'Visible','on');
end
set(Break_param2,'Visible','off');

depth_F      =findobj('Tag','edit_depth_F');
depth_O      =findobj('Tag','pushbutton_load_file_bathy');
Slope     =findobj('Tag','edit_slope');
set(depth_F ,'tooltipString','Input depth');
set(depth_F ,'Enable','on');
set(depth_O ,'Visible','off');
set(depth_O ,'Enable','on');
set(Slope ,'Enable','off');

set(handles.edit_Xwall ,'Visible','off');
set(handles.popupmenu_WallType,'Visible','off');
set(handles.edit_Refl_percent,'Visible','off');
set(handles.togglebutton_inputEq,'Visible','off');
set(handles.text68,'Visible','off');
set(handles.text83,'Visible','off');
set(handles.text84,'Visible','off');
set(handles.text86,'Visible','off');


set(handles.pushbutton1,'enable','on')
set(handles.pushbutton_run,'enable','off')
set(handles.togglebutton_logfile,'enable','off')

set(handles.popupmenu_IVP,'Value',1);
set(handles.edit_IVP_A,'enable','off');
set(handles.edit_IVP_A,'userdata',[]);
set(handles.edit_IVP_lambda,'enable','off');
set(handles.edit_IVP_lambda,'userdata',[]);
set(handles.edit_IVP_x0,'enable','off');
set(handles.edit_IVP_x0,'userdata',[]);

set(handles.pushbutton_load_file_IVP,'visible','off')
set(handles.pushbutton_load_file_IVP,'userdata',[])

set(handles.togglebutton_inputEq,'userdata','1-min(f.^2,0.5)');
DispEq='sign(k).*sqrt(9.81*k.*tanh(d.*k))';
set(handles.togglebutton_InputDispersion,'userdata',DispEq);
UgEq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)';
set(handles.togglebutton_InputGroupVel,'userdata',UgEq);


if Flag==1
    statusbarObj.setText('Clear.');
end

function save_handles(handles,ProjectDir)
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');

GUImain.proj.path=handles.input.pathnow;
GUImain.proj.module =handles.input.module;
GUImain.proj.projhist=handles.input.projecthist;
Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
GUImain.proj.projdir=Curdir{get(handles.popupmenu_CurDir,'Value')};
GUImain.proj.name=get(handles.edit_savename,'String');
GUImain.proj.workdir=[GUImain.proj.projdir,'\',GUImain.proj.name,'\'];
GUImain.proj.note=get(handles.edit_note,'String');

GUImain.IVP.type=get(handles.popupmenu_IVP,'value');
IVPname         =cellstr(get(handles.popupmenu_IVP,'String'));
GUImain.IVP.typename=IVPname{get(handles.popupmenu_IVP,'value')};
GUImain.IVP.file=get(handles.pushbutton_load_file_IVP,'Userdata');
GUImain.IVP.A=get(handles.edit_IVP_A,'Userdata');
GUImain.IVP.lambda=get(handles.edit_IVP_lambda,'Userdata');
GUImain.IVP.x0=get(handles.edit_IVP_x0,'Userdata');

wavetype=cellstr(get(handles.WaveType,'String'));
GUImain.WT=wavetype{get(handles.WaveType,'Value')};

if get(handles.WaveType,'Value')==2
    GUImain.Hs= get(handles.edit_Hs,'Userdata')*2; % harmonic input is amplitude
elseif get(handles.WaveType,'Value')==3
    GUImain.Hs= get(handles.edit_Hs,'Userdata');
end

GUImain.Tp= get(handles.edit_Tp,'Userdata');
GUImain.Jsg= get(handles.edit_Jsg,'Userdata');
GUImain.WT_O_data=get(handles.WT_Other_data,'Userdata');

GUImain.Filter= get(handles.checkbox_Filter,'Value');
GUImain.Filter_LFHF=get(handles.edit_LFHF,'Userdata');
GUImain.ramp.check=get(handles.checkbox_ramp,'Value');
GUImain.ramp.val=get(handles.edit_ramp,'userdata');


modeltype=cellstr(get(handles.dynModel,'String'));
GUImain.evol=modeltype{get(handles.dynModel,'Value')};
GUImain.Break= get(handles.checkbox_Breaking,'Value');
GUImain.Break_param= get(handles.Breaking_param,'Userdata');
GUImain.BreakDef= get(handles.checkbox_def_paramB,'Value');
GUImain.Break_paramDef= get(handles.edit_paramB_UF_T,'Userdata');

disptype=cellstr(get(handles.popupmenu_dispersion,'String'));
GUImain.disp=disptype{get(handles.popupmenu_dispersion,'Value')};
GUImain.disp_O_file=get(handles.togglebutton_InputDispersion,'Userdata');
GUImain.groupVel_O_file=get(handles.togglebutton_InputGroupVel,'Userdata');
Influxtype=cellstr(get(handles.popupmenu_influx,'String'));
GUImain.inftype=Influxtype{get(handles.popupmenu_influx,'Value')};
Influxdir=cellstr(get(handles.popupmenu_influx_direction,'String'));
GUImain.infdir=Influxdir{get(handles.popupmenu_influx_direction,'Value')};

GUImain.tinterv=get(handles.edit_time,'Userdata');
GUImain.tstep=get(handles.edit_timestep,'Userdata');
GUImain.xinterv=get(handles.edit_xinterv,'Userdata');
GUImain.Xinflux=get(handles.edit_Xinflux,'Userdata');
GUImain.dx=get(handles.edit_p,'Userdata');
GUImain.cutfracK=get(handles.edit_cutfracK,'Userdata');
GUImain.FBL=get(handles.edit_FBL,'Userdata');
GUImain.NonlinAdj=get(handles.edit_Nonlin_Adj,'Userdata');
bathytype=cellstr(get(handles.popupmenu_bathy,'String'));
GUImain.bathy=bathytype{get(handles.popupmenu_bathy,'Value')};

GUImain.depth=get(handles.edit_depth_F,'Userdata');
GUImain.slope=get(handles.edit_slope,'Userdata');
GUImain.bathyO=get(handles.pushbutton_load_file_bathy,'Userdata');
GUImain.bottomfriction_check=get(handles.checkbox_bottom_friction,'Value');
GUImain.frictiondata=get(handles.togglebutton_friction,'Userdata');

GUImain.wind.check=get(handles.checkbox_wind,'value');
GUImain.wind.coef=get(handles.edit_coef_wind,'userdata');
GUImain.wind.tinterv=get(handles.edit_wind_tinterv,'userdata');
GUImain.wind.xinterv=get(handles.edit_wind_xinterv,'userdata');

GUImain.Wall= get(handles.checkbox_wall,'Value');
GUImain.Wall_position= get(handles.edit_Xwall,'Userdata');
GUImain.Wall_Coef= get(handles.edit_Refl_percent,'Userdata');
GUImain.Wall_method=get(handles.popupmenu_wall_method,'value');
GUImain.Wall_type= get(handles.popupmenu_WallType,'Value');
GUImain.Wall_reflEq=get(handles.togglebutton_inputEq,'Userdata');
GUImain.InteriorProp=get(handles.checkbox_interior,'Value');
GUImain.IP_time=get(handles.edit_Interior_Time,'Userdata');
GUImain.IP_dt=get(handles.edit_dtIP,'Userdata');
GUImain.partition=get(handles.edit_partition,'userdata');
GUImain.partitionId=get(handles.checkbox_partition,'value');
GUImain.partitionCombines=get(handles.checkbox_combineFiles,'Value');
GUImain.shippar=get(handles.pushbutton_ship_setupGui,'userdata');
GUImain.shippar.check=get(handles.checkbox_ship,'value');

save_state(GUImain);
statusbarObj.setText(['Project saved as abproj_',GUImain.proj.name,'.mat']);


% --- Executes on button press in checkbox_partition.
function checkbox_partition_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_partition
Id=get(hObject,'Value');
if Id==0
    set(handles.edit_partition,'enable','on');
    set(handles.checkbox_combineFiles,'enable','on');
else
    set(handles.edit_partition,'enable','off');
    set(handles.checkbox_combineFiles,'enable','off');
end

function edit_partition_Callback(hObject, eventdata, handles)
% hObject    handle to edit_partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_partition as text
%        str2double(get(hObject,'String')) returns contents of edit_partition as a double
param=str2num(get(hObject,'String'));
set(handles.edit_partition,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_partition);
else %if correct
    statusbarObj.setText('');
end

if fix(param)~=param
    statusbarObj.setText('Error: not an integer input');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_partition);
elseif isinteger(param) %if correct
    statusbarObj.setText('');
end



% --- Executes during object creation, after setting all properties.
function edit_partition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Xinflux_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xinflux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xinflux as text
%        str2double(get(hObject,'String')) returns contents of edit_Xinflux as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Xinflux,'Userdata',param);

[statusbarObj]=JavaFrame_handling();

if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Xinflux);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_Xinflux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xinflux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nonlin_Adj_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Nonlin_Adj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Nonlin_Adj as text
%        str2double(get(hObject,'String')) returns contents of edit_Nonlin_Adj as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Nonlin_Adj,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Nonlin_Adj);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_Nonlin_Adj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nonlin_Adj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_bottom_friction.
function checkbox_bottom_friction_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bottom_friction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_bottom_friction
Id=get(hObject,'Value');
if Id==1
    set(handles.togglebutton_friction,'visible','on');
else
    set(handles.togglebutton_friction,'visible','off');
    delete(input_friction_interval)
end

function edit_bottom_friction_coef_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bottom_friction_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bottom_friction_coef as text
%        str2double(get(hObject,'String')) returns contents of edit_bottom_friction_coef as a double

% --- Executes during object creation, after setting all properties.
function edit_bottom_friction_coef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bottom_friction_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_friction_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_friction_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_friction_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_friction_interval as a double

% --- Executes during object creation, after setting all properties.
function edit_friction_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_friction_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cutfracK_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cutfracK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cutfracK as text
%        str2double(get(hObject,'String')) returns contents of edit_cutfracK as a double
param=str2num(get(hObject,'String'));
set(handles.edit_cutfracK,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_cutfracK);
else %if correct
    statusbarObj.setText('');
end

if fix(param)~=param
    statusbarObj.setText('Error: not an integer input.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
end

if param<0
    statusbarObj.setText('Error: negative value.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_cutfracK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cutfracK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Hs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Hs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Hs as text
%        str2double(get(hObject,'String')) returns contents of edit_Hs as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Hs,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Hs);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Hs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Hs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tp as text
%        str2double(get(hObject,'String')) returns contents of edit_Tp as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Tp,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Tp);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Jsg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Jsg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Jsg as text
%        str2double(get(hObject,'String')) returns contents of edit_Jsg as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Jsg,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Jsg);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Jsg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Jsg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dtIP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dtIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dtIP as text
%        str2double(get(hObject,'String')) returns contents of edit_dtIP as a double
param=str2num(get(hObject,'String'));
set(handles.edit_dtIP,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_dtIP);
else %if correct
    statusbarObj.setText('');
end

if fix(param)~=param
    statusbarObj.setText('Error: not an integer input');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_dtIP);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_dtIP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dtIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Refl_percent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Refl_percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Refl_percent as text
%        str2double(get(hObject,'String')) returns contents of edit_Refl_percent as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Refl_percent,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Refl_percent);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Refl_percent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Refl_percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Xwall_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xwall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xwall as text
%        str2double(get(hObject,'String')) returns contents of edit_Xwall as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Xwall,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Xwall);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Xwall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xwall (see GCBO)
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
set(handles.edit_slope,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_slope);
else %if correct
    statusbarObj.setText('');
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


% --- Executes on selection change in popupmenu_WallType.
function popupmenu_WallType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_WallType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_WallType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_WallType
Id=get(handles.popupmenu_WallType,'value');
global IDdemo

if IDdemo==1
    if Id==2
        Id=1;
        set(handles.popupmenu_WallType,'value',1)
        [statusbarObj]=JavaFrame_handling();
        statusbarObj.setText('Frequency dependent wall is not available in the demo version');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
    end
end

if Id==1
    set(handles.text84,'Visible','on');
    set(handles.edit_Refl_percent,'Visible','on');
    set(handles.togglebutton_inputEq,'Visible','off');
    set(handles.edit_Refl_percent,'tooltipString','Specify the reflection coefficient in [0,1]');
else
    set(handles.edit_Refl_percent,'Visible','off');
    set(handles.text84,'Visible','off');
    set(handles.togglebutton_inputEq,'Visible','on');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_WallType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_WallType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_def_porosity.
function checkbox_def_porosity_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_def_porosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_def_porosity


% --- Executes on button press in pushbutton_freq_porosity_file.
function pushbutton_freq_porosity_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_freq_porosity_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in popupmenu_IVP.
function popupmenu_IVP_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_IVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_IVP contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_IVP
Id=get(handles.popupmenu_IVP,'value');
if Id==1
    set(handles.edit_IVP_A,'enable','off');
    set(handles.edit_IVP_lambda,'enable','off');
    set(handles.edit_IVP_x0,'enable','off');
    set(handles.pushbutton_load_file_IVP,'visible','off')
elseif Id==2 ||Id==3
    set(handles.edit_IVP_A,'enable','on');
    set(handles.edit_IVP_lambda,'enable','on');
    set(handles.edit_IVP_x0,'enable','on');
    set(handles.pushbutton_load_file_IVP,'visible','off')
elseif Id==4
    set(handles.edit_IVP_A,'enable','off');
    set(handles.edit_IVP_lambda,'enable','off');
    set(handles.edit_IVP_x0,'enable','off');
    set(handles.pushbutton_load_file_IVP,'visible','on')
end

% --- Executes during object creation, after setting all properties.
function popupmenu_IVP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_IVP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IVP_A_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IVP_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IVP_A as text
%        str2double(get(hObject,'String')) returns contents of edit_IVP_A as a double
param=str2num(get(hObject,'String'));
set(handles.edit_IVP_A,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_IVP_A);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_IVP_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IVP_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IVP_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IVP_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IVP_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_IVP_lambda as a double
param=str2num(get(hObject,'String'));
set(handles.edit_IVP_lambda,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_IVP_lambda);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_IVP_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IVP_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IVP_x0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IVP_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IVP_x0 as text
%        str2double(get(hObject,'String')) returns contents of edit_IVP_x0 as a double
param=str2num(get(hObject,'String'));
set(handles.edit_IVP_x0,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_IVP_x0);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_IVP_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IVP_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IVP_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IVP_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IVP_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_IVP_mu as a double
param=str2num(get(hObject,'String'));
set(handles.edit_IVP_mu,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_IVP_mu);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_IVP_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IVP_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ramp.
function checkbox_ramp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ramp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ramp
Id=get(hObject,'Value');

if Id==1
    set(handles.edit_ramp,'visible','on')
    set(handles.text100,'visible','on')
    set(handles.text101,'visible','on')
else
    set(handles.edit_ramp,'visible','off')
    set(handles.text100,'visible','off')
    set(handles.text101,'visible','off')
end


function edit_ramp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ramp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ramp as text
%        str2double(get(hObject,'String')) returns contents of edit_ramp as a double
param=str2num(get(hObject,'String'));
set(handles.edit_ramp,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_ramp);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_ramp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ramp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_combineFiles.
function checkbox_combineFiles_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_combineFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_combineFiles


% --- Executes on button press in pushbutton_inputEq.
function pushbutton_inputEq_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_inputEq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in togglebutton_inputEq.
function togglebutton_inputEq_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_inputEq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_inputEq
[statusbarObj]=JavaFrame_handling();

Id=get(hObject,'Value');

def.eq=get(handles.togglebutton_inputEq,'userdata');
if isempty(def.eq)
    def.eq='1-min(f.^2,0.5)';
end
def.hf=handles.togglebutton_inputEq;
if Id==1
    input_formula(def);
else
    [hf,yy]=  input_formula;
    set(handles.togglebutton_inputEq,'userdata',yy);
    delete(input_formula)
end
str=get(handles.togglebutton_inputEq,'userdata');
if isempty(strfind(str,'f'))
    statusbarObj.setText('Wrong format of the reflection coefficient. Use f as variable');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.togglebutton_inputEq);
    return;
else
    statusbarObj.setText('');
end



% --- Executes on button press in togglebutton_InputDispersion.
function togglebutton_InputDispersion_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_InputDispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_InputDispersion
[statusbarObj]=JavaFrame_handling();
Id=get(hObject,'Value');

def.eq=get(handles.togglebutton_InputDispersion,'userdata');
if isempty(def.eq)
    def.eq='sign(k).*sqrt(9.81*k.*tanh(d.*k))';
end

def.hf=handles.togglebutton_InputDispersion;

if Id==1
    input_formula_Disp(def);
else
    [hf,yy]=  input_formula_Disp;
    set(handles.togglebutton_InputDispersion,'userdata',yy);
    delete(input_formula_Disp)
end

str=get(handles.togglebutton_InputDispersion,'userdata');
if isempty(strfind(str,'k'))
    statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.togglebutton_InputDispersion);
    return;
else
    statusbarObj.setText('');
end
if isempty(strfind(str,'d'))
    statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.togglebutton_InputDispersion);
    return;
else
    statusbarObj.setText('');
end


% --- Executes on button press in togglebutton_InputGroupVel.
function togglebutton_InputGroupVel_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_InputGroupVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_InputGroupVel
[statusbarObj]=JavaFrame_handling();
Id=get(hObject,'Value');

def.eq=get(handles.togglebutton_InputGroupVel,'userdata');
if isempty(def.eq)
    def.eq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)';
end

def.hf=handles.togglebutton_InputGroupVel;
if Id==1
    input_formula_Vg(def);
else
    [hf,yy]=  input_formula_Vg;
    set(handles.togglebutton_InputGroupVel,'userdata',yy);
    delete(input_formula_Vg)
end
str=get(handles.togglebutton_InputGroupVel,'userdata');
if isempty(strfind(str,'k'))
    statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.togglebutton_InputGroupVel);
    return;
else
    statusbarObj.setText('');
end
if isempty(strfind(str,'d'))
    statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.togglebutton_InputGroupVel);
    return;
else
    statusbarObj.setText('');
end


% --- Executes on button press in togglebutton_friction.
function togglebutton_friction_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_friction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_friction
[statusbarObj]=JavaFrame_handling();
Id=get(hObject,'Value');

def.str=get(handles.togglebutton_friction,'userdata');
assignin('base','str',def.str)
def.hf=handles.togglebutton_friction;
if Id==1
    input_friction_interval(def);
else
    [hf,yy]=  input_friction_interval;
    set(handles.togglebutton_friction,'userdata',yy);
    delete(input_friction_interval)
end






% --------------------------------------------------------------------
function Activation_Callback(hObject, eventdata, handles)
% hObject    handle to Activation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
About_AB;


% --- Executes on button press in checkbox_wind.
function checkbox_wind_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_wind
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_coef_wind,'enable','on');
    set(handles.edit_wind_tinterv,'enable','on');
    set(handles.edit_wind_xinterv,'enable','on');
else
    set(handles.edit_coef_wind,'enable','off');
    set(handles.edit_coef_wind,'string','');
    set(handles.edit_wind_tinterv,'enable','off','string','');
    set(handles.edit_wind_xinterv,'enable','off','string','');
end


function edit_coef_wind_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coef_wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coef_wind as text
%        str2double(get(hObject,'String')) returns contents of edit_coef_wind as a double
param=str2num(get(hObject,'String'));
set(handles.edit_coef_wind,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_coef_wind);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit_coef_wind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coef_wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ship.
function checkbox_ship_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ship (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ship
Id=get(hObject,'Value');
if Id==1
    set(handles.pushbutton_ship_setupGui,'enable','on')
else
    set(handles.pushbutton_ship_setupGui,'enable','off')
end

% --- Executes on button press in pushbutton_ship_setupGui.
function pushbutton_ship_setupGui_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ship_setupGui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Curdir   =cellstr(get(handles.popupmenu_CurDir,'String'));
path=Curdir{get(handles.popupmenu_CurDir,'Value')};
datainput=get(hObject,'userdata');
[statusbarObj]=JavaFrame_handling();
if ~isfield(datainput,'FlagCloseGui')
    datainput.FlagCloseGui=0;
end
if datainput.FlagCloseGui==0
    statusbarObj.setText('');
    datainput.FlagCloseGui=1;
    set(hObject,'userdata',datainput);
    shipparam=GUI_Ship_Setup(path,datainput);
    set(hObject,'userdata',shipparam);
else
    statusbarObj.setText('Ship setup panel is already opened');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
end



function edit_wind_tinterv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wind_tinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wind_tinterv as text
%        str2double(get(hObject,'String')) returns contents of edit_wind_tinterv as a double
param=str2num(get(hObject,'String'));
set(hObject,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_coef_wind);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_wind_tinterv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wind_tinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_wind_xinterv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wind_xinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wind_xinterv as text
%        str2double(get(hObject,'String')) returns contents of edit_wind_xinterv as a double
param=str2num(get(hObject,'String'));
set(hObject,'Userdata',param);
[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_coef_wind);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_wind_xinterv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wind_xinterv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_depth_reff.
function popupmenu_depth_reff_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_depth_reff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_depth_reff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_depth_reff
[statusbarObj]=JavaFrame_handling();
Id=get(handles.popupmenu_bathy,'value');
if Id==3
    set(hObject,'value',2);
    statusbarObj.setText('Only 3 depths reference for run-up case');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.popupmenu_depth_reff);
elseif Id==4
    bath=get(handles.pushbutton_load_file_bathy,'userdata');
    if ~isempty(bath.data)
        if any(bath.data(:,2)>0)
            set(hObject,'value',2);
            statusbarObj.setText('Only 3 depths reference for run-up case');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.popupmenu_depth_reff);
        end
    end
else
    statusbarObj.setText('');
end
IdD=get(hObject,'value');
if IdD==2
    set(handles.edit_Dmidref,'enable','on');
    set(handles.text118,'enable','on');
    set(handles.text119,'enable','on');
else
    set(handles.edit_Dmidref,'enable','off');
    set(handles.text118,'enable','off');
    set(handles.text119,'enable','off');
end



% --- Executes during object creation, after setting all properties.
function popupmenu_depth_reff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_depth_reff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17
Id=get(hObject,'Value');
if Id==1
    set(handles.togglebutton_bdy_assim,'visible','on');
else
    set(handles.togglebutton_bdy_assim,'visible','off');
    delete(input_bdyassim1d)
end


% --- Executes on button press in togglebutton_bdy_assim.
function togglebutton_bdy_assim_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_bdy_assim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_bdy_assim
Id=get(hObject,'Value');

def.userdata=get(hObject,'userdata');
def.projdir=handles.input.projectdirectory;
def.hf=hObject;
if Id==1
    input_bdyassim1d(def);
else
    [hf,yy]=  input_bdyassim1d;
    set(hObject,'userdata',yy);
    delete(input_bdyassim1d)
end





function edit_Dmidref_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Dmidref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Dmidref as text
%        str2double(get(hObject,'String')) returns contents of edit_Dmidref as a double
param=str2num(get(hObject,'String'));
set(hObject,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit_Dmidref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Dmidref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_wall_method.
function popupmenu_wall_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_wall_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_wall_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_wall_method


% --- Executes during object creation, after setting all properties.
function popupmenu_wall_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_wall_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
