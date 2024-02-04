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

function varargout = PostProc_GUI(varargin)
% POSTPROC_GUI MATLAB code for PostProc_GUI.fig
%      POSTPROC_GUI, by itself, creates a new POSTPROC_GUI or raises the existing
%      singleton*.
%
%      H = POSTPROC_GUI returns the handle to a new POSTPROC_GUI or the handle to
%      the existing singleton*.
%
%      POSTPROC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTPROC_GUI.M with the given input arguments.
%
%      POSTPROC_GUI('Property','Value',...) creates a new POSTPROC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PostProc_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PostProc_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PostProc_GUI

% Last Modified by GUIDE v2.5 29-Jan-2020 14:40:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PostProc_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @PostProc_GUI_OutputFcn, ...
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


% --- Executes just before PostProc_GUI is made visible.
function PostProc_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PostProc_GUI (see VARARGIN)

% Choose default command line output for PostProc_GUI

%reset_handles(handles,0);


global IDdemoPP

handles.output = hObject;
handles.input=varargin;

%set(hObject,'units','normalized','Position',[0.6 0.3 0.31 0.48])
% set(0,'ScreenPixelsPerInch',96)

pathstr=handles.input{1}.Proj.path;


%javaFrame = get(hObject,'JavaFrame');
try   
    axes(handles.Logo1);
    if ~isdeployed
        imshow([pathstr,'\Toolbox\lib1d\Misc\logo\hawassi_black.jpg']); %
        %   javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
    else
        imshow('hawassi_black.jpg');
        %   javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'))
        
    end
catch
end
    %  warning('off','all');

IDdemoPP=0;
if IDdemoPP==0
set(handles.figure1,'Name',['HAWASSI AB1: Post-Processing']);
else
set(handles.figure1,'Name',['HAWASSI AB1 demo: Post-Processing']);    
end
% set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
% set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')

set(handles.popupmenu1,'Value',1)
set(handles.pushbutton_load_simul,'Enable','on'); 
set(handles.pushbutton_load_simul,'Userdata',[]);

if ~isfield(handles.input{1},'time')
    additemtolistbox(handles.popupmenu_CurDir,handles.input{1}.Proj.projdir);
elseif isfield(handles.input{1},'time')
    time=handles.input{1}.time;
    eta=handles.input{1}.eta;
    x=handles.input{1}.x;
    bathy=handles.input{1}.bathy;
    projdir=handles.input{1}.Proj.projdir;
    workdir=handles.input{1}.Proj.Dir;
    Xinflux=handles.input{1}.Xinflux;
    minbath=-min(bathy);
    IndtI=1;
    IndtF=length(time);
    tcoarse=1;
    Fbdy=handles.input{1}.Fbdy;
    additemtolistbox(handles.popupmenu_CurDir,projdir);
    
    set(handles.popupmenu1,'Value',2);
    set(handles.pushbutton_load_simul,'enable','off');
    
    
    PP_X=[ceil(x(1));floor(x(end))];
    
    set(handles.edit_Profile_Position,'Userdata',PP_X);
    set(handles.edit_Profile_Position,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
    set(handles.edit8,'Userdata',PP_X);
    set(handles.edit8,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
    
    PP_T=[ceil(time(IndtI));floor(time(IndtF))];
    set(handles.edit_Signal_time,'Userdata',PP_T);
    set(handles.edit_Signal_time,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
    set(handles.edit7,'Userdata',PP_T);
    set(handles.edit7,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
    set(handles.edit16,'Userdata',PP_T);
    set(handles.edit16,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
    
    
    savename=handles.input{1}.savename;
    set(handles.edit_projectname,'string',savename);
end

set(handles.edit_scale,'Userdata',1);
set(handles.edit_scale,'String',[num2str(1),'/',num2str(1)]);
set(handles.edit24,'Userdata',1);
set(handles.edit24,'String',[num2str(1),'/',num2str(1)]);
set(handles.edit23,'Userdata',1);
set(handles.edit23,'String',num2str(1));
set(handles.edit11,'Userdata',1);
set(handles.edit11,'String',num2str(1));
set(handles.edit_GIF_setting,'Userdata',[0.01;inf]);
set(handles.edit_GIF_setting,'String','1/100;inf');


Data_Sim_load =findobj('Tag','pushbutton_load_simul');
%text_filter=findobj('Tag','text11');
edit_filter=findobj('Tag','edit12');
my_data.Load=0;

set(handles.pushbutton_load_simul,'Userdata',my_data);

SimData=cellstr(get(handles.popupmenu1,'String'));
SimData_pp=SimData{get(handles.popupmenu1,'Value')};

if strcmp(SimData_pp,'Other')==0
    set(Data_Sim_load,'Enable','off');
else
    set(Data_Sim_load,'Enable','on'); 
end

%set(text_filter,'Enable','off');
set(edit_filter,'Visible','off');
set(handles.text11,'Visible','off');
set(handles.text33,'Visible','off');


%%%Init condition of Radio Button
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);
set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton1,'Enable','On');
set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton3,'Enable','On');
set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton4,'Enable','On');
set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'off')
%%%%
set(handles.pushbutton2,'visible','on')
set(handles.pushbutton7,'visible','on')
set(handles.pushbutton3,'Visible','On');

%%%%'tooltipString'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.popupmenu_CurDir,'tooltipString','The directory list')
set(handles.pushbutton_browse_Dir,'tooltipString','Specify the directory')
set(handles.popupmenu1,'tooltipString','Data status: Other: data need to be provided, please load.')
set(handles.checkbox1,'tooltipString','Option to restrict the frequency band of the time signals')
set(handles.edit12,'tooltipString','Provides the desired frequency band. Input format [Low Freq; High Freq] ')
set(handles.radiobutton1,'tooltipString','Option to show the wave profile(s) and/or the time trace(s)')
set(handles.checkbox_Profile,'tooltipString','Option to show the wave profile(s)')
set(handles.edit_Profile_time,'tooltipString','Provide the time at t: t1 or at t1;t2;t3')
set(handles.edit_Profile_Position,'tooltipString','Specify the horizontal spatial interval. [x_start;x_end]')
set(handles.checkbox_MTA,'tooltipString','Option to show Maximum Temporal Amplitude (MTA)')
set(handles.checkbox10,'tooltipString','Option to show bathymetry')
set(handles.edit_scale,'tooltipString','Show the bathymetry in a scale, eg. 1 for full scale, or 1/20 in a scale')
set(handles.checkbox_signal,'tooltipString','Option to show the time trace(s) or spectrum.')
set(handles.edit_Signal_Position,'tooltipString','Provide the buoy(s) position at x: x1 or x1;x2;x3')
set(handles.edit_Signal_time,'tooltipString','Specify the time interval. [t_start;t_end]')
set(handles.checkbox12,'tooltipString','Option to show the time trace(s)')
set(handles.checkbox13,'tooltipString','Option to show the spectrum')
set(handles.checkbox_Ampl_Sp,'tooltipString','Option to show the Amplitude spectrum')
set(handles.checkbox_Energy_Sp,'tooltipString','Option to show the Energy spectrum')
set(handles.edit23,'tooltipString','Provide the smoothing factor of the spectrum. 1: No smoother applied, >1 smoother applied.')

set(handles.radiobutton4,'tooltipString','Option to show the animation of the simulation')
set(handles.edit7,'tooltipString','Specify the time interval. [t_start;t_end]')
set(handles.edit8,'tooltipString','Specify the horizontal spatial interval. [x_start;x_end]')
set(handles.checkbox5,'tooltipString','Option to show Maximum Temporal Amplitude (MTA)')
set(handles.checkbox17,'tooltipString','Option to show bathymetry')
set(handles.edit24,'tooltipString','Show the bathymetry in a scale, eg. 1 for full scale, or 1/20 in a scale')
set(handles.checkbox_GIF,'tooltipString','Option to set the output .gif file')
set(handles.edit_GIF_setting,'tooltipString','Adjust the delay time and number of loop. Input format [delay_time;#loop], default: [1/100;inf]')

set(handles.radiobutton3,'tooltipString','Option to show comparison between the simulation and the measurement')
set(handles.pushbutton3,'tooltipString','Provide the measurement data');
set(handles.edit10,'tooltipString','This is the buoy(s) position that loaded from the file');
set(handles.edit16,'tooltipString','Specify the time interval. [t_start;t_end]')
set(handles.checkbox32,'tooltipString','Option to shift the simulation in order to get the best correlation with the measurement')
set(handles.checkbox33,'tooltipString','Automatic shift to get the best correlation.')
set(handles.edit34,'tooltipString','Manual shift. Provide the shifting factor. Input format: a *dt, a in [1,2,3, ...].')
set(handles.checkbox16,'tooltipString','Option to show the time trace(s)')
set(handles.checkbox_Val_Spectra,'tooltipString','Option to show the spectra')
set(handles.checkbox_Amplitude_SP,'tooltipString','Option to show the Amplitude spectrum')
set(handles.checkbox__P3_energySp,'tooltipString','Option to show the Energy spectrum')
set(handles.edit11,'tooltipString','Provide the smoothing factor of the spectrum. 1: No smoother applied, >1 smoother applied.')
set(handles.checkbox_Val_quant,'tooltipString','Option to show the quantitative results, e.g. Correlation, Variance, Kurtosis , Significant Wave height ')
set(handles.checkbox_skew_assym,'tooltipString','Option to show the plot of the significant Height (Hs), skewness (Sk), Asymmetry (As) and Kurtosis (K)')
set(handles.checkbox_MTA_MWL,'tooltipString','Option to show the plot of Maximum Temporal Amplitude (MTA) and Mean Water Level (MWL)')

set(handles.checkbox_xlim,'tooltipString','Option to adjust the horizontal interval')
set(handles.edit28,'tooltipString','Specify the horizontal interval. [x_start;x_end]')
set(handles.checkbox28,'tooltipString','Option to adjust the vertical interval')
set(handles.edit29,'tooltipString','Specify the vertical interval. [z_start;z_end]')
set(handles.checkbox29,'tooltipString','Option to coarse the data in spatial')
set(handles.edit30,'tooltipString','Specify the coarsing factor. [1,2,...]')
set(handles.checkbox30,'tooltipString','Option to coarse the data in time')
set(handles.edit31,'tooltipString','Specify the coarsing factor. [1,2,...]')
set(handles.checkbox_MTA_time_interval,'tooltipString','Option to adjust the time interval for the MTA plot')
set(handles.edit_MTA_time,'tooltipString','Specify the time interval. [t_start;t_end]')
set(handles.checkbox35,'tooltipString','Option to save the figures.')
set(handles.popupmenu5,'tooltipString','Choose an output file type.')
set(handles.checkbox34,'tooltipString','Option to save the postprocessing data.')

set(handles.pushbutton2,'tooltipString','Start to process.')
set(handles.pushbutton7,'tooltipString','Closes all active figures.')

% Update handles structure
guidata(hObject, handles);

global PPdata
PPdata=[];
PPdata.workdir='';

% UIWAIT makes PostProc_GUI wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PostProc_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1} =[];


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
param=str2num(get(hObject,'String'));
set(handles.edit3,'Userdata',param);


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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
param=str2num(get(hObject,'String'));
set(handles.edit5,'Userdata',param);


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



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[pathstr,~,~] = fileparts(mfilename('fullpath'));
% if ~isdeployed
% cd (pathstr);
% addpath(genpath(pathstr));
% end
% GUIpp.path=pathstr;

[statusbarObj]=JavaFrame_handling();
CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
GUIpp.projdir=CurDir{get(handles.popupmenu_CurDir,'Value')};

SimData=cellstr(get(handles.popupmenu1,'String'));
GUIpp.SimData=SimData{get(handles.popupmenu1,'Value')};
GUIpp.MySim=get(handles.pushbutton_load_simul,'Userdata');

if strcmp(GUIpp.SimData,'Other')
    statusbarObj.setText('No simulation data, please load!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end


if isfield(GUIpp.MySim,'Load')
    Load_data_Id=GUIpp.MySim.Load;
else
    Load_data_Id=1;
end

if isempty(handles.input) ||Load_data_Id==1
    if isfield(GUIpp.MySim,'model') %just indicates from the new version of saved data
        GUIpp.simul_file=[GUIpp.MySim.Proj.savename,'_simul.mat'];
        GUIpp.savename  =GUIpp.MySim.Proj.savename;
    else
        GUIpp.simul_file=[GUIpp.MySim.savename,'_simul.mat'];
        GUIpp.savename  =GUIpp.MySim.savename;
    end
    timesimul=GUIpp.MySim.output.time;
    xref=[round(GUIpp.MySim.par.x(1)) round(GUIpp.MySim.par.x(end))];
    tref=[floor(timesimul(1)) round(timesimul(end))];
else
    GUIpp.savename=handles.input{1}.savename;
    GUIpp.simul_file=[GUIpp.savename,'_simul.mat'];
    
    timesimul=handles.input{1}.time;
    xsimul=handles.input{1}.x;
    xref=[round(xsimul(1)) round(xsimul(end))];
    tref=[floor(timesimul(1)) round(timesimul(end))];
end

GUIpp.Filter=get(handles.checkbox1,'Value');
GUIpp.LFHF=get(handles.edit12,'Userdata');


GUIpp.PP1=get(handles.radiobutton1,'Value');
GUIpp.PP1_Prof=get(handles.checkbox_Profile,'Value');
GUIpp.PP1_Prof_T=get(handles.edit_Profile_time,'Userdata');
GUIpp.PP1_Prof_X=get(handles.edit_Profile_Position,'Userdata');
GUIpp.PP1_Prof_MTA=get(handles.checkbox_MTA,'Value');
GUIpp.PP1_Prof_Bathy=get(handles.checkbox10,'Value');
GUIpp.PP1_Prof_Bathy_scale=get(handles.edit_scale,'Userdata');

GUIpp.PP1_Buoy=get(handles.checkbox_signal,'Value');
GUIpp.PP1_Buoy_T=get(handles.edit_Signal_time,'Userdata');
GUIpp.PP1_Buoy_X=get(handles.edit_Signal_Position,'Userdata');
GUIpp.PP1_Buoy_Signal=get(handles.checkbox12,'Value');
GUIpp.PP1_Buoy_Spectrum=get(handles.checkbox13,'Value');
GUIpp.PP1_Buoy_Spectrum_Smooth=get(handles.edit23,'Userdata');
GUIpp.PP1_Buoy_Spectrum_Energy=get(handles.checkbox_Energy_Sp,'Value');
GUIpp.PP1_Buoy_Spectrum_Ampl=get(handles.checkbox_Ampl_Sp,'Value');
GUIpp.PP1_Buoy_exceedance=get(handles.checkbox_buoy_exceedance,'Value');
GUIpp.PP1_Energy=get(handles.checkbox_Energy,'Value');
GUIpp.PP1_Energy_in_x=get(handles.checkbox_Energy_in_x,'Value');
GUIpp.PP1_Spectra=get(handles.checkbox_PP1_spectra,'Value');
GUIpp.PP1_breakevents=get(handles.checkbox_breakevents,'Value');
GUIpp.PP1_HsAsSkK=get(handles.checkbox_PP1_Hs_Sk,'Value');
GUIpp.PP1_MTA_MWL=get(handles.checkbox_PP1_MTA_MWL,'Value');
GUIpp.PP1_Quant=get(handles.checkbox_PP1_buoy_quant,'Value');

GUIpp.PP1_kinematic=get(handles.checkbox_kinematic,'value');

GUIpp.PP1_Ship_XZThe=get(handles.checkbox_Ship_XZThe,'Value');
GUIpp.PP1_Ship_Forces=get(handles.checkbox_Ship_forces,'Value');

GUIpp.PP2=get(handles.radiobutton4,'Value');
GUIpp.PP2_T=get(handles.edit7,'Userdata');
GUIpp.PP2_X=get(handles.edit8,'Userdata');
GUIpp.PP2_MTA=get(handles.checkbox5,'Value');
GUIpp.PP2_Bathy=get(handles.checkbox17,'Value');
GUIpp.PP2_Bathy_scale=get(handles.edit24,'Userdata');
GUIpp.PP2_GIF=get(handles.checkbox_GIF,'Value');
GUIpp.PP2_GIF_setting=get(handles.edit_GIF_setting,'Userdata');

GUIpp.PP3=get(handles.radiobutton3,'Value');
GUIpp.PP3_X=get(handles.edit10,'Userdata');
GUIpp.PP3_T=get(handles.edit16,'Userdata');
GUIpp.PP3_Signals=get(handles.checkbox16,'Value');
GUIpp.PP3_Spectra=get(handles.checkbox_Val_Spectra,'Value');
GUIpp.PP3_Spectra_smooth=get(handles.edit11,'Userdata');
GUIpp.PP3_Spectra_Energy=get(handles.checkbox__P3_energySp,'Value');
GUIpp.PP3_Spectra_Ampl=get(handles.checkbox_Amplitude_SP,'Value');
GUIpp.PP3_Quant=get(handles.checkbox_Val_quant,'Value');
GUIpp.PP3_AsSk =get(handles.checkbox_skew_assym,'value');
GUIpp.PP3_MTAMWL=get(handles.checkbox_MTA_MWL,'value');
GUIpp.PP3MeasDat=get(handles.pushbutton3,'Userdata');
GUIpp.PP3_timeshift=get(handles.checkbox32,'value');
GUIpp.PP3_timeshiftBest=get(handles.checkbox33,'value');
GUIpp.PP3_timeshiftCustom=get(handles.edit34,'Userdata');
GUIpp.PPsetting_xlim=get(handles.checkbox_xlim,'value');
GUIpp.PPsetting_xlimVal=get(handles.edit28,'Userdata');
GUIpp.PPsetting_ylim=get(handles.checkbox28,'value');
GUIpp.PPsetting_ylimVal=get(handles.edit29,'Userdata');
GUIpp.PPsetting_xcoarse=get(handles.checkbox29,'value');
GUIpp.PPsetting_xcoarseVal=get(handles.edit30,'Userdata');
GUIpp.PPsetting_tcoarse=get(handles.checkbox30,'value');
GUIpp.PPsetting_tcoarseVal=get(handles.edit31,'Userdata');
GUIpp.PPsetting_MTA_time=get(handles.checkbox_MTA_time_interval,'value');
GUIpp.PPsetting_MTA_timeVal=get(handles.edit_MTA_time,'Userdata');
GUIpp.PPsetting_save_data=get(handles.checkbox34,'value');
GUIpp.save_figure=get(handles.checkbox35,'value');
GUIpp.save_animation=get(handles.checkbox_saveanimation,'value');

filetype=cellstr(get(handles.popupmenu5, 'string'));
GUIpp.save_figure_type=filetype{get(handles.popupmenu5,'value')};

GUIpp.workdir=[GUIpp.projdir,'\',GUIpp.savename,'\'];


%%%%Error Handling of User Input
%clc;
 statusbarObj.setText('');

if GUIpp.PP1==0 &&GUIpp.PP2==0 &&GUIpp.PP3==0
    statusbarObj.setText('Choose one of post-processing modes!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if GUIpp.PP1==1
    if GUIpp.PP1_Prof+ GUIpp.PP1_Buoy+GUIpp.PP1_Energy+GUIpp.PP1_Energy_in_x+GUIpp.PP1_Spectra+GUIpp.PP1_breakevents...
         +GUIpp.PP1_HsAsSkK+GUIpp.PP1_MTA_MWL+GUIpp.PP1_kinematic+GUIpp.PP1_Ship_XZThe...
         +GUIpp.PP1_Ship_Forces==0
        statusbarObj.setText('Choose an option!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if GUIpp.PP1_Prof==1
        if isempty(GUIpp.PP1_Prof_T)
            statusbarObj.setText('Input a single or multiple time!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_time);
            return;
        end
        
        if any(GUIpp.PP1_Prof_T<tref(1))
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(tref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_time);
            return;
        end
        
        if any(GUIpp.PP1_Prof_T>tref(2))
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(tref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_time);
            return;
        end
        
        if isempty(GUIpp.PP1_Prof_X)
            statusbarObj.setText('Input spatial interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if  length(GUIpp.PP1_Prof_X)~=2
            statusbarObj.setText('Wrong input format of the spatial interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(1)>GUIpp.PP1_Prof_X(2)
            statusbarObj.setText('Error:Wrong input format! [x_start;x_end]');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(1)<xref(1)
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(xref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(2)>xref(2)
            statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(xref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
    end
    
    if GUIpp.PP1_Buoy==1
        
        if isempty(GUIpp.PP1_Buoy_X)
            statusbarObj.setText('Input Buoy position(s)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_Position);
            return;
        end
        
        if any(GUIpp.PP1_Buoy_X<xref(1))
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(xref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_Position);
            return;
        end
        
        if any(GUIpp.PP1_Buoy_X>xref(2))
            statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(xref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_Position);
            return;
        end
        
        
        
        if isempty(GUIpp.PP1_Buoy_T)
            statusbarObj.setText('Input time interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if  length(GUIpp.PP1_Buoy_T)~=2
            statusbarObj.setText('Wrong input format of the time interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(1)>GUIpp.PP1_Buoy_T(2)
            statusbarObj.setText('Error:Wrong input format! [t_start;t_end]');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(1)<tref(1)
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(tref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(2)>tref(2)
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(tref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_Spectrum+GUIpp.PP1_Buoy_Signal+GUIpp.PP1_Quant==0
            statusbarObj.setText('Choose an option: signal, spectrum plot or quantitative info!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.checkbox_Ampl_Sp);
            return;
        end
        
        if GUIpp.PP1_Buoy_Spectrum==1
            if GUIpp.PP1_Buoy_Spectrum_Energy==0 && ...
                    GUIpp.PP1_Buoy_Spectrum_Ampl==0
                statusbarObj.setText('Choose a type of spectrum plot!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.checkbox_Ampl_Sp);
                return;
            end
            
            if isempty(GUIpp.PP1_Buoy_Spectrum_Smooth)
                statusbarObj.setText('Specify a smoother factor!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit23);
                return;
            end
            
            if length(GUIpp.PP1_Buoy_Spectrum_Smooth)~=1
                statusbarObj.setText('Wrong input format of the smoother factor!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit23);
                return;
            end
            
            if GUIpp.PP1_Buoy_Spectrum_Smooth<0
                statusbarObj.setText('Error: Negative smoother factor!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit23);
                return;
            end
        end
        
    end
    
    
    if GUIpp.PP1_Energy==1 || GUIpp.PP1_Spectra==1 ||GUIpp.PP1_Energy_in_x==1
        if isempty(GUIpp.PP1_Buoy_T)
            statusbarObj.setText('Input time interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if  length(GUIpp.PP1_Buoy_T)~=2
            statusbarObj.setText('Wrong input format of the time interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(1)>GUIpp.PP1_Buoy_T(2)
            statusbarObj.setText('Error:Wrong input format! [t_start;t_end]');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(1)<tref(1)
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(tref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
        if GUIpp.PP1_Buoy_T(2)>tref(2)
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(tref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Signal_time);
            return;
        end
        
         if isempty(GUIpp.PP1_Prof_X)
            statusbarObj.setText('Input spatial interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if  length(GUIpp.PP1_Prof_X)~=2
            statusbarObj.setText('Wrong input format of the spatial interval!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(1)>GUIpp.PP1_Prof_X(2)
            statusbarObj.setText('Error:Wrong input format! [x_start;x_end]');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(1)<xref(1)
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(xref(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        if GUIpp.PP1_Prof_X(2)>xref(2)
            statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(xref(2)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            uicontrol(handles.edit_Profile_Position);
            return;
        end
        
        
    end
    
elseif GUIpp.PP2==1
    if isempty(GUIpp.PP2_T)
        statusbarObj.setText('Input time interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit7);
        return;
    end
    
    if length(GUIpp.PP2_T)~=2
        statusbarObj.setText('Wrong input format of the time interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit7);
        return;
    end
    
    if GUIpp.PP2_T(1)>GUIpp.PP2_T(2)
        statusbarObj.setText('Error: Input value t_start>t_end');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit7);
        return;
    end
    
    if GUIpp.PP2_T(1)<tref(1)
        statusbarObj.setText(['Error: Input value t< minimum value (t=',num2str(tref(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit7);
        return;
    end
    
    if GUIpp.PP2_T(2)>tref(2)
        statusbarObj.setText(['Error: Input value t> maximum value (t=',num2str(tref(2)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit7);
        return;
    end
    
    if isempty(GUIpp.PP2_X)||length(GUIpp.PP2_X)~=2
        statusbarObj.setText('Input spatial interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit8);
        return;
    end
    
    
    if GUIpp.PP2_X(1)>GUIpp.PP2_X(2)
        statusbarObj.setText('Error: Wrong input format! [x_start;x_end]');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit8);
        return;
    end
    
    if any(GUIpp.PP2_X<xref(1))
        statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(xref(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit8);
        return;
    end
    
    if any(GUIpp.PP2_X>xref(2))
        statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(xref(2)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit8);
        return;
    end
    
    
elseif GUIpp.PP3==1
    if isempty(GUIpp.PP3MeasDat)
        statusbarObj.setText('load measurement data!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.pushbutton3);
        return;
    end
    
    if isempty(GUIpp.PP3_X)
        statusbarObj.setText('Input measurement positions!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit10);
        return;
    end
    
    if isempty(GUIpp.PP3_T)
        statusbarObj.setText('Input time interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit16);
        return;
    end
    
    if  length(GUIpp.PP3_T)~=2
        statusbarObj.setText('Wrong input format of the time interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit16);
        return;
    end
    
    if GUIpp.PP3_T(1)>GUIpp.PP3_T(2)
        statusbarObj.setText('Error:Wrong input format! [t_start;t_end]');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit16);
        return;
    end
    
    if GUIpp.PP3_T(1)<tref(1)
        statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(tref(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit16);
        return;
    end
    
    if GUIpp.PP3_T(2)>tref(2)
        statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(tref(2)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit16);
        return;
    end
    
    if GUIpp.PP3_timeshift==1
       if GUIpp.PP3_timeshiftBest==0
            if  isempty(GUIpp.PP3_timeshiftCustom)
                statusbarObj.setText('Input time step for shifting!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit34);
                return;
            end
             if  length(GUIpp.PP3_timeshiftCustom)~=1
                statusbarObj.setText('Wrong input format of the time step for shifting!');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                uicontrol(handles.edit34);
                return;
            end
       end
    end
    
    if GUIpp.PP3_AsSk+GUIpp.PP3_MTAMWL+GUIpp.PP3_Quant+GUIpp.PP3_Spectra+GUIpp.PP3_Signals==0
        statusbarObj.setText('Choose an option to process!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
        return;
    end
end

if GUIpp.PPsetting_xlim==1
    if isempty(GUIpp.PPsetting_xlimVal)
        statusbarObj.setText('Input horizontal axes limit');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit28);
        return;
    end
    
    if  length(GUIpp.PPsetting_xlimVal)~=2
        statusbarObj.setText('Wrong input format of the horizontal axes limit');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit28);
        return;
    end
    
    if GUIpp.PPsetting_xlimVal(1)>GUIpp.PPsetting_xlimVal(2)
        statusbarObj.setText('Error: Wrong input format! [x_start;x_end]');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit28);
        return;
    end
    
end
if GUIpp.PPsetting_ylim==1
    if isempty(GUIpp.PPsetting_ylimVal)
        statusbarObj.setText('Input vertical axes limit');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit29);
        return;
    end
    
    if  length(GUIpp.PPsetting_ylimVal)~=2
        statusbarObj.setText('Wrong input format of the vertical axes limit');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit29);
        return;
    end
    
    if GUIpp.PPsetting_ylimVal(1)>GUIpp.PPsetting_ylimVal(2)
        statusbarObj.setText('Error: Wrong input format! [z_bottom;z_top]');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit29);
        return;
    end
end

if GUIpp.PPsetting_xcoarse==1
    if isempty(GUIpp.PPsetting_xcoarseVal)
        statusbarObj.setText('Input x coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit30);
        return;
    end
    
    if  length(GUIpp.PPsetting_xcoarseVal)~=1
        statusbarObj.setText('Wrong input format of the x coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit30);
        return;
    end
    
    if fix(GUIpp.PPsetting_xcoarseVal)~=GUIpp.PPsetting_xcoarseVal
        statusbarObj.setText('Input an integer value of the x coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit30);
        return;
    end
    if GUIpp.PPsetting_xcoarseVal<0
        statusbarObj.setText('Error: negative x coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit30);
        return;
    end
end

if GUIpp.PPsetting_tcoarse==1
    if isempty(GUIpp.PPsetting_tcoarseVal)
        statusbarObj.setText('Input t coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit31);
        return;
    end
    if  length(GUIpp.PPsetting_tcoarseVal)~=1
        statusbarObj.setText('Wrong input format of the t coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit31);
        return;
    end
    
    if fix(GUIpp.PPsetting_tcoarseVal)~=GUIpp.PPsetting_tcoarseVal
        statusbarObj.setText('Input an integer value of the t coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit31);
        return;
    end
    if GUIpp.PPsetting_tcoarseVal<0
        statusbarObj.setText('Error: negative t coarse');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit31);
        return;
    end
end

if GUIpp.PPsetting_MTA_time==1
    if isempty(GUIpp.PPsetting_MTA_timeVal)
        statusbarObj.setText('Input MTA time interval');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_MTA_time);
        return;
    end
    
    if  length(GUIpp.PPsetting_MTA_timeVal)~=2
        statusbarObj.setText('Wrong input format of the MTA time interval');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_MTA_time);
        return;
    end
    
    if GUIpp.PPsetting_MTA_timeVal(1)>GUIpp.PPsetting_MTA_timeVal(2)
        statusbarObj.setText('Error:Wrong input format! [t_start;t_end]');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_MTA_time);
        return;
    end
    
    if GUIpp.PPsetting_MTA_timeVal(1)<tref(1)
        statusbarObj.setText(['Error:Input value t in MTA < minimum value (t=',num2str(tref(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_MTA_time);
        return;
    end
    
    if GUIpp.PPsetting_MTA_timeVal>tref(2)
        statusbarObj.setText(['Error:Input value t in MTA > maximum value (t=',num2str(tref(2)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(handles.edit_MTA_time);
        return;
    end
    
end
statusbarObj.setText('');

%uiresume(gcbf);
%save_state(GUIpp);
%run([GUIpp.path,'/Program/PostProc/','PostProc_HaWaSSI.m']);
PostProc_source();


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
But1=get(hObject,'Value');
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);

set_panel1(But1,handles);
set_panel11(But1,handles);




% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
But3=get(hObject,'Value');
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton4,'Value',0);


set_panel3(But3,handles);
set_panel11(But3,handles);



% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
But4=get(hObject,'Value');

set(handles.radiobutton1,'Value',0);
set(handles.radiobutton3,'Value',0);


set_panel2(But4,handles);
set_panel11(But4,handles);



function set_panel11(Id,handles)
if Id==1
    set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'on')
    if get(handles.checkbox_xlim,'value')==1, set(handles.edit28,'enable','on');
    else set(handles.edit28,'enable','off'); end;
    if get(handles.checkbox28,'value')==1, set(handles.edit29,'enable','on');
    else set(handles.edit29,'enable','off'); end;
    if get(handles.checkbox29,'value')==1, set(handles.edit30,'enable','on');
    else set(handles.edit30,'enable','off'); end;
    if get(handles.checkbox30,'value')==1, set(handles.edit31,'enable','on');
    else set(handles.edit31,'enable','off'); end;
    if get(handles.checkbox_MTA_time_interval,'value')==1, set(handles.edit_MTA_time,'enable','on');
    else set(handles.edit_MTA_time,'enable','off'); end;
    if get(handles.checkbox35,'value')==1, set(handles.popupmenu5,'enable','on');
    else set(handles.popupmenu5,'enable','off');end;
    if get(handles.radiobutton4,'value')==1
        set(handles.checkbox_saveanimation,'enable','on')
        set(handles.checkbox35,'enable','off');
        set(handles.popupmenu5,'enable','off');
    else
        set(handles.checkbox_saveanimation,'enable','off')
    end
    
else
    set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'off')
end


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
param=str2num(get(hObject,'String'));
set(handles.edit10,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
param=str2num(get(hObject,'String'));
set(handles.edit11,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
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
Input_gui = guidata(hObject);

[statusbarObj]=JavaFrame_handling();

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathnow=CurDir{get(handles.popupmenu_CurDir,'Value')};
[file_name,directory]=uigetfile([pathnow,'/','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
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
        set(Input_gui.pushbutton3,'Userdata',my_data);
        
        myBuoy_X=my_data(1,2:end);
        set(handles.edit10,'Userdata',myBuoy_X);
        
        set(handles.edit10,'String','')
        
        NX=length(myBuoy_X);
        
        if NX~=0
            if NX~=1
                for i=1:NX-1
                    Buoy_X{i}=[num2str(myBuoy_X(i)),';',];
                end
            end
            Buoy_X{NX}=num2str(myBuoy_X(NX));
            set(handles.edit10,'String',[Buoy_X{1:NX}]);
        end;
        
        statusbarObj.setText([file_name,' has been loaded']);
    else
        statusbarObj.setText('wrong input file!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    end
else
    statusbarObj.setText('');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
param=str2num(get(hObject,'String'));
set(handles.edit7,'Userdata',param);


[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit7);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
param=str2num(get(hObject,'String'));
set(handles.edit8,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit8);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
param=str2num(get(hObject,'String'));
set(handles.edit12,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit12);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
Id=get(hObject,'Value');
edit_filter=findobj('Tag','edit12');
if Id==1
    set(edit_filter,'Visible','on');
    set(handles.text11,'Visible','on');set(handles.text33,'Visible','on');
else
    set(edit_filter,'Visible','off');
    set(handles.text11,'Visible','off');set(handles.text33,'Visible','off');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
param=str2num(get(hObject,'String'));
set(handles.edit16,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit16);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton_load_simul.
function pushbutton_load_simul_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_simul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp('uploading data ...')
% Input_gui = guidata(hObject);
[statusbarObj]=JavaFrame_handling();

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathnow=CurDir{get(handles.popupmenu_CurDir,'Value')};

[file_name,directory]=uigetfile([pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
if directory~=0
    %reset_handles(handles,1)
    %cd(directory);
    statusbarObj.setText('loading...');
    my_data=load([directory,file_name]);
    my_data.Load=1;
    %     my_data=get(handles.pushbutton_load_simul,'Userdata');
    if isempty(my_data)
        set(handles.popupmenu1,'Value',1)
        statusbarObj.setText('No data loaded');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    elseif ~isfield(my_data,'par')||~isfield(my_data,'influx')
        set(handles.popupmenu1,'Value',1)
        statusbarObj.setText('Wrong input file.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
    else
        reset_handles(handles,1);
        [statusbarObj]=JavaFrame_handling();
        statusbarObj.setText('');
        set(handles.popupmenu1,'Value',2)
        set(handles.pushbutton_load_simul,'Enable','off')
        
        
        if isfield(my_data,'Proj')  %just indicates from the new version of saved data
            time=my_data.output.time;
            x=my_data.output.x;
            my_data.savename=my_data.Proj.savename;
            set(handles.edit_projectname,'string',my_data.Proj.savename);
        else
            time=my_data.time;
            x=my_data.x;
        end
        
        Fbdy=my_data.par.bf0;
        IndtI=1;%my_data.IndtI;
        IndtF=length(time);%my_data.IndtF;
        
%         if length(Fbdy)==2
%             PP_X=[floor(x(1)+Fbdy(1));floor(x(end)-Fbdy(2))];
%         else
%             PP_X=[floor(x(1)+Fbdy);floor(x(end)-Fbdy)];
%         end
        PP_X=[ceil(x(1));floor(x(end))];
        set(handles.edit_Profile_Position,'Userdata',PP_X);
        set(handles.edit_Profile_Position,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
        set(handles.edit8,'Userdata',PP_X);
        set(handles.edit8,'String',[num2str(PP_X(1)),';',num2str(PP_X(2))]);
        
        PP_T=[ceil(time(IndtI));floor(time(IndtF))];
        set(handles.edit_Signal_time,'Userdata',PP_T);
        set(handles.edit_Signal_time,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
        set(handles.edit7,'Userdata',PP_T);
        set(handles.edit7,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
        set(handles.edit16,'Userdata',PP_T);
        set(handles.edit16,'String',[num2str(PP_T(1)),';',num2str(PP_T(2))]);
        
        set(handles.pushbutton_load_simul,'Userdata',my_data);
                 
        statusbarObj.setText([file_name,' has been loaded']);
    end
 global PPdata
PPdata=[];
PPdata.workdir='';   
    
end




function save_state(state)
state.MySim=[]; %to reduce size and increase speed
state.PP4MeasDat=[];
state.gui='PostProc_GUI';
if ~isdir(state.workdir)
    mkdir(state.workdir);
end
save([state.workdir,'abproj1d_postproc.mat'], 'state')


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
%close all;

function keep_guis
fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
for fh = fig_h
    uih = findobj( fh, 'Type', 'uicontrol' );
    if isempty( uih )
        delete( fh );
    end
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
% --- Executes on button press in pushbutton_load_simul.

contents = cellstr(get(hObject,'String'));
set(handles.popupmenu1,'value',1);
Cont=contents{1};

Data_Sim_load =findobj('Tag','pushbutton_load_simul');

if strcmp(Cont,'Other')
    set(Data_Sim_load,'Enable','on');
else
    set(Data_Sim_load,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox_GIF.
function checkbox_GIF_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_GIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_GIF

if get(handles.checkbox_GIF,'value')==1
    set(handles.edit_GIF_setting,'enable','on')
else
    set(handles.edit_GIF_setting,'enable','off')
end

function edit_GIF_setting_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GIF_setting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_GIF_setting as text
%        str2double(get(hObject,'String')) returns contents of edit_GIF_setting as a double
param=str2num(get(hObject,'String'));
set(handles.edit_GIF_setting,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_GIF_setting);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_GIF_setting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GIF_setting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Profile.
function checkbox_Profile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Profile
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_Profile_time,'enable','on');
    set(handles.edit_Profile_Position,'enable','on');
    set(handles.checkbox_MTA,'enable','on');
    set(handles.checkbox10,'enable','on');
else
    set(handles.edit_Profile_time,'enable','off');
    set(handles.edit_Profile_Position,'enable','off');
    set(handles.checkbox_MTA,'enable','off');
    set(handles.checkbox10,'enable','off');
end
if get(handles.checkbox10,'value')==1
    set(handles.edit_scale,'enable','on');
else
    set(handles.edit_scale,'enable','off');
end
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);

Idbuoy=get(handles.checkbox_signal,'value');
if Idbuoy==0
set(handles.edit_Signal_time,'enable','off');
else
set(handles.edit_Signal_time,'enable','on');    
end


% --- Executes on button press in checkbox_signal.
function checkbox_signal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_signal




if get(handles.checkbox_signal,'value')==1
    set(handles.edit_Signal_time,'enable','on');
    set(handles.edit_Signal_Position,'enable','on');
    set(handles.checkbox12,'enable','on');
    set(handles.checkbox13,'enable','on');
    set(handles.checkbox_PP1_buoy_quant,'enable','on')
    
    
    if get(handles.checkbox13,'value')==1
        set(handles.edit23,'enable','on');
        set(handles.checkbox_Energy_Sp,'enable','on');
        set(handles.checkbox_Ampl_Sp,'enable','on');
    else
        set(handles.edit23,'enable','off');
        set(handles.checkbox_Energy_Sp,'enable','off');
        set(handles.checkbox_Ampl_Sp,'enable','off');
    end
    
    set(handles.checkbox_Energy,'value',0);
    set(handles.checkbox_Energy_in_x,'value',0);
    set(handles.checkbox_breakevents,'value',0);
    set(handles.checkbox_PP1_Hs_Sk,'value',0);
    set(handles.checkbox_PP1_MTA_MWL,'value',0);
    set(handles.checkbox_kinematic,'value',0);
    set(handles.checkbox_PP1_spectra,'value',0);
    
    Idbuoy=get(handles.checkbox_Profile,'value');
    if Idbuoy==0
        set(handles.edit_Profile_Position,'enable','off');
    else
        set(handles.edit_Profile_Position,'enable','on');
    end
   set(handles.checkbox_buoy_exceedance,'enable','on');
    
else
    set(handles.edit_Signal_time,'enable','off');
    set(handles.edit_Signal_Position,'enable','off');
    set(handles.checkbox12,'enable','off');
    set(handles.checkbox13,'enable','off');
    set(handles.checkbox_PP1_buoy_quant,'enable','off')
    set(handles.checkbox_Energy_Sp,'enable','off');
    set(handles.checkbox_Ampl_Sp,'enable','off');
    set(handles.checkbox_buoy_exceedance,'enable','off');
end



function edit_Signal_Position_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Signal_Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Signal_Position as text
%        str2double(get(hObject,'String')) returns contents of edit_Signal_Position as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Signal_Position,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit_Signal_Position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Signal_Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Signal_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Signal_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Signal_time as text
%        str2double(get(hObject,'String')) returns contents of edit_Signal_time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Signal_time,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Signal_time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Signal_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Signal_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Profile_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Profile_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Profile_time as text
%        str2double(get(hObject,'String')) returns contents of edit_Profile_time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Profile_time,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit_Profile_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Profile_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Profile_Position_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Profile_Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Profile_Position as text
%        str2double(get(hObject,'String')) returns contents of edit_Profile_Position as a double
param=str2num(get(hObject,'String'));
set(handles.edit_Profile_Position,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_Profile_Position);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_Profile_Position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Profile_Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox_Val_Spectra.
function checkbox_Val_Spectra_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Val_Spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Val_Spectra
if get(handles.checkbox_Val_Spectra,'value')==1
    set(handles.edit11,'enable','on')
    set(handles.checkbox__P3_energySp,'enable','on')
    set(handles.checkbox_Amplitude_SP,'enable','on')
    set(handles.checkbox_Amplitude_SP,'value',1);
    set(handles.checkbox__P3_energySp,'value',0)  
else
    set(handles.edit11,'enable','off')
    set(handles.checkbox__P3_energySp,'enable','off')
    set(handles.checkbox_Amplitude_SP,'enable','off')
end

% --- Executes on button press in checkbox_Val_quant.
function checkbox_Val_quant_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Val_quant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Val_quant


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17
if get(handles.checkbox17,'value')==1
    set(handles.edit24,'enable','on')
else
    set(handles.edit24,'enable','off')
end


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double
param=str2num(get(hObject,'String'));
set(handles.edit24,'Userdata',param);



% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_MTA.
function checkbox_MTA_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MTA


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
if get(handles.checkbox10,'value')==1
    set(handles.edit_scale,'enable','on');
else
    set(handles.edit_scale,'enable','off');
end



function edit_scale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scale as text
%        str2double(get(hObject,'String')) returns contents of edit_scale as a double
param=str2num(get(hObject,'String'));
set(handles.edit_scale,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13
if get(handles.checkbox13,'value')==1
    set(handles.edit23,'enable','on');
    set(handles.checkbox_Energy_Sp,'enable','on');
    set(handles.checkbox_Ampl_Sp,'enable','on');
    set(handles.checkbox_Energy_Sp,'value',0);
    set(handles.checkbox_Ampl_Sp,'value',1);
else
    set(handles.edit23,'enable','off');
    set(handles.checkbox_Energy_Sp,'enable','off');
    set(handles.checkbox_Ampl_Sp,'enable','off');
end

% --- Executes on button press in checkbox_Energy_Sp.
function checkbox_Energy_Sp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Energy_Sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Energy_Sp
if get(handles.checkbox_Energy_Sp,'value')==1
   set(handles.checkbox_Ampl_Sp,'value',0) 
else
   set(handles.checkbox_Ampl_Sp,'value',1)  
end

% --- Executes on button press in checkbox_Ampl_Sp.
function checkbox_Ampl_Sp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Ampl_Sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Ampl_Sp
if get(handles.checkbox_Ampl_Sp,'value')==1
   set(handles.checkbox_Energy_Sp,'value',0) 
else
   set(handles.checkbox_Energy_Sp,'value',1)  
end


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double
param=str2num(get(hObject,'String'));
set(handles.edit23,'Userdata',param);


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
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
function OpenProject_Callback(hObject, eventdata, handles)
% hObject    handle to OpenProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(findall(0, 'type', 'figure'));
HAWASSI_AB_startpage;

% --------------------------------------------------------------------
function SaveProject_Callback(hObject, eventdata, handles)
% hObject    handle to SaveProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
pathNow=CurDir{get(handles.popupmenu_CurDir,'Value')};
[file,browsedir] = uiputfile([pathNow,'\','abproj1d_postproc.mat'],'Save project file');
if browsedir~=0
save_handles(handles,browsedir);
end

% --------------------------------------------------------------------
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reset_handles(handles,1);

% --------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
quit;

% --------------------------------------------------------------------
function Calculator_Callback(hObject, eventdata, handles)
% hObject    handle to Calculator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Calculator();

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Documentation_Callback(hObject, eventdata, handles)
% hObject    handle to Documentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('2015_08_02_Manual_AB.pdf');

function [statusbarObj]=JavaFrame_handling
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
jFrame    = get(handle(gcf),'JavaFrame');
try
    jRootPane = jFrame.fHG1Client.getWindow;
catch
    jRootPane = jFrame.fHG2Client.getWindow;    %>=2014b
end
statusbarObj = com.mathworks.mwswing.MJStatusBar;
jRootPane.setStatusBar(statusbarObj);



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


% --- Executes on button press in checkbox__P3_energySp.
function checkbox__P3_energySp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox__P3_energySp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox__P3_energySp
Id=get(hObject,'Value');
if Id==1
    set(handles.checkbox_Amplitude_SP,'value',0);
else
    set(handles.checkbox_Amplitude_SP,'value',1);
end

% --- Executes on button press in checkbox_Amplitude_SP.
function checkbox_Amplitude_SP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Amplitude_SP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Amplitude_SP
Id=get(hObject,'Value');
if Id==1
    set(handles.checkbox__P3_energySp,'value',0);
else
    set(handles.checkbox__P3_energySp,'value',1);
end

function set_panel1(But1,handles)
if But1==1
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'on')
end

if get(handles.checkbox_Profile,'value')==1
    set(handles.edit_Profile_time,'enable','on');
    set(handles.edit_Profile_Position,'enable','on');
    set(handles.checkbox_MTA,'enable','on');
    set(handles.checkbox10,'enable','on');
else
    set(handles.edit_Profile_time,'enable','off');
    set(handles.edit_Profile_Position,'enable','off');
    set(handles.checkbox_MTA,'enable','off');
    set(handles.checkbox10,'enable','off');
end
if get(handles.checkbox10,'value')==1
    set(handles.edit_scale,'enable','on');
else
    set(handles.edit_scale,'enable','off');
end

if get(handles.checkbox_signal,'value')==1
    set(handles.edit_Signal_time,'enable','on');
    set(handles.edit_Signal_Position,'enable','on');
    set(handles.checkbox12,'enable','on');
    set(handles.checkbox13,'enable','on');
    set(handles.checkbox_PP1_buoy_quant,'enable','on')
    set(handles.checkbox_buoy_exceedance,'enable','on');
else
    set(handles.edit_Signal_time,'enable','off');
    set(handles.edit_Signal_Position,'enable','off');
    set(handles.checkbox12,'enable','off');
    set(handles.checkbox13,'enable','off');
    set(handles.checkbox_PP1_buoy_quant,'enable','off')
    set(handles.checkbox_buoy_exceedance,'enable','off');
end

if get(handles.checkbox13,'value')==1
    set(handles.edit23,'enable','on');
    set(handles.checkbox_Energy_Sp,'enable','on');
    set(handles.checkbox_Ampl_Sp,'enable','on');
else
    set(handles.edit23,'enable','off');
    set(handles.checkbox_Energy_Sp,'enable','off');
    set(handles.checkbox_Ampl_Sp,'enable','off');
end

if get(handles.checkbox_signal,'value')==1 ||...
   get(handles.checkbox_Profile,'value')==1 
set(handles.checkbox_Energy,'value',0)
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_breakevents,'value',0)
set(handles.checkbox_PP1_MTA_MWL,'value',0)
set(handles.checkbox_PP1_Hs_Sk,'value',0)
end
IdKin=get(handles.checkbox_kinematic,'value');
IdEn=get(handles.checkbox_Energy,'value');
if IdEn==1 ||IdKin==1||...
    get(handles.checkbox_breakevents,'value')==1
if IdEn==1||IdKin==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');    
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');     
end





set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
    
end


if But1==0
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
end
set(handles.radiobutton1,'Enable','On');
set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton3,'Enable','On');
set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton4,'Enable','On');

function set_panel2(But4,handles)
if But4==1
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'on');
end

if get(handles.checkbox17,'value')==1
    set(handles.edit24,'enable','on')
else
    set(handles.edit24,'enable','off')
end


if get(handles.checkbox_GIF,'value')==1
    set(handles.edit_GIF_setting,'enable','on')
else
    set(handles.edit_GIF_setting,'enable','off')
end


if But4==0
    set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
end
set(handles.radiobutton4,'Enable','On');


set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton3,'Enable','On');
set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton1,'Enable','On');

function set_panel3(But3,handles)

if But3==1
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'on')
end

if get(handles.checkbox32,'value')==1
    set(handles.checkbox33,'enable','on')
    if get(handles.checkbox33,'value')==1
        set(handles.edit34,'enable','off')
    else
        set(handles.edit34,'enable','on')
    end
else
    set(handles.checkbox33,'enable','off')
    set(handles.edit34,'enable','off')
end


if get(handles.checkbox_Val_Spectra,'value')==1
    set(handles.edit11,'enable','on')
    set(handles.checkbox__P3_energySp,'enable','on')
    set(handles.checkbox_Amplitude_SP,'enable','on')
else
    set(handles.edit11,'enable','off')
    set(handles.checkbox__P3_energySp,'enable','off')
    set(handles.checkbox_Amplitude_SP,'enable','off')
end

if But3==0
    set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
end
set(handles.radiobutton3,'Enable','On');
set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton1,'Enable','On');
set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton4,'Enable','On');

function reset_handles(handles,Flag)
if Flag==1
    [statusbarObj]=JavaFrame_handling();
    statusbarObj.setText('');
    
    JavProgressBar;
    jProgressBar.setVisible(0);
end

if isfield(handles,'input')
    handles.input=[];
end

set(handles.edit_projectname,'string','');  
set(handles.pushbutton_load_simul,'Userdata',[]);
set(handles.pushbutton_load_simul,'Enable','on');
set(handles.popupmenu1,'Value',1);
set(handles.checkbox1,'Value',0);
set(handles.edit12,'Userdata',[])  ;
set(handles.edit12,'String','')  ;
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton4,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.checkbox_Profile,'Value',0);
set(handles.edit_Profile_time,'Userdata',[]);
set(handles.edit_Profile_time,'String','');
set(handles.edit_Profile_Position,'Userdata',[]);
set(handles.edit_Profile_Position,'String','');

set(handles.checkbox_MTA,'Value',0);
set(handles.checkbox10,'Value',0);
set(handles.edit_scale,'Userdata',[]);
set(handles.edit_scale,'String','');
set(handles.checkbox_signal,'Value',0);
set(handles.edit_Signal_Position,'Userdata',[]);
set(handles.edit_Signal_Position,'String','');

set(handles.edit_Signal_time,'Userdata',[]);
set(handles.edit_Signal_time,'String','');
set(handles.checkbox12,'Value',0);
set(handles.checkbox13,'Value',0);
set(handles.edit23,'Value',[]);
set(handles.edit23,'String','');
set(handles.checkbox_Energy_Sp,'Value',0);
set(handles.checkbox_Ampl_Sp,'Value',0);

%PP=2
set(handles.edit7,'Userdata',[]);
set(handles.edit7,'String','');
set(handles.edit8,'Userdata',[]);
set(handles.edit8,'String','');
set(handles.checkbox5,'Value',0);
set(handles.checkbox17,'Value',0);
set(handles.edit24,'Userdata',0);
set(handles.edit24,'String','');
set(handles.checkbox_GIF,'Value',0);
set(handles.edit_GIF_setting,'Userdata',[]);
set(handles.edit_GIF_setting,'String','');
set(handles.edit10,'Userdata','');
set(handles.edit10,'String','');
set(handles.edit16,'Userdata',[]);
set(handles.edit16,'String','');
set(handles.checkbox16,'Value',0);
set(handles.checkbox_Val_Spectra,'Value',0);
set(handles.edit11,'Userdata', []);
set(handles.edit11,'String', '')
set(handles.checkbox__P3_energySp,'Value',0);
set(handles.checkbox_Amplitude_SP,'Value',0);
set(handles.checkbox_Val_quant,'Value',0);
set(handles.pushbutton3,'Userdata',[]);
set(handles.edit_GIF_setting,'Userdata',[0.1;1]);
set(handles.edit_GIF_setting,'String',[num2str(1),'/',num2str(10),';',num2str(1)]);

set(handles.checkbox32,'value',0);
set(handles.checkbox33,'value',0);
set(handles.edit34,'Userdata',[]);
set(handles.edit34,'String','');
set(handles.checkbox_xlim,'value',0);
set(handles.edit28,'Userdata',[]);
set(handles.edit28,'String','');
set(handles.checkbox28,'value',0);
set(handles.edit29,'Userdata',[]);
set(handles.edit29,'String','');
set(handles.checkbox29,'value',0);
set(handles.edit30,'Userdata',0);
set(handles.edit30,'String','');
set(handles.checkbox30,'value',0);
set(handles.edit31,'Userdata',[]);
set(handles.edit31,'String','');
set(handles.checkbox_MTA_time_interval,'value',0);
set(handles.edit_MTA_time,'Userdata',[]);
set(handles.edit_MTA_time,'String','');


set(handles.radiobutton1,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);
set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton1,'Enable','On');
set(findall(handles.uipanel5, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton3,'Enable','On');
set(findall(handles.uipanel4, '-property', 'enable'), 'enable', 'off')
set(handles.radiobutton4,'Enable','On');
set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'off')

set(handles.edit_scale,'Userdata',1);
set(handles.edit_scale,'String',[num2str(1),'/',num2str(1)]);
set(handles.edit24,'Userdata',1);
set(handles.edit24,'String',[num2str(1),'/',num2str(1)]);
set(handles.edit23,'Userdata',1);
set(handles.edit23,'String',num2str(1));
set(handles.edit11,'Userdata',1);
set(handles.edit11,'String',num2str(1));
set(handles.edit_GIF_setting,'Userdata',[0.01;inf]);
set(handles.edit_GIF_setting,'String','1/100;inf');




if Flag==1
    statusbarObj.setText('Clear.');
end

function save_handles(handles,path)
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');

GUIpp.path=path;
SimData=cellstr(get(handles.popupmenu1,'String'));
GUIpp.SimData=SimData{get(handles.popupmenu1,'Value')};
GUIpp.MySim=get(handles.pushbutton_load_simul,'Userdata');
GUIpp.savename=get(handles.edit_projectname,'string');
if isempty(handles.input)
    if isfield(GUIpp.MySim,'savename')
        GUIpp.simul_file=[GUIpp.MySim.savename,'_simul.mat'];
    else
        GUIpp.simul_file=[GUIpp.savename,'_simul.mat'];
    end
else
    savename=handles.input{1}.savename;
    GUIpp.simul_file=[savename,'_simul.mat'];
end

    


GUIpp.Filter=get(handles.checkbox1,'Value');
GUIpp.LFHF=get(handles.edit12,'Userdata');


GUIpp.PP1=get(handles.radiobutton1,'Value');
GUIpp.PP1_Prof=get(handles.checkbox_Profile,'Value');
GUIpp.PP1_Prof_T=get(handles.edit_Profile_time,'Userdata');
GUIpp.PP1_Prof_X=get(handles.edit_Profile_Position,'Userdata');
GUIpp.PP1_Prof_MTA=get(handles.checkbox_MTA,'Value');
GUIpp.PP1_Prof_Bathy=get(handles.checkbox10,'Value');
GUIpp.PP1_Prof_Bathy_scale=get(handles.edit_scale,'Userdata');

GUIpp.PP1_Buoy=get(handles.checkbox_signal,'Value');
GUIpp.PP1_Buoy_T=get(handles.edit_Signal_time,'Userdata');
GUIpp.PP1_Buoy_X=get(handles.edit_Signal_Position,'Userdata');
GUIpp.PP1_Buoy_Signal=get(handles.checkbox12,'Value');
GUIpp.PP1_Buoy_Spectrum=get(handles.checkbox13,'Value');
GUIpp.PP1_Buoy_Spectrum_Smooth=get(handles.edit23,'Userdata');
GUIpp.PP1_Buoy_Spectrum_Energy=get(handles.checkbox_Energy_Sp,'Value');
GUIpp.PP1_Buoy_Spectrum_Ampl=get(handles.checkbox_Ampl_Sp,'Value');

GUIpp.PP2=get(handles.radiobutton4,'Value');
GUIpp.PP2_T=get(handles.edit7,'Userdata');
GUIpp.PP2_X=get(handles.edit8,'Userdata');
GUIpp.PP2_MTA=get(handles.checkbox5,'Value');
GUIpp.PP2_Bathy=get(handles.checkbox17,'Value');
GUIpp.PP2_Bathy_scale=get(handles.edit24,'Userdata');
GUIpp.PP2_GIF=get(handles.checkbox_GIF,'Value');
GUIpp.PP2_GIF_setting=get(handles.edit_GIF_setting,'Userdata');

GUIpp.PP3=get(handles.radiobutton3,'Value');
GUIpp.PP3_X=get(handles.edit10,'Userdata');
GUIpp.PP3_T=get(handles.edit16,'Userdata');
GUIpp.PP3_Signals=get(handles.checkbox16,'Value');
GUIpp.PP3_Spectra=get(handles.checkbox_Val_Spectra,'Value');
GUIpp.PP3_Spectra_smooth=get(handles.edit11,'Userdata');
GUIpp.PP3_Spectra_Energy=get(handles.checkbox__P3_energySp,'Value');
GUIpp.PP3_Spectra_Ampl=get(handles.checkbox_Amplitude_SP,'Value');
GUIpp.PP3_Quant=get(handles.checkbox_Val_quant,'Value');
GUIpp.PP3MeasDat=get(handles.pushbutton3,'Userdata');

GUIpp.PPsetting_xlim=get(handles.checkbox_xlim,'value');
GUIpp.PPsetting_xlimVal=get(handles.edit28,'Userdata');
GUIpp.PPsetting_ylim=get(handles.checkbox28,'value');
GUIpp.PPsetting_ylimVal=get(handles.edit29,'Userdata');
GUIpp.PPsetting_xcoarse=get(handles.checkbox29,'value');
GUIpp.PPsetting_xcoarseVal=get(handles.edit30,'Userdata');
GUIpp.PPsetting_tcoarse=get(handles.checkbox30,'value');
GUIpp.PPsetting_tcoarseVal=get(handles.edit31,'Userdata');
GUIpp.PPsetting_MTA_time=get(handles.checkbox_MTA_time_interval,'value');
GUIpp.PPsetting_MTA_timeVal=get(handles.edit_MTA_time,'Userdata');
GUIpp.PPsetting_save_data=get(handles.checkbox34,'value');
GUIpp.save_figure=get(handles.checkbox35,'value');
filetype=cellstr(get(handles.popupmenu5, 'string'));
GUIpp.save_figure_type=filetype{get(handles.popupmenu5,'value')};

GUIpp.workdir=path;

save_state(GUIpp);
statusbarObj.setText('Project saved as GUIpp_state.mat');


% --- Executes on button press in checkbox_skew_assym.
function checkbox_skew_assym_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_skew_assym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_skew_assym


% --- Executes on button press in checkbox_MTA_MWL.
function checkbox_MTA_MWL_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MTA_MWL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MTA_MWL


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28
Id=get(hObject,'Value');
if Id==1
    set(handles.edit29,'enable','on')
else
    set(handles.edit29,'enable','off')
end


function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double
param=str2num(get(hObject,'String'));
set(handles.edit29,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit29);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_xlim.
function checkbox_xlim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_xlim
Id=get(hObject,'Value');
if Id==1
    set(handles.edit28,'enable','on')
else
    set(handles.edit28,'enable','off')
end


function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
param=str2num(get(hObject,'String'));
set(handles.edit28,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit28);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29
Id=get(hObject,'Value');
if Id==1
    set(handles.edit30,'enable','on')
else
    set(handles.edit30,'enable','off')
end


function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double
param=str2num(get(hObject,'String'));
set(handles.edit30,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit30);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30
Id=get(hObject,'Value');
if Id==1
    set(handles.edit31,'enable','on')
else
    set(handles.edit31,'enable','off')
end


function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double
param=str2num(get(hObject,'String'));
set(handles.edit31,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=1
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit31);
else %if correct
    statusbarObj.setText('');
end


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_MTA_time_interval.
function checkbox_MTA_time_interval_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MTA_time_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MTA_time_interval
Id=get(hObject,'Value');
if Id==1
    set(handles.edit_MTA_time,'enable','on')
else
    set(handles.edit_MTA_time,'enable','off')
end


function edit_MTA_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MTA_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MTA_time as text
%        str2double(get(hObject,'String')) returns contents of edit_MTA_time as a double
param=str2num(get(hObject,'String'));
set(handles.edit_MTA_time,'Userdata',param);

[statusbarObj]=JavaFrame_handling();
if length(param)~=2
    statusbarObj.setText('wrong input format');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(handles.edit_MTA_time);
else %if correct
    statusbarObj.setText('');
end

% --- Executes during object creation, after setting all properties.
function edit_MTA_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MTA_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32
if get(handles.checkbox32,'value')==1
    set(handles.checkbox33,'enable','on')
    set(handles.checkbox33,'value',1)
    set(handles.edit34,'enable','off')
else
    set(handles.checkbox33,'enable','off')
    set(handles.edit34,'enable','off')
end

% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox33
if get(handles.checkbox33,'value')==1
    set(handles.edit34,'enable','off')
else
    set(handles.edit34,'enable','on')
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
if length(param)~=1
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


% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox34


% --- Executes on button press in checkbox35.
function checkbox35_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox35
if get(handles.checkbox35,'value')==1, set(handles.popupmenu5,'enable','on');
else set(handles.popupmenu5,'enable','off');end;

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keep_guis;


% --- Executes on button press in checkbox_saveanimation.
function checkbox_saveanimation_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveanimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_saveanimation



function edit_projectname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_projectname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_projectname as text
%        str2double(get(hObject,'String')) returns contents of edit_projectname as a double


% --- Executes during object creation, after setting all properties.
function edit_projectname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_projectname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Energy.
function checkbox_Energy_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Energy
Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);


% --- Executes on button press in checkbox_breakevents.
function checkbox_breakevents_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_breakevents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_breakevents
set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);

set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');  
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);


% --- Executes on button press in checkbox_PP1_buoy_quant.
function checkbox_PP1_buoy_quant_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_buoy_quant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_buoy_quant
if get(hObject,'value')==1
   set(handles.edit23,'enable','on'); 
else
   set(handles.edit23,'enable','off');
end


% --- Executes on button press in checkbox_PP1_Hs_Sk.
function checkbox_PP1_Hs_Sk_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_Hs_Sk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_Hs_Sk
Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);

% --- Executes on button press in checkbox_PP1_MTA_MWL.
function checkbox_PP1_MTA_MWL_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_MTA_MWL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_MTA_MWL

Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);

% --- Executes on button press in checkbox_kinematic.
function checkbox_kinematic_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_kinematic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_kinematic
Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);


% --- Executes on button press in checkbox_Ship_XZThe.
function checkbox_Ship_XZThe_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Ship_XZThe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Ship_XZThe
if get(hObject,'Value')==1
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_Ship_forces,'value',0)
set(handles.checkbox_PP1_spectra,'value',0);

% --- Executes on button press in checkbox_Ship_forces.
function checkbox_Ship_forces_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Ship_forces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Ship_forces
if get(hObject,'Value')==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_Ship_XZThe,'value',0)
set(handles.checkbox_PP1_spectra,'value',0);


% --- Executes on button press in checkbox_buoy_exceedance.
function checkbox_buoy_exceedance_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_buoy_exceedance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_buoy_exceedance


% --- Executes on button press in checkbox_PP1_spectra.
function checkbox_PP1_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PP1_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PP1_spectra
Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');   
end    

set(handles.checkbox_Energy,'value',0);
set(handles.checkbox_Energy_in_x,'value',0);
set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);


% --- Executes on button press in checkbox_Energy_in_x.
function checkbox_Energy_in_x_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Energy_in_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Energy_in_x
Id=get(hObject,'Value');
if Id==1
set(handles.edit_Profile_Position,'Enable','on');
set(handles.edit_Signal_time,'Enable','on');
else
set(handles.edit_Profile_Position,'Enable','off');
set(handles.edit_Signal_time,'Enable','off');    
end    

set(handles.checkbox_Profile,'value',0);
set(handles.checkbox_signal,'value',0);
set(handles.edit_Profile_time,'enable','off');
set(handles.checkbox_MTA,'enable','off');
set(handles.checkbox10,'enable','off');
set(handles.edit_scale,'enable','off');

set(handles.edit_Signal_Position,'Enable','off');
set(handles.checkbox12,'enable','off');
set(handles.checkbox13,'enable','off');
set(handles.edit23,'enable','off');
set(handles.checkbox_Energy_Sp,'enable','off');
set(handles.checkbox_Ampl_Sp,'enable','off');
set(handles.checkbox_PP1_buoy_quant,'enable','off');

set(handles.checkbox_breakevents,'value',0);
set(handles.checkbox_PP1_Hs_Sk,'value',0);
set(handles.checkbox_PP1_MTA_MWL,'value',0);
set(handles.checkbox_kinematic,'value',0);
set(handles.checkbox_PP1_spectra,'value',0);
set(handles.checkbox_Energy,'value',0);