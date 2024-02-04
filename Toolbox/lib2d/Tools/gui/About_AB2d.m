function varargout = About_AB2d(varargin)
% ABOUT_AB2D MATLAB code for About_AB2d.fig
%      ABOUT_AB2D, by itself, creates a new ABOUT_AB2D or raises the existing
%      singleton*.
%
%      H = ABOUT_AB2D returns the handle to a new ABOUT_AB2D or the handle to
%      the existing singleton*.
%
%      ABOUT_AB2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUT_AB2D.M with the given input arguments.
%
%      ABOUT_AB2D('Property','Value',...) creates a new ABOUT_AB2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before About_AB2d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to About_AB2d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help About_AB2d

% Last Modified by GUIDE v2.5 15-Nov-2017 09:12:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @About_AB2d_OpeningFcn, ...
                   'gui_OutputFcn',  @About_AB2d_OutputFcn, ...
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


% --- Executes just before About_AB2d is made visible.
function About_AB2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to About_AB2d (see VARARGIN)

% Choose default command line output for About_AB2d
global IDdemoAbout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes2);     % Show the LabMath.

%[pathstr,~,~] = fileparts(mfilename('fullpath'));
%javaFrame = get(hObject,'JavaFrame');
%if ~isdeployed
%     pathstr=pwd;
    imshow(['hawassi_black.jpg']); %
%   %  javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
% else
%     imshow('hawassi_black.jpg');
%  %   javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'))
% end  
IDdemoAbout=0;
if IDdemoAbout==0
set(handles.figure1,'Name',['HAWASSI AB2: About']);
set(handles.text3,'string','AB 2D v.1.1 (01-01-2018)')
else
set(handles.figure1,'Name',['HAWASSI AB2 demo: About']);
set(handles.text3,'string','AB 2D v.1.1 demo (01-01-2018)')
end
log_string{1}='HAWASSI stands for';
log_string{2}='Hamiltonian Wave-Ship-Structure Interaction';
log_string{4}='HAWASSI–AB is a spatial-spectral implementation' ;
log_string{5}='           of the Analytic Boussinesq Model (AB)';
log_string{6}='Presently the code is for simulation of wave-structure interactions,';
log_string{7}='coupled wave-ship interaction is foreseen in future releases.';



set(handles.listbox2,'string',{log_string{:}})
% UIWAIT makes About_AB2d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = About_AB2d_OutputFcn(hObject, eventdata, handles) 
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
