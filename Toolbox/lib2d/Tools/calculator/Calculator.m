function varargout = Calculator(varargin)
% CALCULATOR MATLAB code for Calculator.fig
%      CALCULATOR, by itself, creates a new CALCULATOR or raises the existing
%      singleton*.
%
%      H = CALCULATOR returns the handle to a new CALCULATOR or the handle to
%      the existing singleton*.
%
%      CALCULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCULATOR.M with the given input arguments.
%
%      CALCULATOR('Property','Value',...) creates a new CALCULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Calculator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Calculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Calculator

% Last Modified by GUIDE v2.5 22-Jul-2015 22:14:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Calculator_OpeningFcn, ...
                   'gui_OutputFcn',  @Calculator_OutputFcn, ...
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


% --- Executes just before Calculator is made visible.
function Calculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Calculator (see VARARGIN)

% Choose default command line output for Calculator
global IDdemoCal
IDdemoCal=0;

handles.output = hObject;
handles.input  = varargin;
% Update handles structure
guidata(hObject, handles);

set(findobj(handles.figure1, '-property', 'FontName'),'FontName', 'default')
set(findobj(handles.figure1, '-property', 'Units'),'Units', 'normalized')
set(handles.edit2,'value',0);  
set(handles.edit2,'enable','off'); 
set(handles.edit3,'value',0);
set(handles.edit3,'enable','off');    

[pathstr,~,~] = fileparts(mfilename('fullpath'));

%javaFrame = get(hObject,'JavaFrame');
if ~isdeployed
cd (pathstr); cd ..;
pathstr=pwd;
addpath(genpath(pathstr))
%javaFrame.setFigureIcon(javax.swing.ImageIcon([pathstr,'\Misc\logo\icon_16.png']));
else
%javaFrame.setFigureIcon(javax.swing.ImageIcon('icon_16.png'));    
end

Disp_O       =findobj('Tag','togglebutton_disp');
Group_vel_O  =findobj('Tag','togglebutton_groupvel');
set(Disp_O,'Visible','off');  set(Group_vel_O,'Visible','off'); 
DispEq='sign(k).*sqrt(9.81*k.*tanh(d.*k))';
set(Disp_O,'userdata',DispEq);
UgEq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)';
set(Group_vel_O,'userdata',UgEq);


set(handles.figure1,'Name',['HAWASSI AB: Calculator']);


% UIWAIT makes Calculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Calculator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[statusbarObj]=JavaFrame_handling();

Param1_check=get(handles.popupmenu1,'value');
Param1_value=get(handles.edit1,'userdata');
Param2_check=get(handles.checkbox1,'value');
Depth       =get(handles.edit2,'userdata');
Param3_check=get(handles.checkbox2,'value');
Amplitude     =get(handles.edit3,'userdata');
disptype=cellstr(get(handles.popupmenu2,'String'));
disp=disptype{get(handles.popupmenu2,'Value')};
Id_disp=get(handles.popupmenu2,'Value');
if strcmp(disp,'Exact')
    fundisp='funCal_OmExact';
    funUg='funCal_UgExact';
    omAdd='';ugAdd='';
elseif strcmp(disp,'SWE')
    fundisp='funCal_OmSWE';
    funUg='funCal_UgSWE';
    omAdd='';ugAdd='';  
elseif strcmp(disp,'KdV')
    fundisp='funCal_OmKdV';
    funUg='funCal_UgKdV';
    omAdd='';ugAdd='';
elseif strcmp(disp,'BBM')
    fundisp='funCal_OmBBM';
    funUg='funCal_UgBBM';
    omAdd='';ugAdd='';
else
    fundisp='funCal_OmUser';
    funUg='funCal_UgUser';
    omAdd=get(handles.togglebutton_disp,'Userdata');
    ugAdd=get(handles.togglebutton_groupvel,'Userdata');

    if isempty(omAdd) 
     statusbarObj.setText('Input a dispersion formula');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
    if isempty(strfind(omAdd,'k')) 
     statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
    if isempty(strfind(omAdd,'d')) 
     statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
    if isempty(ugAdd)
     statusbarObj.setText('Input a group velocity formula');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
    if isempty(strfind(ugAdd,'k')) 
     statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
    if isempty(strfind(ugAdd,'d')) 
     statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
     statusbarTxt = statusbarObj.getComponent(0);
     statusbarTxt.setForeground(java.awt.Color.red);   
     return;
    end
    
end 

Om            = str2func(fundisp);
Ug            = str2func(funUg);


if isempty(Param1_value)
       statusbarObj.setText('Input a value');
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red); 
     return;
else
    statusbarObj.setText('');
end 

if Param2_check==1
    if isempty(Depth)
       statusbarObj.setText('Input a value');
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red); 
     return
    
    else
    statusbarObj.setText('');
    end
end 

if Param3_check==1
    if isempty(Amplitude )
       statusbarObj.setText('Input a value');
       statusbarTxt = statusbarObj.getComponent(0);
       statusbarTxt.setForeground(java.awt.Color.red); 
     return
    
    else
    statusbarObj.setText('');
    end
end 
   

if Param2_check==0 && Param3_check==0
    if Param1_check==1
        Tp=Param1_value;
        Freq_s=1./Tp;
        Freq_rad=2*pi*Freq_s;
        
        string{1}=['Frequency=',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                 =',num2str(Freq_rad,3),' [rad/s]'];
        
    elseif Param1_check==2
        Freq_s=Param1_value;
        Freq_rad=2*pi*Freq_s;
        Tp=1/Freq_s;
         string{1}=['Frequency=',num2str(Freq_rad,3),' [rad/s]'];
         string{2}=['Period       =',num2str(Tp,3),' [s]'];
    elseif Param1_check==3
        Freq_rad=Param1_value;
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
         string{1}=['Frequency=',num2str(Freq_s,3),' [1/s]'];
         string{2}=['Period       =',num2str(Tp,3),' [s]'];
    else
        lambda=Param1_value;
        k=2*pi/lambda;
        string{1}=['Wave-number=',num2str(k,3)];
    end
elseif Param2_check==1 && Param3_check==0
     if Param1_check==1
        Tp=Param1_value;
        Freq_s=1./Tp;
        Freq_rad=2*pi*Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        kp_h    =kp*Depth;
        
        lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(kp)|| kp<0
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
    
        string{4}=['------------------------------------------------------------'];
        string{5}=['There is no more info that can be generated.'];
        string{6}=['The dispersion is not reliable'];
        string{7}=['for the input wave frequency.'];
        string{8}=['Consider to increase period'];
        string{9}=['or decrease depth or change the dispersion.'];
        string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{9}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{11}=['Dispersion                                   = ',disp];
        end
    elseif Param1_check==2
        Freq_s=Param1_value;
        Freq_rad=2*pi*Freq_s;
        Tp=1/Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        kp_h    =kp*Depth;
          lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(kp)|| kp<0
        string{1}=['Frequency                                   = ',num2str(Freq_rad,3),' [rad/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];    
        string{4}=['------------------------------------------------------------'];
        string{5}=['There is no more info that can be generated.'];
        string{6}=['The dispersion is not reliable'];
        string{7}=['for the input wave frequency.'];
        string{8}=['Consider to decrease frequency or '];
        string{9}=['decrease depth or change the dispersion.'];
        string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_rad,3),' [rad/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{9}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{11}=['Dispersion                                   = ',disp];
        end
        
    elseif Param1_check==3
        Freq_rad=Param1_value;
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
        
          lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(kp)|| kp<0
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];    
        string{4}=['------------------------------------------------------------'];
        string{5}=['There is no more info that can be generated.'];
        string{6}=['The dispersion is not reliable'];
        string{7}=['for the input wave frequency.'];
        string{8}=['Consider to decrease frequency  or'];
        string{9}=['decrease depth or change the dispersion.'];
        string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{9}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{11}=['Dispersion                                   = ',disp];
        end
        
        
    else
        lambda=Param1_value;
        kp=2*pi/lambda;
        Freq_rad=Om(kp,Depth,omAdd);
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
        
          lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(Freq_rad)|| Freq_rad<0
        string{1}=['Wave number (k)                        = ',num2str(kp,3)];
        string{2}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{3}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{4}=['                     (',category,')'];
        string{6}=['------------------------------------------------------------'];
        string{7}=['There is no more info that can be generated.'];
        string{8}=['The dispersion is not reliable for kh>1.'];
        string{9}=['Consider to increase wave-length  or'];
        string{10}=['decrease depth or change the dispersion.'];
        string{11}=['------------------------------------------------------------'];
        else
        
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
        string{3}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{4}=['Wave number (k)                        = ',num2str(kp,3)];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{9}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{11}=['Dispersion                                   = ',disp];
        end
    end

elseif Param2_check==0 && Param3_check==1
    if Param1_check==1
        Tp=Param1_value;
        Freq_s=1./Tp;
        Freq_rad=2*pi*Freq_s;
        
        string{1}=['Frequency=',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                 =',num2str(Freq_rad,3),' [rad/s]'];
        
    elseif Param1_check==2
        Freq_s=Param1_value;
        Freq_rad=2*pi*Freq_s;
        Tp=1/Freq_s;
         string{1}=['Frequency=',num2str(Freq_rad,3),' [rad/s]'];
         string{2}=['Period       =',num2str(Tp,3),' [s]'];
    elseif Param1_check==3
        Freq_rad=Param1_value;
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
         string{1}=['Frequency=',num2str(Freq_s,3),' [1/s]'];
         string{2}=['Period       =',num2str(Tp,3),' [s]'];
    else
        lambda=Param1_value;
        k=2*pi/lambda;
        ka=k*Amplitude;
        string{1}=['Wave-number=',num2str(k,3)];
        string{2}=['Steepness     =',num2str(ka,3)];
    end
    
else
    if Param1_check==1
        Tp=Param1_value;
        Freq_s=1./Tp;
        Freq_rad=2*pi*Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
        lambda_per_h=lambda/Depth;
        kp_a     =kp*Amplitude;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
  
        if isnan(kp)|| kp<0
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
        string{4}=['------------------------------------------------------------'];
        string{5}=['There is no more info that can be generated.'];
        string{6}=['The dispersion is not reliable'];
        string{7}=['for the input wave frequency.'];
        string{8}=['Consider to increase period'];
        string{9}=['or decrease depth or change the dispersion.'];
        string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Steepness (k*a)                          = ',num2str(kp_a,3)];
        string{9}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{10}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{12}=['Dispersion                                   = ',disp];
        end
        

        
        
    elseif Param1_check==2
        Freq_s=Param1_value;
        Freq_rad=2*pi*Freq_s;
        Tp=1/Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
          lambda_per_h=lambda/Depth;
           kp_a     =kp*Amplitude;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(kp)|| kp<0
            string{1}=['Frequency                                   = ',num2str(Freq_rad,3),' [rad/s]'];
            string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];
            string{4}=['------------------------------------------------------------'];
            string{5}=['There is no more info that can be generated.'];
            string{6}=['The dispersion is not reliable'];
            string{7}=['for the input wave frequency.'];
            string{8}=['Consider to decrease frequency or '];
            string{9}=['decrease depth or change the dispersion.'];
            string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_rad,3),' [rad/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Steepness (k*a)                          = ',num2str(kp_a,3)];
        string{9}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{10}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{12}=['Dispersion                                   = ',disp];
        end
        
    elseif Param1_check==3
        Freq_rad=Param1_value;
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
        kp      =funCal_invOm(Freq_rad,Depth,Om,omAdd);
        lambda  =2*pi./kp;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
         kp_a     =kp*Amplitude;
        
          lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
         
        if isnan(kp)|| kp<0
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];    
        string{4}=['------------------------------------------------------------'];
        string{5}=['There is no more info that can be generated.'];
        string{6}=['The dispersion is not reliable'];
        string{7}=['for the input wave frequency.'];
        string{8}=['Consider to decrease frequency  or'];
        string{9}=['decrease depth or change the dispersion.'];
        string{10}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{3}=['Wave number (k)                        = ',num2str(kp,3)];
        string{4}=['Wave length  (lambda)                = ',num2str(lambda,3), '[m]'];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Steepness (k*a)                          = ',num2str(kp_a,3)];
        string{9}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{10}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{12}=['Dispersion                                   = ',disp];
        end
        
        
    
    else
        lambda=Param1_value;
        kp=2*pi/lambda;
        Freq_rad=Om(kp,Depth,omAdd);
        Freq_s=Freq_rad./(2*pi);
        Tp=1/Freq_s;
        Cp      =Freq_rad./kp;
        Vg      =Ug(kp,Depth,ugAdd);
        
        kp_h    =kp*Depth;
         kp_a     =kp*Amplitude;
        
          lambda_per_h=lambda/Depth;
        
        if lambda_per_h<2
            category='Deep water';
        elseif lambda_per_h>20
            category='Shallow water';
        else
            category='Intermediate depth';
        end
        
        if isnan(Freq_rad)|| Freq_rad<0
        string{1}=['Wave number (k)                        = ',num2str(kp,3)];
        string{2}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{3}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{4}=['                     (',category,')'];
        string{6}=['------------------------------------------------------------'];
        string{7}=['There is no more info that can be generated.'];
        string{8}=['The dispersion is not reliable forthe input wave-length.'];
        string{9}=['Consider to increase wave-length  or'];
        string{10}=['decrease depth or change the dispersion.'];
        string{11}=['------------------------------------------------------------'];
        else
        string{1}=['Frequency                                   = ',num2str(Freq_s,3),' [1/s]'];
        string{2}=['                                                    = ',num2str(Freq_rad,3),' [rad/s]'];
        string{3}=['Period                                          = ',num2str(Tp,3),' [s]'];
        string{4}=['Wave number (k)                        = ',num2str(kp,3)];
        string{5}=['Relative wave-length (lambda/h)= ',num2str(lambda_per_h,3)];
        string{6}=['                                   (k*h)         = ',num2str(kp_h,3)];
        string{7}=['                     (',category,')'];
        string{8}=['Steepness (k*a)                          = ',num2str(kp_a,3)];
        string{9}=['Phase speed                               = ',num2str(Cp,3), '[m/s]'];
        string{10}=['Group speed                               = ',num2str(Vg,3), '[m/s]'];  
        
        string{12}=['Dispersion                                   = ',disp];
        end
    end
       
    
end


set(handles.listbox1,'string',string);

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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
val=get(hObject,'Value');
if val==1
set(handles.edit2,'enable','on');    
else
set(handles.edit2,'enable','off');    
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
param=str2double(get(hObject,'String'));
set(handles.edit1,'Userdata',param);

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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
param=str2double(get(hObject,'String'));
set(handles.edit2,'Userdata',param);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
val=get(hObject,'Value');
if val==1
set(handles.edit3,'enable','on');    
else
set(handles.edit3,'enable','off');    
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
param=str2double(get(hObject,'String'));
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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
contents = cellstr(get(hObject,'String'));
Disp=contents{get(hObject,'Value')};

global IDdemoCal

if IDdemoCal==1
    if get(hObject,'Value')~=1
        set(handles.popupmenu2,'Value',1);
        [statusbarObj]=JavaFrame_handling();
        statusbarObj.setText('Only exact dispersion is available.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.blue);
    end
    
else

Disp_O       =findobj('Tag','togglebutton_disp');
Group_vel_O  =findobj('Tag','togglebutton_groupvel');

if strcmp(Disp,'User-defined')
   set(Disp_O,'Visible','on'); set(Group_vel_O,'Visible','on'); 
else
     set(Disp_O,'Visible','off');  set(Group_vel_O,'Visible','off'); 
end
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in togglebutton_disp.
function togglebutton_disp_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_disp
Id=get(hObject,'Value');

def.eq=get(handles.togglebutton_disp,'userdata');
if isempty(def.eq)
def.eq='sign(k).*sqrt(9.81*k.*tanh(d.*k))'; 
end
def.hf=handles.togglebutton_disp;
if Id==1
     input_formula_Disp(def);
else
 [hf,yy]=  input_formula_Disp;
 set(handles.togglebutton_disp,'userdata',yy);
 delete(input_formula_Disp)
end

% --- Executes on button press in togglebutton_groupvel.
function togglebutton_groupvel_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_groupvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_groupvel
Id=get(hObject,'Value');

def.eq=get(handles.togglebutton_groupvel,'userdata');
if isempty(def.eq)
def.eq='sign(k).*sqrt(9.81)/2./(k.*tanh(d.*k))^(1/2).*(tanh(d.*k)+k*(1-tanh(d*k)^2).*d)'; 
end

def.hf=handles.togglebutton_groupvel;
if Id==1
     input_formula_Vg(def);
else
 [hf,yy]=  input_formula_Vg;
 set(handles.togglebutton_groupvel,'userdata',yy);
 delete(input_formula_Vg)
end
