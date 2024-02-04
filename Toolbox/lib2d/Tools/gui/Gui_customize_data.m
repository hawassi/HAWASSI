function Gui_customize_data(dat,savedir)

hh.fig = figure('Name', strcat('Customize output data'),'unit','normalized', 'Position',...
    [0.5 0.4 0.2 0.3],'NumberTitle', 'Off','MenuBar', 'None', 'Toolbar', 'None',...
    'Color', [.94 .94 .94]);

panel = uipanel('Parent', hh.fig);
hh.spatial_domain_title  = uicontrol('Style','text','String','Customize output data', ...
    'Parent',panel, 'Unit', 'Normalized',...
    'Position', [.01 .93 .98 .06], 'Backgroundcolor', [0.8 0.8 0.8], ...
    'Fontunits','normalized','FontSize', 0.8, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'left');
hh.spatial_text0 = uicontrol('Style','text','String','Interval:',...
    'HorizontalAlignment', 'left', 'Parent', panel,...
    'Unit', 'Normalized', 'Position', [.022 .74 .1 .1],'Fontunits',...
    'normalized','FontSize', 0.3);
hh.spatial_text1 = uicontrol('Style','text','String','Xmin:',...
    'HorizontalAlignment', 'left', 'Parent', panel,...
    'Unit', 'Normalized', 'Position', [.05 .65 .1 .1],'Fontunits',...
    'normalized','FontSize', 0.3);
hh.spatial_edit_xmin = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.12 .705 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text2 = uicontrol('Style','text','String','[m]', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.25 .65 .1 .1],'Fontunits', 'normalized','FontSize', 0.3);
hh.spatial_text2a = uicontrol('Style','text','String','Xmax:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.3 .65 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_xmax = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel, 'Unit', 'Normalized',...
    'Position', [.38 .705 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text3 = uicontrol('Style','text','String','[m]', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.51 .65 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_text3a = uicontrol('Style','text','String','Coarse:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.58 .65 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_xcoarse = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel, 'Unit', 'Normalized',...
    'Position', [.68 .705 .12  .05],'Fontunits','normalized','FontSize', 0.7);

hh.spatial_text4 = uicontrol('Style','text','String','Ymin:', ...
    'HorizontalAlignment', 'left', 'Parent',panel ,...
    'Unit', 'Normalized', 'Position', [.05 .55 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_ymin = uicontrol('Style', 'edit', 'String','',...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w',...
    'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.12 .605 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text5 = uicontrol('Style','text','String','[m]', 'HorizontalAlignment',...
    'left', 'Parent',panel , 'Unit', 'Normalized',...
    'Position', [.25 .55 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_text5a = uicontrol('Style','text','String','Ymax:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.3 .55 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_ymax = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Parent',...
    panel , 'Unit', 'Normalized',...
    'Position', [.38 .605 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text6 = uicontrol('Style','text','String','[m]', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized', ...
    'Position', [.51 .55 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_text7 = uicontrol('Style','text','String','Coarse:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.58 .55 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_ycoarse = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel, 'Unit', 'Normalized',...
    'Position', [.68 .605 .12  .05],'Fontunits','normalized','FontSize', 0.7);

hh.spatial_text8 = uicontrol('Style','text','String','Tmin:', ...
    'HorizontalAlignment', 'left', 'Parent',panel ,...
    'Unit', 'Normalized', 'Position', [.05 .45 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_tmin = uicontrol('Style', 'edit', 'String','',...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w',...
    'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.12 .505 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text9 = uicontrol('Style','text','String','[s]', 'HorizontalAlignment',...
    'left', 'Parent',panel , 'Unit', 'Normalized',...
    'Position', [.25 .45 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_text10 = uicontrol('Style','text','String','Tmax:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.3 .45 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_tmax = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Parent',...
    panel , 'Unit', 'Normalized',...
    'Position', [.38 .505 .12  .05],'Fontunits','normalized','FontSize', 0.7);
hh.spatial_text11 = uicontrol('Style','text','String','[s]', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized', ...
    'Position', [.51 .45 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_text12 = uicontrol('Style','text','String','Coarse:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.58 .45 .1 .1],'Fontunits','normalized','FontSize', 0.3);
hh.spatial_edit_tcoarse = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel, 'Unit', 'Normalized',...
    'Position', [.68 .505 .12  .05],'Fontunits','normalized','FontSize', 0.7);

hh.options_outputvar_text = uicontrol('Style','text','String','Output variables:', ...
    'HorizontalAlignment', 'left', 'Parent', panel , ...
    'Unit', 'Normalized', 'Position', [.04 .32 .3 .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.3);
hh.options_outputvar_popup = uicontrol('Style','popupmenu','String',{'elevation','elevation & potential','all'}, ...
    'HorizontalAlignment', 'left', 'Parent', panel, ...
    'Unit', 'Normalized', 'Position', [.3 .33 .4 .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.3,'backgroundcolor','w');
hh.options_outputformat_text = uicontrol('Style','text','String','data format:', ...
    'HorizontalAlignment', 'left', 'Parent', panel , ...
    'Unit', 'Normalized', 'Position', [.04 .22 .3 .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.3);
hh.options_outputformat_popup = uicontrol('Style','popupmenu','String',{'.mat','.txt'}, ...
    'HorizontalAlignment', 'left', 'Parent', panel, ...
    'Unit', 'Normalized', 'Position', [.3 .23 .4 .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.3,'backgroundcolor','w');
hh.options_text_savename = uicontrol('Style','text','String','Save name:', 'HorizontalAlignment',...
    'left', 'Parent', panel , 'Unit', 'Normalized',...
    'Position', [.04 .12 .2 .1],'Fontunits','normalized','FontSize', 0.3);
hh.options_edit_savename = uicontrol('Style', 'edit', 'String','', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
    'Parent', panel, 'Unit', 'Normalized',...
    'Position', [.3 .17 .4  .06],'Fontunits','normalized','FontSize', 0.7);

hh.pp_project_button_save= uicontrol('Style','pushbutton','String','save data',...
    'HorizontalAlignment', 'left', 'Parent', panel,...
    'Unit', 'Normalized', 'Position', [.3 .08 .3 .07],'Fontunits',...
    'normalized','FontSize', 0.5,'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1=dat.output.X(1);X2=dat.output.X(end);
Y1=dat.output.Y(1);Y2=dat.output.Y(end);
T1=dat.output.time(1);T2=dat.output.time(end);
set(hh.spatial_edit_xmin,'userdata',X1,'string',num2str(X1));
set(hh.spatial_edit_xmax,'userdata',X2,'string',num2str(X2));
set(hh.spatial_edit_xcoarse,'userdata',1,'string','1');
set(hh.spatial_edit_ymin,'userdata',Y1,'string',num2str(Y1));
set(hh.spatial_edit_ymax,'userdata',Y2,'string',num2str(Y2));
set(hh.spatial_edit_ycoarse,'userdata',1,'string','1');
set(hh.spatial_edit_tmin,'userdata',T1,'string',num2str(T1));
set(hh.spatial_edit_tmax,'userdata',T2,'string',num2str(T2));
set(hh.spatial_edit_tcoarse,'userdata',1,'string','1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Callback functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(hh.spatial_edit_xmin,'callback',{@callback_edit_param})
set(hh.spatial_edit_xmax,'callback',{@callback_edit_param})
set(hh.spatial_edit_xcoarse,'callback',{@callback_edit_param})
set(hh.spatial_edit_ymin,'callback',{@callback_edit_param})
set(hh.spatial_edit_ymax,'callback',{@callback_edit_param})
set(hh.spatial_edit_ycoarse,'callback',{@callback_edit_param})
set(hh.spatial_edit_tmin,'callback',{@callback_edit_param})
set(hh.spatial_edit_tmax,'callback',{@callback_edit_param})
set(hh.spatial_edit_tcoarse,'callback',{@callback_edit_param})

set(hh.options_outputvar_popup,'callback',{@callback_popup_var,hh});
set(hh.options_outputformat_popup,'callback',{@callback_popup_format,hh});

set(hh.pp_project_button_save,'callback',{@callback_save_data,hh,dat,savedir});

function callback_popup_var(hObj,~,hh)
 Id=get(hObj,'value');
 if Id==3
 set(hh.options_outputformat_popup,'value',1); 
 end
end

function callback_popup_format(hObj,~,hh)
 Idvar=get(hh.options_outputvar_popup,'value');
 if Idvar==3
 set(hObj,'value',1); 
 end
end

function callback_edit_param(hObj,eventdata)
par=get(hObj,'string');
set(hObj,'userdata',str2num(par));
end

function callback_save_data(hObj,eventdat,hh,dat,savedir)
x1=get(hh.spatial_edit_xmin,'userdata');
x2=get(hh.spatial_edit_xmax,'userdata');
stepx=get(hh.spatial_edit_xcoarse,'userdata');
y1=get(hh.spatial_edit_ymin,'userdata');
y2=get(hh.spatial_edit_ymax,'userdata');
stepy=get(hh.spatial_edit_ycoarse,'userdata');
t1=get(hh.spatial_edit_tmin,'userdata');
t2=get(hh.spatial_edit_tmax,'userdata');
stept=get(hh.spatial_edit_tcoarse,'userdata');
Idvar=get(hh.options_outputvar_popup,'value');
Idformat=get(hh.options_outputformat_popup,'value');
savename=get(hh.options_edit_savename,'string');
[statusbarObj]=funGui_JavaFrame_handling();
statusbarObj.setText('');
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.black);

if isempty(x1)
    statusbarObj.setText('Specify xmin.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmin);
    return;
end
if isempty(x2)
    statusbarObj.setText('Specify xmax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end

if x1>=x2
    statusbarObj.setText('xmin must be smaller than xmax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;     
end
if isempty(stepx)
    statusbarObj.setText('Specify coarse fact. for x.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end

if isempty(y1)
    statusbarObj.setText('Specify ymin.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmin);
    return;
end
if isempty(y2)
    statusbarObj.setText('Specify ymax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end


if y1>=y2
    statusbarObj.setText('ymin must be smaller than ymax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;     
end

if isempty(stepy)
    statusbarObj.setText('Specify coarse fact. for y.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end

if isempty(t1)
    statusbarObj.setText('Specify tmin.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmin);
    return;
end
if isempty(t2)
    statusbarObj.setText('Specify tmax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end


if t1>=t2
    statusbarObj.setText('tmin must be smaller than tmax.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;     
end

if isempty(stepy)
    statusbarObj.setText('Specify coarse fact. for time.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.spatial_edit_xmax);
    return;
end

if Idvar==2
    if isfield(dat.output,'phi')
    if isempty(dat.output.phi)
      statusbarObj.setText('there is no potential data . . .');
      set(hh.options_outputvar_popup,'value',1);
      statusbarTxt = statusbarObj.getComponent(0);
      statusbarTxt.setForeground(java.awt.Color.red);
      return; 
    end   
    else
        statusbarObj.setText('there is no potential data . . .');
        set(hh.options_outputvar_popup,'value',1);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
end


if strcmpi(savename,'')
    statusbarObj.setText('Specify a savename');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(hh.options_edit_savename);
    return;
end

statusbarObj.setText('saving...');
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.black);

etadat=dat.output.eta;
xx=dat.output.X;
yy=dat.output.Y;
tt=dat.output.time;
indt1=funC_closest(tt,t1);
indt2=funC_closest(tt,t2);
indx1=funC_closest(xx,x1);
indx2=funC_closest(xx,x2);
indy1=funC_closest(yy,y1);
indy2=funC_closest(yy,y2);
output.X=xx(indx1:stepx:indx2);
output.Y=yy(indy1:stepy:indy2);
output.time=tt(indt1:stept:indt2);
output.eta=etadat(indt1:stept:indt2,indy1:stepy:indy2,indx1:stepx:indx2);
Nt=length(output.time);
Ny=length(output.Y);
Nx=length(output.X);
    
if Idformat==1
    if Idvar==1
        save('-v7.3',[savedir,'/',savename,'.mat'],'output')
    elseif Idvar==2
        phidat=dat.output.phi;
        output.phi=phidat(indt1:stept:indt2,indy1:stepy:indy2,indx1:stepx:indx2);
        save ('-v7.3',[savedir,'/',savename,'.mat'],'output')
    elseif Idvar==3
        if isfield(dat.output,'phi')
            phidat=dat.output.phi;
            output.phi=phidat(indt1:stept:indt2,indy1:stepy:indy2,indx1:stepx:indx2);
        end
        output.break_nodes=dat.output.break_nodes;
        output.break_crest=dat.output.break_crest;
        output.break_speed=dat.output.break_speed;
        output.comtime=dat.output.comtime;
        input=dat.input;par=dat.par;model=dat.model;bath=dat.bath;
        dom=dat.dom;ivp=dat.ivp;bdyassim=dat.bdyassim;influx=dat.influx;
        dom.X=output.X;dom.Nx=Nx;
        dom.Y=output.Y;dom.Ny=Ny;
        timeSimul.t_init=output.time(1);timeSimul.t_end=output.time(2);
        timeSimul.dt=output.time(2)-output.time(1);
        timeSimul.interval=output.time;
        timeSimul.Nt=Nt;
        
        if input.bdyassim.option==1
            tempchar=bdyassim.charupdate;bdyassim.charupdate=[];
            bdyassim.charupdate=tempchar(indy1:stepy:indy2,indx1:stepx:indx2);  
        end
        if input.wall.option==1
            tempchar=dom.wall.char;dom.wall.char=[];
            dom.wall.char=tempchar(indy1:stepy:indy2,indx1:stepx:indx2);
        end

        Proj=dat.Proj;timeSimul=dat.timeSimul;
        
        save ('-v7.3',[savedir,'/',savename,'.mat'], ...
                    'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul');
    end
else
    

    outputeta(1,1:3)=[tt(1) tt(end) Nt];
    outputeta(2,1:3)=[xx(indx1) xx(indx2) Nx];
    outputeta(3,1:3)=[yy(indx1) yy(indx2) Ny];
    outputphi(1,1:3)=[tt(1) tt(end) Nt];
    outputphi(2,1:3)=[xx(indx1) xx(indx2) Nx];
    outputphi(3,1:3)=[yy(indx1) yy(indx2) Ny];
   
    Iy=4;
    for i=indt1:stept:indt2
        outputeta(Iy:Iy+Ny-1,1:Nx)=etadat(i,indy1:stepy:indy2,indx1:stepx:indx2);
        if Idvar==2
            phidat=dat.output.phi;
            outputphi(Iy:Iy+Ny-1,1:Nx)=phidat(i,indy1:stepy:indy2,indx1:stepx:indx2);
        end
         Iy=Iy+Ny;
    end
    
     save([savedir,'/',savename,'.txt'],'outputeta','-ascii')
     if Idvar==2
     save([savedir,'/',savename,'.txt'],'outputphi','-ascii')
     end
end
statusbarObj.setText('saved');
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.black);



end

       
end

