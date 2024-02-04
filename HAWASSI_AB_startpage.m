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
% !   HAWASSI AB startpage
% !
% !--------------------------------------------------------------------------------------
% !

function HAWASSI_AB_startpage

h.fig = figure('Name', 'Start Page','Unit', 'normalized',...
    'Position', [0.2 0.3 0.5 0.4], 'NumberTitle', 'Off', 'MenuBar', 'none');

if ~isdeployed
    [pathstr,~,~] = fileparts(mfilename('fullpath'));
    cd (pathstr);
    addpath(genpath(pathstr))
    h.pathnow=pathstr;
else
    projdir= getenv('USERPROFILE');
    h.pathnow=projdir;
end

userfold=getenv('LOCALAPPDATA');
projhistfile=[userfold,'\projhistAB.mat'];
if exist(projhistfile,'file')
    load(projhistfile, 'project');
    FlagEmptyProj=0;
    if isempty(project)
    project={'','','',userfold,''};   
    FlagEmptyProj=1;
    end
else
    project={'','','',userfold,'' };
end

h.projecthist = project;  

Panel = uipanel('Parent', h.fig);
try
    AxPic = axes('Parent',Panel, 'Position', [0.01 0.7 .96 .28], 'Unit', 'normalized');
    if ~isdeployed
        [pict,map] = imread([h.pathnow,'\Toolbox\lib2d\Tools\gui\images\logo\hawassi_black.jpg']);
    else
        [pict,map] = imread('hawassi_black.jpg');        
    end
imshow(pict, map, 'Parent', AxPic, 'XData',[0 1000], 'YData',[200 400])    
catch
end

h.Open_text=uicontrol('Style','text','String','Open an Existing Project',...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Parent', Panel, ...
    'Unit', 'Normalized', 'Position', [0.025 .62 .4 .05],'fontunits','normalized',...
    'fontsize',0.5);

h.table_hist= uitable('Parent', Panel, 'Unit', ...
    'normalized', 'Position', [0.025 0.27 .95 .35], 'ColumnName',...
    {'Name','Created','Modified','Location'}, 'RowName', ([]),'fontunits','normalized',...
    'fontsize',0.1);

if exist(projhistfile,'file')
if FlagEmptyProj==0
set(h.table_hist,'data', h.projecthist(:,1:4));
else
set(h.table_hist,'data', {'','','','',''});    
end
else
set(h.table_hist,'data', {'','','','',''});    
end

if exist(projhistfile,'file')
    Nsaved=20;
    if length(h.projecthist(:,1))>Nsaved
        set(h.table_hist,'data', h.projecthist(length(h.projecthist)-Nsaved-1:end,1:4));
    end
end

TableSize = getpixelposition(h.table_hist);
set(h.table_hist,'ColumnWidth', {1/5*TableSize(3), 1/5*TableSize(3), 1/5*TableSize(3), 2/5*TableSize(3)});
h.desc_box = uicontrol('Style', 'edit', 'String', '', 'Unit', 'Normalized',...
    'BackgroundColor', 'w', 'Position', [0.028 0.14 .945 .11], 'Enable', 'off',...
    'fontunits','normalized','fontsize',0.3);
h.NPbutton = uicontrol('Style','pushbutton','String','New Project', ...
    'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized',...
    'Position', [0.56 .06 .13 .05],'fontunits','normalized','fontsize',0.5);
h.OPbutton = uicontrol('Style','pushbutton','String','Open Project', ...
    'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized',...
    'Position', [0.7 .06 .13 .05],'fontunits','normalized','fontsize',0.5,...
    'Enable', 'off');
h.DPbutton = uicontrol('Style','pushbutton','String','Delete Project', ...
    'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized',...
    'Position', [0.84 .06 .13 .05],'fontunits','normalized','fontsize',0.5,...
    'Enable', 'off');
try
    if ~isdeployed
        [icon_Browse,~]= imread([h.pathnow,'\Toolbox\lib2d\Tools\gui\images\folder.png']);
    else
        [icon_Browse,~]= imread('folder.png');
    end
    
    %browse transparency
    icon_Browse = double(icon_Browse)/255;
    index1 = icon_Browse(:,:,1) == 1;
    index2 = icon_Browse(:,:,2) == 1;
    index3 = icon_Browse(:,:,3) == 1;
    
    indexWhite = index1+index2+index3==3;
    
    for idx = 1 : 3
        rgb = icon_Browse(:,:,idx);     % extract part of the image
        rgb(indexWhite) = NaN;          % set the white portion of the image to NaN
        icon_Browse(:,:,idx) = rgb;     % substitute the update values
    end
    
    h.Browsebutton = uicontrol('Style','pushbutton','String','', ...
        'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized',...
        'Position', [0.9 .64 .035 .07],'cdata', icon_Browse);
catch
    h.Browsebutton = uicontrol('Style','pushbutton','String','Browse', ...
        'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized',...
        'Position', [0.88 .64 .1 .07]);
end
h.info=uimenu(h.fig,'Label','Info');
h.licenceCheck=uimenu(h.info,'Label','Check licence');
h.licenceActivate=uimenu(h.info,'Label','Input licence');
    
[IdLicence,h.InfoLicence]=licence_check_AB(h,1);

if IdLicence==0
 set(h.DPbutton,'enable','off');
 set(h.NPbutton,'enable','off');
 set(h.OPbutton,'enable','off');
 set(h.Browsebutton,'enable','off');
else
 set(h.DPbutton,'enable','on');
 set(h.NPbutton,'enable','on');
 set(h.OPbutton,'enable','on');
 set(h.Browsebutton,'enable','on');  
end

set(h.fig,'ResizeFcn',{@f1ResizeFcn,h});
set(h.table_hist,'CellSelectionCallback',{@tableSelection,h})
set(h.OPbutton, 'Callback', {@OPbutton_callback,h})
set(h.NPbutton, 'Callback', {@NPbutton_callback,h})
set(h.DPbutton, 'Callback', {@DPbutton_callback,h})
set(h.Browsebutton, 'Callback', {@Browsebutton_callback,h})
set(h.licenceActivate,'Callback',{@licence_activate_callback,h})
set(h.licenceCheck,'Callback',{@licence_check_callback,h})

guidata(h.fig, h);
end

function licence_check_callback(hObj,eventd,h)

end

function licence_activate_callback(hObj,eventd,h)

end

function f1ResizeFcn(~,~,h)
TableSize = getpixelposition(h.table_hist);
set(h.table_hist,'ColumnWidth', {1/5*TableSize(3), 1/5*TableSize(3), 1/5*TableSize(3), 2/5*TableSize(3)});
end

function NPbutton_callback(~,~,h)
close(h.fig);
output.pathnow=h.pathnow;
userfold=getenv('LOCALAPPDATA');
projhistfile=[userfold,'\projhistAB.mat'];

if exist(projhistfile,'file')
    load(projhistfile, 'project');
    if isempty(project)
    project={'','','',userfold,''};
    end  
else
    project={'','','',userfold,''};
end

h.projecthist=project;
output.projecthist=h.projecthist;
output.InfoLicence=h.InfoLicence;
GUI_Newproject(output);
end

function tableSelection(hObj,eventdata,h)
index = eventdata.Indices;
[IdLicence,~]=licence_check_AB(h,0);

if any(index)
    h.tablehist_row = index(1);
    h.tablehist_col = index(2);
    
    set(hObj, 'userdata',h.tablehist_row);
    notes = h.projecthist(h.tablehist_row,5);
    notes = notes{1};
    
    if ~isempty(cell2mat(h.projecthist(h.tablehist_row,4)))
        set(h.desc_box, 'String', notes, 'Enable','on');
        set(h.OPbutton, 'Enable', 'on');
        set(h.DPbutton, 'Enable', 'on');
    end
end
if IdLicence==0
 set(h.OPbutton,'enable','off');
 set(h.DPbutton,'enable','off');
 set(h.NPbutton,'enable','off');
 set(h.OPbutton,'enable','off');
 set(h.Browsebutton,'enable','off');
else
 set(h.OPbutton,'enable','on');
 set(h.DPbutton,'enable','on');
 set(h.NPbutton,'enable','on');
 set(h.OPbutton,'enable','on');
 set(h.Browsebutton,'enable','on');  
end


guidata(h.fig, h);
end

function OPbutton_callback(hObj,eventdata,h)
h=guidata(hObj);
[statusbarObj]=funGui_JavaFrame_handling();
statusbarObj.setText('Loading a project. . .');

userfold=getenv('LOCALAPPDATA');
projhistfile=[userfold,'\projhistAB.mat'];
if exist(projhistfile,'file')
    load(projhistfile, 'project');
else
    project={'','','',userfold,''};
end
IdClose=1;
try 
    datahist = get(h.table_hist,'Data');
    try
    load([cell2mat(datahist(h.tablehist_row,4)),'abproj_',cell2mat(datahist(h.tablehist_row,1)),'.mat']); 
    catch
    load([cell2mat(datahist(h.tablehist_row,4)),'abprojSA_',cell2mat(datahist(h.tablehist_row,1)),'.mat']); 
    end
    statusbarObj.setText('');
    statusbarTxt = statusbarObj.getComponent(0);
catch
    IdClose=0;
    statusbarObj.setText('file is not found');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

h.projecthist=project;
InfoLicence=h.InfoLicence;
output.flagOpenProj=1;
output.pathnow=h.pathnow;
output.projecthist=h.projecthist;

output.module=GUIinput.proj.module;
output.usernote=GUIinput.proj.note;
output.projectdirectory=GUIinput.proj.projdir;
output.projectname=GUIinput.proj.name;
output.GUIinput=GUIinput;


statusbarObj.setText('');

if ~strcmp(InfoLicence.module,'wave 1D & 2D') && ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
if strcmp(output.module,'AB 2D - Wave') && ~strcmp(InfoLicence.module,'wave 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'AB 1D - Wave') && ~strcmp(InfoLicence.module,'wave 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && (strcmp(InfoLicence.module,'wave 1D')||strcmp(InfoLicence.module,'wave 2D'))
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
end

if ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
if strcmp(output.module,'AB 2D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if strcmp(output.module,'AB 1D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
end

output.IdDemo=InfoLicence.demo;
output.restrict=InfoLicence.restrict;


statusbarObj.setText(['>>Opening a gui, please wait...']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.blue);
  
if strcmp(output.module,'AB 2D - Wave')
    output.ShipFlag=0;
    GuiMain2dWaveModule(output);
elseif strcmp(output.module,'AB 1D - Wave')
    output.ShipFlag=0;
    HaWaSSI_main(output);
elseif strcmp(output.module,'AB 2D - Wave Ship')
    output.ShipFlag=1;
    GuiMain2dWaveModule(output);
elseif strcmp(output.module,'AB 1D - Wave Ship')
    output.ShipFlag=1;
    HaWaSSI_main(output);
elseif strcmp(output.module,'Analytic - Wave Ship')
    GUIMainShipAnalytic(output);
else
    statusbarObj.setText('>>Radar module is not available yet.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    IdClose=0;
end
statusbarObj.setText('');
statusbarTxt = statusbarObj.getComponent(0);

if IdClose==1
close(h.fig);
end

end

function DPbutton_callback(hObj,eventdata,h)
h=guidata(hObj);
[statusbarObj]=funGui_JavaFrame_handling();
data = get(h.table_hist,'Data');
rows = get(h.table_hist,'UserData');

if ~isempty(rows)
  if ~isempty(data)&& rows<= length(data(:,1))
    temp = cell2mat(data(rows,1));
    mask = (1:length(data(:,1)))';
    mask(rows) = [];
    data = data(mask,:);
    h.projecthist=h.projecthist(mask,:);
    project = h.projecthist;
    userfold=getenv('LOCALAPPDATA');
    projhistfile=[userfold,'\projhistAB.mat'];
    
    if exist(projhistfile,'file')
        set(h.table_hist,'Data',data);
        save(projhistfile,'project')
    else
        set(h.table_hist,'Data', {'','','',''});
    end
    
    
    statusbarObj.setText(sprintf('Project %s has been deleted',temp));
    statusbarTxt = statusbarObj.getComponent(0);
    set(h.desc_box, 'String', '', 'Enable','off');
    h.tablehist_row=[];
    guidata(h.fig, h);
  else
    statusbarObj.setText(sprintf('Please select a project'));
    statusbarTxt = statusbarObj.getComponent(0);
  end
else
    statusbarObj.setText(sprintf('Please select a project'));
    statusbarTxt = statusbarObj.getComponent(0);
end
end

function Browsebutton_callback(hObj,eventdata,h)
[statusbarObj]=funGui_JavaFrame_handling();
statusbarObj.setText('Loading a project. . .');

userfold=getenv('LOCALAPPDATA');
projhistfile=[userfold,'\projhistAB.mat'];
if exist(projhistfile,'file')
    load(projhistfile, 'project');
    if isempty(project)
    project={'','','',userfold,''};   
    end
else
    project={'','','',userfold,''};
end

projhist=project;
projdir=cell2mat(projhist(1,4));

if ~isempty(cell2mat(projhist(1,1)))
Ns=length(cell2mat(projhist(1,1)))+1;
projdir=projdir(1:end-Ns);
end

[file_name,directory]=uigetfile([projdir,'\','*.mat'],'Load a project file');

if directory~=0
    load([directory,file_name]);
    if exist('state','var')
        if strcmp(state.gui,'HaWaSSI_Main') %adjusted for AB1D
         GUIinput=state;
         GUIinput.proj.name=state.SaveName;
         GUIinput.proj.workdir='--';
         GUIinput.proj.projdir='--';
         GUIinput.proj.note=state.Note;
         GUIinput.proj.module='AB 1D - Wave';
        end
    end
    if exist('GUIinput','var')
        if ~isfield(GUIinput,'proj')
            statusbarObj.setText('wrong project file.');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return;
        end
    else
        statusbarObj.setText('wrong project file.');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    statusbarObj.setText('Project file loaded.');
else
    % the user clicked the Cancel button.
    % user should provide project data!!
    % otherwise it will not continue
    statusbarObj.setText('No project file loaded.');
    return;
end
GUIinput.proj.workdir='--';
GUIinput.proj.projdir='--';
         
output.flagOpenProj=2;
output.pathnow=h.pathnow;
output.projecthist=projhist;
output.module=GUIinput.proj.module;
output.usernote=GUIinput.proj.note;
output.projectdirectory=GUIinput.proj.projdir;
output.projectname=GUIinput.proj.name;
output.GUIinput=GUIinput;

statusbarObj.setText('');
IdClose=1;
InfoLicence=h.InfoLicence;
if ~strcmp(InfoLicence.module,'wave 1D & 2D') && ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
if strcmp(output.module,'AB 2D - Wave') && ~strcmp(InfoLicence.module,'wave 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if strcmp(output.module,'AB 1D - Wave') && ~strcmp(InfoLicence.module,'wave 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && (strcmp(InfoLicence.module,'wave 1D')||strcmp(InfoLicence.module,'wave 2D'))
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
end

if ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
if strcmp(output.module,'AB 2D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end

if strcmp(output.module,'AB 1D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
if strcmp(output.module,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
    statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return;
end
end

statusbarObj.setText(['>>Opening a gui, please wait...']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.blue);

output.IdDemo=InfoLicence.demo;
output.restrict=InfoLicence.restrict;

if strcmp(output.module,'AB 2D - Wave')
    output.ShipFlag=0;
    GuiMain2dWaveModule(output);
elseif strcmp(output.module,'AB 1D - Wave')
    output.ShipFlag=0;
    HaWaSSI_main(output);
elseif strcmp(output.module,'AB 2D - Wave Ship')
    output.ShipFlag=1;
    GuiMain2dWaveModule(output);
elseif strcmp(output.module,'AB 1D - Wave Ship')
    output.ShipFlag=1;
    HaWaSSI_main(output);
elseif strcmp(output.module,'Analytic - Wave Ship')
    GUIMainShipAnalytic(output);

else
    statusbarObj.setText('>>Radar module is not available yet.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    IdClose=0;
end
if IdClose==1
close(h.fig);
end

end

%%%%%%%%%%%%%%%%%%%%GUI new project %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUI_Newproject(inputdata)
h.fig = figure('Name', 'New Project','Unit', 'normalized',...
    'Position', [0.2 0.3 0.4 0.4], 'NumberTitle',...
    'Off', 'MenuBar', 'none');%, 'Resize', 'off');
h.inputdata=inputdata;

Panel = uipanel('Parent', h.fig);
h.producttypes_text = uicontrol('Style','text','String','Product Types:', ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Parent', Panel ,...
    'Unit', 'Normalized', 'Position', [0.02 .9 .4 .05],'fontunits','normalized',...
    'fontsize',0.5);
h.modules_text = uicontrol('Style','text','String','Modules:', 'FontWeight',...
    'bold', 'HorizontalAlignment', 'left', 'Parent', Panel , 'Unit', ...
    'Normalized', 'Position', [0.55 .9 .3 .05],'fontunits','normalized',...
    'fontsize',0.5);
h.userselection = uicontrol('Style','edit','String','AB 1D - Wave', ...
    'HorizontalAlignment', 'left', 'Parent', Panel , 'Unit', 'Normalized',...
    'Position', [0.02 .45 .935 .05],'fontunits','normalized',...
    'fontsize',0.5);
h.okbutton = uicontrol('Style','pushbutton','String','OK', 'HorizontalAlignment',...
    'center', 'Parent', Panel , 'Unit', 'Normalized', 'Position', [0.68 .06 .12 .05],...
    'fontunits','normalized','fontsize',0.5);
h.cancelbutton = uicontrol('Style','pushbutton','String','Cancel',...
    'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized', ...
    'Position', [0.835 .06 .12 .05],'fontunits','normalized','fontsize',0.5);
h.productlist = uicontrol('Style', 'listbox', 'String', {'AB 1D','AB 2D','Analytic'}, ...
    'Unit', 'Normalized', 'Position', [0.022 .55 .4 .35], 'Value', 1,...
    'fontunits','normalized','fontsize',0.1);
h.modulelist = uicontrol('Style', 'listbox', 'String', {'Wave', 'Wave Ship', ...
    'Wave Radar'}, 'Unit', 'Normalized', 'Position', [0.55 .55 .4 .35], ...
    'Value', 1,'fontunits','normalized','fontsize',0.1);
h.filename_text = uicontrol('Style','text','String','Project Name ', ...
    'HorizontalAlignment', 'left', 'Parent', Panel , 'Unit', 'Normalized',...
    'Position', [0.02 .35 .4 .05],'fontunits','normalized','fontsize',0.5);
h.filename_edit = uicontrol('Style','edit','String','', 'HorizontalAlignment',...
    'left', 'Parent', Panel ,  'BackgroundColor', 'w', 'Unit', 'Normalized',...
    'Position', [0.125 .36 .83 .05],'fontunits','normalized','fontsize',0.5);
h.workingdir_text = uicontrol('Style','text','String','Directory  ', 'HorizontalAlignment',...
    'left', 'Parent', Panel , 'Unit', 'Normalized', 'Position', [0.02 .26 .3 .05],...
    'fontunits','normalized','fontsize',0.5);
h.workingdir_edit = uicontrol('Style','edit','String','', 'HorizontalAlignment', 'left', ...
    'Parent', Panel ,  'BackgroundColor', 'w', 'Unit', 'Normalized', 'Position', ...
    [0.125 .27 .695 .05],'enable','off','fontunits','normalized','fontsize',0.5);
h.finddirbutton = uicontrol('Style','pushbutton','String','browse', ...
    'HorizontalAlignment', 'center', 'Parent', Panel , 'Unit', 'Normalized', 'Position',...
    [0.85 .27 .1 .05],'fontunits','normalized','fontsize',0.5);
h.usernote_text = uicontrol('Style','text','String','Note        ', ...
    'HorizontalAlignment', 'left', 'Parent', Panel , 'Unit', 'Normalized', ...
    'Position', [0.02 .16 .3 .05],'fontunits','normalized','fontsize',0.5);
h.usernote_edit = uicontrol('Style','edit','String','', 'HorizontalAlignment', ...
    'left', 'Parent', Panel ,  'BackgroundColor', 'w', 'Unit', 'Normalized',...
    'Position', [0.125 .17 .83 .05],'fontunits','normalized','fontsize',0.5);

set(h.okbutton,'Callback', {@OKbutton_callback,h});
set(h.cancelbutton,'Callback',{@cancelbutton_callback,h});
set(h.productlist, 'Callback',{@productlist_callback,h})
set(h.modulelist, 'Callback',{@modulelist_callback,h});
set(h.finddirbutton,'callback',{@finddir_callback,h});
set(h.usernote_edit,'callback',{@usernote_callback,h});
set(h.workingdir_edit,'callback',{@workingdir_callback,h});
set(h.filename_edit,'callback',{@projectname_callback})
guidata(h.fig,h);
end

function productlist_callback(hObj,~,h)
prod = get(h.productlist, 'Value');
prods = get(h.productlist, 'String');
module_ind = get(h.modulelist, 'Value');
module = get(h.modulelist, 'String');

if prod==3
set(h.modulelist, 'Value',2); module_ind=2;   
end

set(h.userselection, 'String', strcat(prods(prod),{' - '}, module(module_ind)));
end

function modulelist_callback(hObj,~,h)
[statusbarObj]=funGui_JavaFrame_handling();

prod = get(h.productlist, 'Value');
prods = get(h.productlist, 'String');
module_ind = get(h.modulelist, 'Value');
module = get(h.modulelist, 'String');
statusbarObj.setText(['']);
if module_ind==3
set(h.modulelist, 'Value',1);  module_ind=1;
statusbarObj.setText(['>>Radar module is not available']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.red);
end

if module_ind==2 && prod==2
set(h.modulelist, 'Value',1); module_ind=1;
statusbarObj.setText(['>> AB2D Wave Ship is not available']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.red);    
end

if prod==3
set(h.modulelist, 'Value',2); module_ind=2;
statusbarObj.setText(['>> Analytic code is for wave-ship module']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.red);    
end

set(h.userselection, 'String', strcat(prods(prod),{' - '}, module(module_ind)));
end

function cancelbutton_callback(~,~,h)
close(h.fig);
HAWASSI_AB_startpage;
end

function OKbutton_callback(~,~,h)
[statusbarObj]=funGui_JavaFrame_handling();

UserSelect=get(h.userselection,'string');

output.flagOpenProj=0;
output.pathnow=h.inputdata.pathnow;
output.projecthist=h.inputdata.projecthist;
output.module=UserSelect;
output.usernote=get(h.usernote_edit,'string');
output.projectdirectory=get(h.workingdir_edit,'string');
output.projectname=get(h.filename_edit,'string');
InfoLicence=h.inputdata.InfoLicence;
if strcmp(output.projectname,'')
    statusbarObj.setText('>>Please input projectname');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(h.filename_edit);
    return;
end

if strcmp(output.projectdirectory,'')
    statusbarObj.setText('>>Please browse a working directory');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(h.finddirbutton);
    return;
end
IdClose=1;
statusbarObj.setText('');

if ~strcmp(InfoLicence.module,'wave 1D & 2D') && ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
    if strcmp(UserSelect,'AB 2D - Wave') && ~strcmp(InfoLicence.module,'wave 2D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if strcmp(UserSelect,'AB 1D - Wave') && ~strcmp(InfoLicence.module,'wave 1D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
     if strcmp(UserSelect,'Analytic - Wave Ship') && (strcmp(InfoLicence.module,'wave 1D')||strcmp(InfoLicence.module,'wave 2D'))
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
     end
    
end

if ~strcmp(InfoLicence.module,'wave ship 1D & 2D')
    if strcmp(UserSelect,'AB 2D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    
    if strcmp(UserSelect,'AB 1D - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
    if strcmp(UserSelect,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 1D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
    end
     if strcmp(UserSelect,'Analytic - Wave Ship') && ~strcmp(InfoLicence.module,'wave ship 2D')
        statusbarObj.setText(['>>The licence is only valid for ',InfoLicence.module]);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return;
     end
end
    
statusbarObj.setText(['>>Opening a gui, please wait...']);
statusbarTxt = statusbarObj.getComponent(0);
statusbarTxt.setForeground(java.awt.Color.blue);

output.IdDemo=InfoLicence.demo;
output.restrict=InfoLicence.restrict;

if strcmp(UserSelect,'AB 2D - Wave')
    output.ShipFlag=0;
    GuiMain2dWaveModule(output);
elseif strcmp(UserSelect,'AB 1D - Wave')
    output.ShipFlag=0;
    HaWaSSI_main(output);
elseif strcmp(UserSelect,'AB 2D - Wave Ship')
    output.ShipFlag=1;
    GuiMain2dWaveModule(output);
elseif strcmp(UserSelect,'AB 1D - Wave Ship')
    output.ShipFlag=1;
    HaWaSSI_main(output);
elseif strcmp(UserSelect,'Analytic - Wave Ship')
    GUIMainShipAnalytic(output);
else
    statusbarObj.setText('>>Radar module is not available yet.');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    IdClose=0;
end
if IdClose==1
close(h.fig);
end

end

function projectname_callback(~,~)
[statusbarObj]=funGui_JavaFrame_handling();
statusbarObj.setText(['']);
end

function finddir_callback(~,~,h)
[statusbarObj]=funGui_JavaFrame_handling();

filename=get(h.filename_edit,'string');

projhist=h.inputdata.projecthist;

projdir=cell2mat(projhist(1,4));
if ~isempty(cell2mat(projhist(1,1)))
Ns=length(cell2mat(projhist(1,1)))+1;
projdir=projdir(1:end-Ns);
end

if strcmp(filename,'')
    statusbarObj.setText(['>>Please input project name']);
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    uicontrol(h.filename_edit);
    return;
else
    workingdir = uigetdir(projdir,'browse a directory');
end
if workingdir ~= 0
    set(h.workingdir_edit,'String',workingdir);
    savedir=[workingdir,'/',filename];
    if isdir(savedir)
        statusbarObj.setText(['>>Warning: Project exists already', ...
            ' it will be overwritten!']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        uicontrol(h.filename_edit);
    end
end
end

function usernote_callback(hObj,~,h)
param=get(hObj,'string');
set(h.usernote_edit,'string',param);
end













