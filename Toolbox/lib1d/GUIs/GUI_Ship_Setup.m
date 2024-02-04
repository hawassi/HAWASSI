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

function varargout=GUI_Ship_Setup(path,varargin)
hObj = figure('Name', 'Ship setup','Unit', 'normalized',...
    'Position', [0.2 0.3 0.4 0.3], 'NumberTitle',...
    'Off', 'MenuBar', 'none');%, 'Resize', 'off');
h= guihandles(hObj);
h.input=varargin;
h.pathnow=path;
varargout{1}.FlagCloseGui= varargin{1}.FlagCloseGui;

panel          = uipanel('Parent',hObj,'Unit', 'normalized','position',[0 0 1 1]);

h.shipsetup_button_addrow = uicontrol('Style', 'pushbutton', 'String','Add Row',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.7 .85 .1  .07], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.5);

h.shipsetup_button_deleterow = uicontrol('Style', 'pushbutton', 'String','Delete Row',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.8 .85 .1  .07], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.5);

h.shipsetup_button_storedata = uicontrol('Style', 'pushbutton', 'String','Input data',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.8 .08 .1  .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.4,'Fontweight','Bold');
h.shipsetup_checkbox_model_linear = uicontrol('Style', 'checkbox', 'String','linear calculation',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.02 .08 .2  .1], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.4,'Fontweight','Bold','value',1);
h.shipsetup_text_EvModes = uicontrol('Style', 'text', 'String','Evan. modes:',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.22 .06 .12  .09], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.4,'Fontweight','Bold');
h.shipsetup_edit_EvModes = uicontrol('Style', 'edit', 'String','',...
    'HorizontalAlignment', 'center', 'Parent', panel, 'Unit',...
    'Normalized', 'Position', [.35 .1 .1  .05], 'Visible', 'on',...
    'Fontunits','normalized','FontSize', 0.8,'Fontweight','Bold');



select_shipsetupshape = {'Barge'; 'Wedge';'Half circle';'User-defined'};
select_shipsetupmotion = {'Fixed'; 'Heave';'Surge'; 'Pitch'; 'Free'};
select_shipmooring = {'no'; 'yes'};
columnname = {'<html><center />Shape</html>', ...
    '<html><center />Motion</html>',...
    '<html><center />Length</html>', ...
    '<html><center />Draft</html>', ...
    '<html><center />Position <br />x [m]</html>', ...
    '<html><center />Position <br />z [m]</html>', ...
    '<html><center />Orientation <br /> &theta [deg]</html>',...
    '<html><center />Gyradius [m]</html>',...
    '<html><center />CoG <br />z [m] </html>',...
    '<html><center />GM</html>',...
    '<html><center />Mooring </html>',...
    '<html><center />T<sub>n<sub>1  </html>',...
    '<html><center />T<sub>n<sub>3  </html>',...
    '<html><center />T<sub>n<sub>5  </html>',...
    '<html><center />C<sub>d<sub>1  </html>',...
    '<html><center />C<sub>d<sub>3  </html>',...
    '<html><center />C<sub>d<sub>5  </html>'};

defdatshipsetup0 =  {'Barge','Fixed', 0, 0, 0, 0, 0, 0, 0, 0,'no','-','-','-',0,0,0};

defdatshipsetup =  defdatshipsetup0;
shipsetupuserdata= [];

if isfield(varargin{1},'linearCal')
    set(h.shipsetup_checkbox_model_linear,'value',varargin{1}.linearCal)
end
if ~isempty(varargin{1}.data)
    defdatshipsetup =  varargin{1}.data;
   if length(defdatshipsetup(1,:))==6 ||length(defdatshipsetup(1,:))==7 ...
           ||length(defdatshipsetup(1,:))==13 ||length(defdatshipsetup(1,:))==16
        defdatshipsetup(:,7)={[0]};
        defdatshipsetup(:,8)={[0]};
        defdatshipsetup(:,9)={[0]};
        defdatshipsetup(:,10)={[0]};
        defdatshipsetup(:,11)={'no'};
        defdatshipsetup(:,12)={'-'};
        defdatshipsetup(:,13)={'-'};
        defdatshipsetup(:,14)={'-'};
        defdatshipsetup(:,15)={[0]};
        defdatshipsetup(:,16)={[0]};
        defdatshipsetup(:,17)={[0]};
   end
    shipsetupuserdata= varargin{1}.user_shapedata;
end


if isfield(varargin{1},'Evmode')
   set(h.shipsetup_edit_EvModes,'userdata',varargin{1}.Evmode); 
   set(h.shipsetup_edit_EvModes,'string',num2str(varargin{1}.Evmode)); 
end

valEvmodes=get(h.shipsetup_edit_EvModes,'userdata');
if isempty(valEvmodes)
set(h.shipsetup_edit_EvModes,'userdata',1,'string','1');
end


columnformat = {{select_shipsetupshape{:}} {select_shipsetupmotion{:}} [] [] [] [] [] [] [] [] {select_shipmooring{:}} [] [] [] [] [] []};
h.shipsetup_table = uitable('Parent', panel, 'Units', 'normalized', ...
    'Position', [0.02 0.2 0.97 0.6], 'ColumnName', columnname, 'ColumnFormat', columnformat,...
    'ColumnEditable', true,'Fontunits','normalized','FontSize', 0.07,'visible','on',...
    'data',defdatshipsetup,'userdata',shipsetupuserdata);
shipsetuptablesize = getpixelposition(h.shipsetup_table);
wts=shipsetuptablesize(3)/9;
set(h.shipsetup_table,'ColumnWidth',{wts wts wts wts wts wts wts wts 1.2*wts wts 1.2*wts wts wts wts wts wts wts})

set(h.shipsetup_table,'CellEditCallback',{@callback_shipsetup_table,h});
set(h.shipsetup_button_addrow,'callback',{@callback_addrow_table,h,h.shipsetup_table,defdatshipsetup0});
set(h.shipsetup_button_deleterow,'callback',{@callback_deleterow_table,h,h.shipsetup_table});
set(h.shipsetup_button_storedata ,'callback',{@callback_storedata,h,hObj});
set(h.shipsetup_edit_EvModes ,'callback',{@callback_editEvmode,h});


    function callback_storedata(~,~,h,hfig)
        shipparam=get(h.shipsetup_table,'data');
        user_shapedata=get(h.shipsetup_table,'userdata');
        [statusbarObj]=JavaFrame_handling();
        Flag=0;
        Ns=length(shipparam(:,1));
        for ii=1:Ns
            if any(isempty(cell2mat(shipparam(ii,3))))||...
                    any(isempty(cell2mat(shipparam(ii,4))))||...
                    any(isempty(cell2mat(shipparam(ii,5))))||...
                    any(isempty(cell2mat(shipparam(ii,6))))||...
                    any(isempty(cell2mat(shipparam(ii,7))))||...
                    any(isempty(cell2mat(shipparam(ii,8))))||...
                    any(isempty(cell2mat(shipparam(ii,15))))||...
                    any(isempty(cell2mat(shipparam(ii,16))))||...
                    any(isempty(cell2mat(shipparam(ii,17))))
                    
                   % any(isnan(cell2mat(shipparam(ii,3:7))))
                statusbarObj.setText('Specify ship parameter');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                Flag=1;
            elseif any(strcmp(shipparam(:,1),'User-defined')) && isempty(user_shapedata)
                statusbarObj.setText('Specify ship shape data...');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                Flag=1;
            end
            
             if strcmp(shipparam(ii,2),'Pitch') ||strcmp(shipparam(ii,2),'Free')
                if isempty(cell2mat(shipparam(ii,9))) 
                statusbarObj.setText('Specify z center rotation');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                Flag=1;
                end
                if isempty(cell2mat(shipparam(ii,10))) 
                statusbarObj.setText('Specify GM');
                statusbarTxt = statusbarObj.getComponent(0);
                statusbarTxt.setForeground(java.awt.Color.red);
                Flag=1;
                end
             end
             if strcmp(shipparam(ii,11),'yes')
                 try
                     if strcmp(shipparam(ii,2),'Surge')
                         if isempty(cell2mat(shipparam(ii,12))) || isnan(cell2mat(shipparam(ii,12)))
                             Flag=1;
                         end
                     elseif strcmp(shipparam(ii,2),'Heave')
                         if isempty(cell2mat(shipparam(ii,13))) || isnan(cell2mat(shipparam(ii,13)))
                             Flag=1;
                         end
                     elseif strcmp(shipparam(ii,2),'Pitch')
                         if isempty(cell2mat(shipparam(ii,14))) || isnan(cell2mat(shipparam(ii,14)))
                             Flag=1;
                         end
                     else
                         if isempty(cell2mat(shipparam(ii,12))) || isnan(cell2mat(shipparam(ii,12))) ...
                                 ||isempty(cell2mat(shipparam(ii,13))) || isnan(cell2mat(shipparam(ii,13)))...
                                 ||isempty(cell2mat(shipparam(ii,14))) || isnan(cell2mat(shipparam(ii,14)))
                             Flag=1;
                         end
                     end
                     
                 catch
                     if strcmp(shipparam(ii,2),'Surge')
                         if isempty(str2num(cell2mat(shipparam(ii,12)))) || isnan(str2num(cell2mat(shipparam(ii,12))))
                             Flag=1;
                         end
                     elseif strcmp(shipparam(ii,2),'Heave')
                         if isempty(str2num(cell2mat(shipparam(ii,13)))) || isnan(str2num(cell2mat(shipparam(ii,13))))
                             Flag=1;
                         end
                     elseif strcmp(shipparam(ii,2),'Pitch')
                         if isempty(str2num(cell2mat(shipparam(ii,14)))) || isnan(str2num(cell2mat(shipparam(ii,14))))
                             Flag=1;
                         end
                     else
                         if isempty(str2num(cell2mat(shipparam(ii,12)))) || isnan(str2num(cell2mat(shipparam(ii,12))))...
                                 ||isempty(str2num(cell2mat(shipparam(ii,13)))) || isnan(str2num(cell2mat(shipparam(ii,13))))...
                                 ||isempty(str2num(cell2mat(shipparam(ii,14)))) || isnan(str2num(cell2mat(shipparam(ii,14))))
                             Flag=1;
                         end
                     end
                     
                 end
                 if Flag==1
                     statusbarObj.setText('Specify natural periode of the ship after mooring applied');
                     statusbarTxt = statusbarObj.getComponent(0);
                     statusbarTxt.setForeground(java.awt.Color.red);
                 end
             end
        end
  
        if Flag==0
            statusbarObj.setText('');
            
            varargout{1}.linearCal=get(h.shipsetup_checkbox_model_linear,'value');
            varargout{1}.data=shipparam;
            varargout{1}.user_shapedata=user_shapedata;
            varargout{1}.FlagCloseGui=0;
            varargout{1}.Evmode=get(h.shipsetup_edit_EvModes,'userdata');
            delete(hfig);
        end
    end


set(hObj,'CloseRequestFcn',{@closefunction,h})

    function closefunction(hObj,~,h)
        varargout{1}.linearCal=get(h.shipsetup_checkbox_model_linear,'value');
        varargout{1}.data=get(h.shipsetup_table,'data');
        varargout{1}.user_shapedata=get(h.shipsetup_table,'userdata');
        varargout{1}.FlagCloseGui=0;
        varargout{1}.Evmode=get(h.shipsetup_edit_EvModes,'userdata');
        delete(hObj);
    end


waitfor(hObj);
end

 function callback_editEvmode(hObj,~,~)
 str=get(hObj,'string');
 set(hObj,'userdata',str2num(str));
 end

function callback_addrow_table(hObj,eventdata,h,htable,defdata)
oldData = get(htable,'Data');
newData = [oldData; defdata];
set(htable,'Data',newData);
end

function callback_deleterow_table(hObj,eventdata,h,htable)
oldData = get(htable,'Data');
Ndata=length(oldData(:,1));
if Ndata>1
    newData = oldData(1:Ndata-1,:);
    set(htable,'Data',newData);
end
end


function callback_shipsetup_table(hObj,editevent,h)
row = editevent.Indices(1);
col = editevent.Indices(2);
config_data = get(h.shipsetup_table,'data');
defaultData=get(h.shipsetup_table,'data');
shipsetupuserdata=get(h.shipsetup_table,'userdata');
if col == 1
    [defaultData,shipsetupuserdata]=modifyPopup1(row,config_data,h,defaultData,shipsetupuserdata);
    set(h.shipsetup_table,'data',defaultData);
    set(h.shipsetup_table,'userdata',shipsetupuserdata)
end


if col == 2
    [defaultData]=modifyPopup2(row,config_data,h,defaultData,shipsetupuserdata);
    set(h.shipsetup_table,'data',defaultData);
end

if col == 8 || col==9 || col==10
    [defaultData]=modifyPopup3(row,config_data,h,defaultData,shipsetupuserdata);
    set(h.shipsetup_table,'data',defaultData);
end


if col == 11
    [defaultData]=modifyPopup4(row,config_data,h,defaultData,shipsetupuserdata);
    set(h.shipsetup_table,'data',defaultData);
end

if col == 12 || col==13 || col==14
    [defaultData]=modifyPopup5(row,config_data,h,defaultData,shipsetupuserdata);
    set(h.shipsetup_table,'data',defaultData);
end




end

function  [defaultData,shipsetupuserdata]=modifyPopup1(row,config_data,h,defaultData,shipsetupuserdata)
newdata = defaultData;
selector = config_data(row,1);
selector = selector{1};

switch selector
    case 'User-defined'
        
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load shipsetup boundary');
        
        if directory~=0
            %             set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            
            newdata(row,1) = {'User-defined'};
            defaultData=newdata;
            %             set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
            shipsetupuserdata(row).bdry=my_data;
            
        else
            newdata(row,1) = {'Barge'};
            %             set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
            defaultData=newdata;
            
        end
end
end

function  [defaultData]=modifyPopup2(row,config_data,h,defaultData,shipsetupuserdata)
newdata = defaultData;
selector = config_data(row,2);
selector = selector{1};

MooredId=config_data(row,11);
switch selector
    case 'Pitch'
        newdata(row,8) = {[0]};
        newdata(row,9) = {[0]};
        newdata(row,10) = {[0]};
        
        if strcmpi(MooredId,'yes')
         newdata(row,14) = {[]};   
         newdata(row,13) = {'-'};
         newdata(row,12) = {'-'};
        end
        newdata(row,15) = {'-'};
        newdata(row,16) = {'-'};
        newdata(row,17) = {[0]};
    case 'Free'
            newdata(row,8) = {[0]};
            newdata(row,9) = {[0]};
            newdata(row,10) = {[0]};
            
            if strcmpi(MooredId,'yes')
                newdata(row,14) = {[]};
                newdata(row,13) = {[]};
                newdata(row,12) = {[]};
            end
             newdata(row,15) = {[0]};
             newdata(row,16) = {[0]};
             newdata(row,17) = {[0]};
    case 'Fixed'
            newdata(row,8) = {'-'};
            newdata(row,9) = {[0]};
            newdata(row,10) = {'-'};
            
                newdata(row,11) = {'no'};
                newdata(row,14) = {'-'};
                newdata(row,13) = {'-'};
                newdata(row,12) = {'-'};
                 newdata(row,15) = {'-'};
                 newdata(row,16) = {'-'};
                 newdata(row,17) = {'-'};
    case 'Heave'
            newdata(row,8) = {'-'}; 
             newdata(row,9) = {'-'};
             newdata(row,10) = {'-'};
            if strcmpi(MooredId,'yes')
                newdata(row,14) = {'-'};
                newdata(row,13) = {[]};
                newdata(row,12) = {'-'};
            else
                newdata(row,14) = {'-'};
                newdata(row,13) = {'-'};
                newdata(row,12) = {'-'};
            end
                 newdata(row,15) = {'-'};
                 newdata(row,16) = {[0]};
                 newdata(row,17) = {'-'};
    case 'Surge'
        newdata(row,8)    = {'-'};
         newdata(row,9) = {'-'};
          newdata(row,10) = {'-'};
        if strcmpi(MooredId,'yes')
            newdata(row,14) = {'-'};
            newdata(row,13) = {'-'};
            newdata(row,12) = {[]};
        else
                newdata(row,14) = {'-'};
                newdata(row,13) = {'-'};
                newdata(row,12) = {'-'};
        end
         newdata(row,15) = {[0]};
         newdata(row,16) = {'-'};
         newdata(row,17) = {'-'};
end
  defaultData=newdata;
end

function  [defaultData]=modifyPopup3(row,config_data,h,defaultData,shipsetupuserdata)
newdata = defaultData;
selector = config_data(row,2);
selector = selector{1};

switch selector
    case 'Heave'
            newdata(row,8) = {'-'}; 
            newdata(row,9) = {'-'}; 
            newdata(row,10) = {'-'}; 
    case 'Surge'
            newdata(row,8) = {'-'}; 
            newdata(row,9) = {'-'}; 
            newdata(row,10) = {'-'}; 
            
            
end
  defaultData=newdata;
end

function  [defaultData]=modifyPopup4(row,config_data,h,defaultData,shipsetupuserdata)
newdata = defaultData;
selector = config_data(row,11);
selector = selector{1};

selector2 = config_data(row,2);

switch selector
    case 'no'
            newdata(row,12) = {'-'};
            newdata(row,13) = {'-'};
            newdata(row,14) = {'-'};
    case 'yes'
        if strcmpi(selector2, 'Surge')
             newdata(row,12) = {[0]};
             newdata(row,13) = {'-'};
             newdata(row,14) = {'-'};
        elseif strcmpi(selector2, 'Heave')
             newdata(row,12) = {'-'};
             newdata(row,13) = {[0]};
             newdata(row,14) = {'-'};
        elseif strcmpi(selector2, 'Pitch')
             newdata(row,12) = {'-'};
             newdata(row,13) = {'-'};
             newdata(row,14) = {[0]}; 
        else
             newdata(row,12) = {[0]};
             newdata(row,13) = {[0]};
             newdata(row,14) = {[0]};    
        end
        
end
  defaultData=newdata;
end

function  [defaultData]=modifyPopup5(row,config_data,h,defaultData,shipsetupuserdata)
newdata = defaultData;
selector = config_data(row,11);
selector = selector{1};
MotionId = config_data(row,2);

switch selector
    case 'no'
            newdata(row,12) = {'-'};
            newdata(row,13) = {'-'};
            newdata(row,14) = {'-'};
    case 'yes'
        if strcmpi(MotionId,'Surge')
            newdata(row,13) = {'-'};
            newdata(row,14) = {'-'};
        elseif strcmpi(MotionId,'Heave')
            newdata(row,12) = {'-'};
            newdata(row,14) = {'-'};
        elseif strcmpi(MotionId,'Pitch')
            newdata(row,12) = {'-'};
            newdata(row,13) = {'-'};
        elseif strcmpi(MotionId,'Fixed')
            newdata(row,12) = {'-'};
            newdata(row,13) = {'-'};
            newdata(row,14) = {'-'};
        end
        
end
  defaultData=newdata;
end


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
end

