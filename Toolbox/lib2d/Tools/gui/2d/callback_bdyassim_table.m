function callback_bdyassim_table(hObj,editevent,h)
row = editevent.Indices(1);
col = editevent.Indices(2);
config_data = get(h.bdyassim_table,'data');
defaultData=get(h.bdyassim_table,'data');
bdyAsimuserdata=get(h.bdyassim_table,'userdata');
if col == 1
    [defaultData,bdyAsimuserdata]=modifyPopup1(row,config_data,h,defaultData,bdyAsimuserdata);
    set(h.bdyassim_table,'data',defaultData);
end

if col == 2 || col == 3 || col==4 || col==5
    newdata = get(h.bdyassim_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'Ring') == 1 || strcmp(select,'Half-ring') == 1
        defaultData = newdata;
    else
        newdata(row,col) = {'-'};
        defaultData = newdata;
    end
    set(h.bdyassim_table,'data',defaultData);
end

if  col == 6 || col == 7 || col == 8 || col== 9
    newdata = get(h.bdyassim_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'Rectangle') == 1
        defaultData = newdata;
    else
        newdata(row,col) = {'-'};
        defaultData = newdata;
    end
    set(h.bdyassim_table,'data',defaultData);
end

if col == 11
    defaultData=modifyPopup2(row,config_data,h,defaultData);
    set(h.bdyassim_table,'data',defaultData);
end

if col == 12
    [defaultData,bdyAsimuserdata]=modifyPopup3(row,config_data,h,defaultData,bdyAsimuserdata);
   
    if ~isempty(bdyAsimuserdata(row).assimdata)
    Nt=bdyAsimuserdata(row).assimdata(1,1);
    time=bdyAsimuserdata(row).assimdata(2,1:Nt);
    tminnow=min(time);
    tmaxnow=max(time);
    dtnow =time(2)-time(1);
    defaultData(row,13)={tminnow};
    defaultData(row,14)={tmaxnow};
    defaultData(row,15)={dtnow};
    end
    set(h.bdyassim_table,'data',defaultData);  
end
    
if col==13 ||col==14||col==15
    
end

set(h.bdyassim_table,'data',defaultData);
set(h.bdyassim_table,'userdata',bdyAsimuserdata)
end

function  [defaultData,bdyAsimuserdata]=modifyPopup1(row,config_data,h,defaultData,bdyAsimuserdata)
newdata = defaultData;
selector = config_data(row,1);
selector = selector{1};

switch selector
    case 'Ring'
        newdata(row,1) = {'Ring'};
        newdata(row,2) = {0};
        newdata(row,3) = {0};
        newdata(row,4) = {0};
        newdata(row,5) = {0};
        newdata(row,6) = {'-'};
        newdata(row,7) = {'-'};
        newdata(row,8) = {'-'};
        newdata(row,9) = {'-'};
        newdata(row,10) = {0};
        newdata(row,11) = {'North'};
        newdata(row,12) = {'-'};
        defaultData=newdata;
    case 'Half ring'
        newdata(row,1) = {'Half-ring'};
        newdata(row,2) = {0};
        newdata(row,3) = {0};
        newdata(row,4) = {0};
        newdata(row,5) = {0};
        newdata(row,6) = {'-'};
        newdata(row,7) = {'-'};
        newdata(row,8) = {'-'};
        newdata(row,9) = {'-'};
        newdata(row,10) = {0};
        newdata(row,11) = {'North'};
        newdata(row,12) = {'-'};
        defaultData=newdata;
    case 'Rectangle'
        newdata(row,1) = {'Rectangle'};
        newdata(row,2) = {'-'};
        newdata(row,3) = {'-'};
        newdata(row,4) = {'-'};
        newdata(row,5) = {'-'};
        newdata(row,6) = {0};
        newdata(row,7) = {0};
        newdata(row,8) = {0};
        newdata(row,9) = {0};
        newdata(row,10) = {0};
        newdata(row,11) = {'North'};
        newdata(row,12) = {'-'};
        defaultData=newdata;
        
    case 'User-defined'
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load a boundary domain');
        if directory~=0
            set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            
            newdata(row,1) = {'User-defined'};
            newdata(row,2) = {'-'};
            newdata(row,3) = {'-'};
            newdata(row,4) = {'data'};
            newdata(row,5) = {'loaded'};
            newdata(row,6) = {'-'};
            newdata(row,7) = {'-'};
            newdata(row,8) = {'-'};
            newdata(row,9) = {'-'};
            newdata(row,10) = {0};
            newdata(row,11) = {'North'};
            newdata(row,12) = {'-'};
            
            defaultData=newdata;
            set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
            bdyAsimuserdata(row).bdry=my_data;
            
        else
            newdata(row,1) = {'Ring'};
            newdata(row,2) = {0};
            newdata(row,3) = {0};
            newdata(row,4) = {0};
            newdata(row,5) = {0};
            newdata(row,6) = {'-'};
            newdata(row,7) = {'-'};
            newdata(row,8) = {'-'};
            newdata(row,9) = {'-'};
            newdata(row,10) = {0};
            newdata(row,11) = {'North'};
            newdata(row,12) = {'-'};
            defaultData=newdata;
            
        end
end
end

function  defaultData=modifyPopup2(row,config_data,h,defaultData)
newdata = defaultData;
selector = config_data(row,11);
selector = selector{1};
switch selector
    case 'North'
        newdata(row,11) = {'North'};
        defaultData=newdata;
    case 'South'
        newdata(row,11) = {'South'};
        defaultData=newdata;
    case 'West'
        newdata(row,11) = {'West'};
        defaultData=newdata;
    case 'East'
        newdata(row,11) = {'East'};
        defaultData=newdata;
end
end

function [defaultData,bdyAsimuserdata]=modifyPopup3(row,config_data,h,defaultData,bdyAsimuserdata)
newdata = defaultData;
selector = config_data(row,12);
selector = selector{1};
switch selector
    case 'load..'
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load assimilation data');
        if directory~=0
            set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            
            newdata(row,12) = {'loaded'};
            
            defaultData=newdata;
            set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
            bdyAsimuserdata(row).assimdata=my_data;
            
        else
            newdata(row,12) = {'-'};
            defaultData=newdata;
            bdyAsimuserdata(row).assimdata=0;
            
        end
    case '-'
        newdata(row,12) = {'-'};
        defaultData=newdata;
        bdyAsimuserdata(row).assimdata=0;
end
end
