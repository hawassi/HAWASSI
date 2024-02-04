function callback_friction_table(hObj,editevent,h)
row = editevent.Indices(1);
col = editevent.Indices(2);
config_data = get(h.bathymetry_table_friction,'data');
defaultData =get(h.bathymetry_table_friction,'data');
fricuserdata=get(h.bathymetry_table_friction,'userdata');
if col == 1
    [defaultData,fricuserdata]=modifyPopup1_damping(row,config_data,h,defaultData,fricuserdata);
    set(h.bathymetry_table_friction,'data',defaultData);
end

if col == 2
    newdata = get(h.bathymetry_table_friction, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else 
        defaultData = newdata;
    end
   set(h.bathymetry_table_friction,'data',defaultData);
end

if col == 3
    newdata = get(h.bathymetry_table_friction, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'data'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end   
set(h.bathymetry_table_friction,'data',defaultData);
end

if col == 4
    newdata = get(h.bathymetry_table_friction, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'loaded'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end   
set(h.bathymetry_table_friction,'data',defaultData);
end

if col == 5
    newdata = get(h.bathymetry_table_friction, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end 
set(h.bathymetry_table_friction,'data',defaultData);
end

if col == 6
newdata = get(h.bathymetry_table_friction, 'data');
    defaultData = newdata;
set(h.bathymetry_table_friction,'data',defaultData);
end

set(h.bathymetry_table_friction,'data',defaultData);
set(h.bathymetry_table_friction,'userdata',fricuserdata);
end

function  [defaultData,fricuserdata]=modifyPopup1_damping(row,config_data,h,defaultData,fricuserdata)
newdata  = defaultData;
selector = config_data(row,1);
selector = selector{1};

switch selector
    case 'Rectangle'
        newdata(row,1) = {'Rectangle'};
        newdata(row,2) = {0};
        newdata(row,3) = {0};
        newdata(row,4) = {0};
        newdata(row,5) = {0};
        newdata(row,6) = {0};
        defaultData=newdata;
    case 'User-defined'
        
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load friction boundary');
        
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
            newdata(row,3) = {'data'};
            newdata(row,4) = {'loaded'};
            newdata(row,5) = {'-'};
            newdata(row,6) = {0};
            defaultData=newdata;
            set(h.monitorbox,'String','>>','foregroundcolor', 'k');
           
            fricuserdata(row).bdry=my_data;
        else
            set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
            newdata(row,1) = {'Rectangle'};
            newdata(row,2) = {0};
            newdata(row,3) = {0};
            newdata(row,4) = {0};
            newdata(row,5) = {0};
            newdata(row,6) = {0};
            defaultData=newdata;
        end
end
end


