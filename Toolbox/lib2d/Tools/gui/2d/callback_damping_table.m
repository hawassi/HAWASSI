function callback_damping_table(hObj,editevent,h)
row = editevent.Indices(1);
col = editevent.Indices(2);
config_data = get(h.damping_table,'data');
defaultData=get(h.damping_table,'data');
dampuserdata=get(h.damping_table,'userdata');
if col == 1
    [defaultData,dampuserdata]=modifyPopup1_damping(row,config_data,h,defaultData,dampuserdata);
    set(h.damping_table,'data',defaultData);
end

if col == 2
    newdata = get(h.damping_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    elseif strcmp(select,'-') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else 
        defaultData = newdata;
    end
   set(h.damping_table,'data',defaultData);
end

if col == 3
    newdata = get(h.damping_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'data'};
        defaultData = newdata;
     elseif strcmp(select,'-') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end   
set(h.damping_table,'data',defaultData);
end

if col == 4
    newdata = get(h.damping_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'loaded'};
        defaultData = newdata;
     elseif strcmp(select,'-') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end   
set(h.damping_table,'data',defaultData);
end

if col == 5
    newdata = get(h.damping_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'User-defined') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
     elseif strcmp(select,'-') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end 
set(h.damping_table,'data',defaultData);
end

if col == 6
    newdata = get(h.damping_table, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'Fourier Bdry.') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
     elseif strcmp(select,'-') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end 
set(h.damping_table,'data',defaultData);
end

set(h.damping_table,'data',defaultData);
set(h.damping_table,'userdata',dampuserdata);
end

function  [defaultData,dampuserdata]=modifyPopup1_damping(row,config_data,h,defaultData,dampuserdata)
newdata  = defaultData;
selector = config_data(row,1);
selector = selector{1};

switch selector
    case 'Fourier Bdry.'
        if row==1
        newdata(row,1) = {'Fourier Bdry.'};
        newdata(row,2) = {0};
        newdata(row,3) = {0};
        newdata(row,4) = {0};
        newdata(row,5) = {0};
        newdata(row,6) = {'-'};
        defaultData=newdata;
        else
        newdata(row,1) = {'-'};
        newdata(row,2) = {'-'};
        newdata(row,3) = {'-'};
        newdata(row,4) = {'-'};
        newdata(row,5) = {'-'};
        newdata(row,6) = {'-'};
        defaultData=newdata;
        end
    case 'User-defined'
        
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load damping boundary');
        
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
           
            dampuserdata(row).bdry=my_data;
        else
            set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
            newdata(row,1) = {'-'};
            newdata(row,2) = {'-'};
            newdata(row,3) = {'-'};
            newdata(row,4) = {'-'};
            newdata(row,5) = {'-'};
            newdata(row,6) = {'-'};
            defaultData=newdata;
        end
        
        case '-'
            newdata(row,1) = {'-'};
            newdata(row,2) = {'-'};
            newdata(row,3) = {'-'};
            newdata(row,4) = {'-'};
            newdata(row,5) = {'-'};
            newdata(row,6) = {'-'};
            defaultData=newdata;
end
end


