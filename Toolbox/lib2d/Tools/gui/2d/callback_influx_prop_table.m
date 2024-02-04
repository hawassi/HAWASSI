function callback_influx_prop_table(hObj,editevent,h)
row = editevent.Indices(1);
col = editevent.Indices(2);
config_data = get(h.waveinput_influx_proptable,'data');
defaultData=get(h.waveinput_influx_proptable,'data');
if col == 1
    [defaultData,userpropData]=modifyPopup1_influx_prop(row,config_data,h,defaultData);
    set(h.waveinput_influx_proptable,'data',defaultData);
    set(h.waveinput_influx_proptable,'userdata',userpropData)
end

if col==2
    newdata = get(h.waveinput_influx_proptable, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'Harmonic') == 1
        defaultData = newdata;
    else
        newdata(row,col) = {'-'};
        defaultData = newdata;
    end
    
    set(h.waveinput_influx_proptable,'data',defaultData);
end


if col == 3 || col == 5 || col == 6
    newdata = get(h.waveinput_influx_proptable, 'data');
    select = config_data(row,1);
    select = select{1};
    if strcmp(select,'Harmonic') == 1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    else
        defaultData = newdata;
    end
    
    set(h.waveinput_influx_proptable,'data',defaultData);
end

if col == 7
    newdata = get(h.waveinput_influx_proptable, 'data');
    select1 = config_data(row,1);
    select1 = select1{1};
    if   strcmp(select1,'User-defined (signal)')
        newdata(row,col) = {'-'};
    end
     defaultData = newdata;
    set(h.waveinput_influx_proptable,'data',defaultData);
end


newdata = get(h.waveinput_influx_proptable, 'data');
select = config_data(row,1);
select = select{1};
if strcmp(select,'User-defined (signal)')
    if col~=1
        newdata(row,col) = {'-'};
        defaultData = newdata;
    end
    
end

set(h.waveinput_influx_proptable,'data',defaultData);
end

function  [defaultData,userpropData]=modifyPopup1_influx_prop(row,config_data,h,defaultData)
newdata = defaultData;
selector = config_data(row,1);
selector = selector{1};

userpropData=get(h.waveinput_influx_proptable,'userdata');

switch selector
    case 'Harmonic'
        newdata(row,1) = {'Harmonic'};
        newdata(row,2) = {0};
        newdata(row,3) = {'-'};
        newdata(row,4) = {0};
        newdata(row,5) = {'-'};
        newdata(row,6) = {'-'};
        newdata(row,7) = {0};
        defaultData=newdata;
    case 'JONSWAP'
        newdata(row,1) = {'JONSWAP'};
        newdata(row,2) = {'-'};
        newdata(row,3) = {0};
        newdata(row,4) = {0};
        newdata(row,5) = {0};
        newdata(row,6) = {0};
        newdata(row,7) = {0};
        defaultData=newdata;
    case 'User-defined (signal)'
        
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load time-signal');
        
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
            
            newdata(row,1) = {'User-defined (signal)'};
            newdata(row,2) = {'-'};
            newdata(row,3) = {'-'};
            newdata(row,4) = {'data'};
            newdata(row,5) = {'loaded'};
            newdata(row,6) = {'-'};
            newdata(row,7) = {'-'};
            defaultData=newdata;
            set(h.monitorbox,'String','>>data loaded','foregroundcolor', 'k');
            
            userpropData(row).inflX=my_data(1,2:end);
            userpropData(row).inflY=my_data(2,2:end);
            userpropData(row).time=my_data(3:end,1);
            userpropData(row).eta=my_data(3:end,2:end);
            set(h.waveinput_influx_proptable, 'userdata', userpropData);
        else
            set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
            newdata(row,1) = {'Harmonic'};
            newdata(row,2) = {0};
            newdata(row,3) = {'-'};
            newdata(row,4) = {0};
            newdata(row,5) = {'-'};
            newdata(row,6) = {'-'};
            newdata(row,7) = {'-'};
            defaultData=newdata;
        end
    case 'User-defined (variance density spectrum)'
        
        [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load variance-density');
        
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
            
            if length(my_data(1,:))~=2
                set(h.monitorbox,'String','>>Wrong input data','foregroundcolor','r');
                newdata(row,1) = {'Harmonic'};
                newdata(row,2) = {0};
                newdata(row,3) = {'-'};
                newdata(row,4) = {0};
                newdata(row,5) = {'-'};
                newdata(row,6) = {'-'};
                newdata(row,7) = {'-'};
                defaultData=newdata;
            else
                if max(any(my_data<0))==1
                    set(h.monitorbox,'String','>>Wrong input data, there is a negative value','foregroundcolor','r');
                    newdata(row,1) = {'Harmonic'};
                    newdata(row,2) = {0};
                    newdata(row,3) = {'-'};
                    newdata(row,4) = {0};
                    newdata(row,5) = {'-'};
                    newdata(row,6) = {'-'};
                    newdata(row,7) = {'-'};
                    defaultData=newdata;
                else
                    newdata(row,1) = {'User-defined (variance density spectrum)'};
                    newdata(row,2) = {'-'};
                    newdata(row,3) = {'data'};
                    newdata(row,4) = {'loaded'};
                    newdata(row,5) = {'-'};
                    newdata(row,6) = {0};
                    newdata(row,7) = {0};
                    defaultData=newdata;
                    set(h.monitorbox,'String','>>data loaded','foregroundcolor', 'k');
                    userpropData(row).varianceDensity=my_data;
                    set(h.waveinput_influx_proptable, 'userdata', userpropData);
                end
            end
        else
            set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
            newdata(row,1) = {'Harmonic'};
            newdata(row,2) = {0};
            newdata(row,3) = {'-'};
            newdata(row,4) = {0};
            newdata(row,5) = {'-'};
            newdata(row,6) = {'-'};
            newdata(row,7) = {'-'};
            defaultData=newdata;
        end
end
set(h.waveinput_influx_proptable, 'data', defaultData);
end

