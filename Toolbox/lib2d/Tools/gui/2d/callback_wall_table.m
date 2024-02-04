function callback_wall_table(hObj,editevent,h)
        row = editevent.Indices(1);
        col = editevent.Indices(2);
        config_data = get(h.wall_table,'data');
        defaultData=get(h.wall_table,'data');
        walluserdata=get(h.wall_table,'userdata');
        if col == 1
            [defaultData,walluserdata]=modifyPopup1(row,config_data,h,defaultData,walluserdata);
             set(h.wall_table,'data',defaultData);
        end
        
        if col == 2 || col == 3 || col==4
            newdata = get(h.wall_table, 'data');
            select = config_data(row,1);
            select = select{1};
            if strcmp(select,'Rectangle') == 1
                newdata(row,col) = {'-'};
                defaultData = newdata;
            else
                defaultData = newdata;
            end
             set(h.wall_table,'data',defaultData);
        end
        
        if  col == 5 || col == 6 || col == 7 || col== 8
            newdata = get(h.wall_table, 'data');
            select = config_data(row,1);
            select = select{1};
            if strcmp(select,'Circle') == 1
                newdata(row,col) = {'-'};
                defaultData = newdata;
            else
                defaultData = newdata;
            end
             set(h.wall_table,'data',defaultData);
        end
        
        if col == 9
            newdata = get(h.wall_table, 'data');
            select = config_data(row,9);
            select = select{1};
            if strcmp(select,'Energy truncation') == 1
                newdata(row,col) = {'Energy truncation'};
                newdata(row,10) = {'Uniform'};
                newdata(row,11) = {0};
                newdata(row,12) = {'-'};
                defaultData = newdata;
            else
                newdata(row,col) = {'Influxing'};
                defaultData = newdata;
            end
            
        end
        
        if col == 10
            defaultData=modifyPopup8(row,config_data,h,defaultData);
             set(h.wall_table,'data',defaultData);
        end
        
        if col == 11
            newdata = get(h.wall_table, 'data');
            select = config_data(row,10);
            select = select{1};
            if strcmp(select,'Freq. dependent') == 1
                newdata(row,11) = {'-'};
                set(h.wall_table, 'data', newdata);
                defaultData = newdata;
            else
                defaultData = newdata;
            end
             set(h.wall_table,'data',defaultData);
        end
        
        if col == 12
            newdata = get(h.wall_table, 'data');
            select = config_data(row,10);
            select = select{1};
            if strcmp(select,'Uniform') == 1
                newdata(row,12) = {'-'};
                defaultData = newdata;
            else
                defaultData = newdata;
            end
             set(h.wall_table,'data',defaultData);
        end
        
        set(h.wall_table,'data',defaultData);
        set(h.wall_table,'userdata',walluserdata)
    end

    function  [defaultData,walluserdata]=modifyPopup1(row,config_data,h,defaultData,walluserdata)
        newdata = defaultData;
        selector = config_data(row,1);
        selector = selector{1};
        
        switch selector
            case 'Circle'
                newdata(row,1) = {'Circle'};
                newdata(row,2) = {0};
                newdata(row,3) = {0};
                newdata(row,4) = {0};
                newdata(row,5) = {'-'};
                newdata(row,6) = {'-'};
                newdata(row,7) = {'-'};
                newdata(row,8) = {'-'};
                defaultData=newdata;
            case 'Rectangle'
                newdata(row,1) = {'Rectangle'};
                newdata(row,2) = {'-'};
                newdata(row,3) = {'-'};
                newdata(row,4) = {'-'};
                newdata(row,5) = {0};
                newdata(row,6) = {0};
                newdata(row,7) = {0};
                newdata(row,8) = {0};
                defaultData=newdata;
            case 'User-defined'
                
                [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load wall boundary');
                
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
                    defaultData=newdata;
                    set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
                    walluserdata(row).bdry=my_data;
                    
                else
                    newdata(row,1) = {'Circle'};
                    newdata(row,2) = {0};
                    newdata(row,3) = {0};
                    newdata(row,4) = {0};
                    newdata(row,5) = {'-'};
                    newdata(row,6) = {'-'};
                    newdata(row,7) = {'-'};
                    newdata(row,8) = {'-'};
                    set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
                    defaultData=newdata;
                    
                end
        end
    end

    function  defaultData=modifyPopup8(row,config_data,h,defaultData)
        newdata = defaultData;
        WallMethod = defaultData(row,9);
        selector = config_data(row,10);
        selector = selector{1};
        switch selector
            case 'Uniform'
                newdata(row,10) = {'Uniform'};
                newdata(row,11) = {0};
                newdata(row,12) = {'-'};
                defaultData=newdata;
            case 'Freq. dependent'
                
               if strcmpi(WallMethod,'Influxing')
               set(h.monitorbox,'String',['>>'],'foregroundcolor', 'k');
                
                newdata(row,10) = {'Freq. dependent'};
                newdata(row,11) = {'-'};
                newdata(row,12) = {'1-min(w.^2,0.5)'};
                %                 newdata(row,8) = {'Freq. dependent'};
                %                 newdata(row,9) = {'-'};
                %                 newdata(row,10) = {'1-min(w.^2,0.5)'};
                defaultData=newdata;
               else
                set(h.monitorbox,'String',['>> Freq. dep. wall is not available.'],'foregroundcolor', 'k');
                newdata(row,10) = {'Uniform'};
                newdata(row,11) = {0};
                newdata(row,12) = {'-'};
                defaultData=newdata;   
               end
        end
    end
