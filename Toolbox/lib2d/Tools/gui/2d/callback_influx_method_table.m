function callback_influx_method_table(hObj,editevent,h)
        row = editevent.Indices(1);
        col = editevent.Indices(2);
        config_data = get(h.waveinput_influx_methodtable,'data');
        defaultData=get(h.waveinput_influx_methodtable,'data');
        
        if col == 2
            defaultData=modifyPopup1_influx_method(row,config_data,defaultData);
            set(h.waveinput_influx_methodtable,'data',defaultData)
        end
        
        if col==3 || col==4 || col==5 || col==6
            newdata = get(h.waveinput_influx_methodtable, 'data');
            select = config_data(row,2);
            select = select{1};
            if strcmp(select,'Circle Arc') == 1
                newdata(row,col) = {'-'};
                defaultData = newdata;
            else
                defaultData = newdata;
            end 
        set(h.waveinput_influx_methodtable,'data',defaultData);
        end
        
        if col==7 || col==8 || col==9 || col==10 || col==11
            newdata = get(h.waveinput_influx_methodtable, 'data');
            select = config_data(row,2);
            select = select{1};
            if strcmp(select,'Circle Arc') == 1
                defaultData = newdata;
            else
                newdata(row,col) = {'-'};
                defaultData = newdata;
            end 
        set(h.waveinput_influx_methodtable,'data',defaultData);
        end

        
       set(h.waveinput_influx_methodtable,'data',defaultData);
       
       tminnow=get(h.time_edit_interval_init,'userdata');
       tmaxnow=get(h.time_edit_interval_end,'userdata');
       dtnow =get(h.time_edit_step,'userdata');
       if isempty(tminnow)
       set(h.time_edit_interval_init,'userdata',cell2mat(defaultData(row,12)),...
           'string',num2str(cell2mat(defaultData(row,12))));
       else
          tminDat=cell2mat(defaultData(row,12));
          if tminnow>tminDat || row==1
            set(h.time_edit_interval_init,'userdata',tminDat,...
           'string',num2str(tminDat));  
          end
        
           
       end
       if isempty(tmaxnow)
       set(h.time_edit_interval_end,'userdata',cell2mat(defaultData(row,13)),...
           'string',num2str(cell2mat(defaultData(row,13))));
       else
           tmaxDat=cell2mat(defaultData(row,13));
          if tmaxnow<tmaxDat ||row==1
            set(h.time_edit_interval_end,'userdata',tmaxDat,...
           'string',num2str(tmaxDat));  
          end
          
         
       end
       if isempty(dtnow) 
       set(h.time_edit_step,'userdata',cell2mat(defaultData(row,14)),...
           'string',num2str(cell2mat(defaultData(row,14))));
       else
             dtDat=cell2mat(defaultData(row,14));
          if dtnow>dtDat || row==1
            set(h.time_edit_step,'userdata',dtDat,...
           'string',num2str(dtDat));  
          end
          
          
       end
       
       
       
end

function defaultData=modifyPopup1_influx_method(row,config_data,defaultData)
       newdata = defaultData;
       selector = config_data(row,2);
       selector = selector{1};
       switch selector
            case 'Straight'
                newdata(row,3) = {0};
                newdata(row,4) = {0};
                newdata(row,5) = {0};
                newdata(row,6) = {0};
                newdata(row,7) = {'-'};
                newdata(row,8) = {'-'};
                newdata(row,9) = {'-'};
                newdata(row,10) = {'-'};
                newdata(row,11) = {'-'};
                newdata(row,12) = {0};
                newdata(row,13) = {0};
                newdata(row,14) = {0};
                defaultData=newdata;
            case 'Circle Arc'
                newdata(row,3) = {'-'};
                newdata(row,4) = {'-'};
                newdata(row,5) = {'-'};
                newdata(row,6) = {'-'};
                newdata(row,7) = {0};
                newdata(row,8) = {0};
                newdata(row,9) = {0};
                newdata(row,10) = {0};
                newdata(row,11) = {0};
                newdata(row,12) = {0};
                newdata(row,13) = {0};
                newdata(row,14) = {0};
               defaultData=newdata;
       end
end

