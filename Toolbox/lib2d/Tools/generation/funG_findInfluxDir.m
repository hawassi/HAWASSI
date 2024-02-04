function InflDir=funG_findInfluxDir(inflLine,angle)
if strcmpi(inflLine,'Vertical')
    if (angle >=90 && angle<=270) || (angle <=-90 && angle>=-270)
      InflDir={'UniL'};  
    else
      InflDir={'UniR'};  
    end
elseif strcmpi(inflLine,'Horizontal')
     if (angle >=180 && angle<=360) || (angle <=0 && angle>=-180)
      InflDir={'UniB'};  
    else
      InflDir={'UniT'};  
    end 
end
   
end

    