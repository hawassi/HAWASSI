function tpartitionAss=funBAssim_tpartition(tintval,tassimInit,tassimEnd,dtAssim,Npartition)

if tassimInit==tintval(1) && tassimEnd>=tintval(end)
    Ntime_part=floor(length(tintval)/Npartition);
    for I=1:Npartition-1
        tpartitionAss(I).time=tintval((I-1)*Ntime_part+1:I*Ntime_part+1);
    end
    tpartitionAss(Npartition).time=tintval((Npartition-1)*Ntime_part+1:end);
else
    if tintval(1)<tassimInit && tintval(end)>tassimEnd
        IndtAsI=funC_closest(tintval,tassimInit);
        IndtAsF=funC_closest(tintval,tassimEnd);
        tpartitionAss(1).time=tintval(1:IndtAsI);
        tintvalTemp=tintval(IndtAsI:IndtAsF);
        Ntime_part=floor(length(tintvalTemp)/(Npartition-2));
        for I=2:Npartition-1
         iter=I-1;
         tpartitionAss(I).time=tintvalTemp((iter-1)*Ntime_part+1:iter*Ntime_part+1);  
        end
        tpartitionAss(Npartition).time=tintval(IndtAsF:end); 
    elseif tintval(1)<tassimInit && tintval(end)==tassimEnd
        IndtAsI=funC_closest(tintval,tassimInit);
        IndtAsF=funC_closest(tintval,tassimEnd);
        tpartitionAss(1).time=tintval(1:IndtAsI);
        tintvalTemp=tintval(IndtAsI:end);
        Ntime_part=floor(length(tintvalTemp)/(Npartition-1));
        for I=2:Npartition-1
         iter=I-1;
         tpartitionAss(I).time=tintvalTemp((iter-1)*Ntime_part+1:iter*Ntime_part+1);  
        end
        tpartitionAss(Npartition).time=tintvalTemp((Npartition-2)*Ntime_part+1:end); 
     elseif tintval(1)==tassimInit && tintval(end)>tassimEnd
       IndtAsF=funC_closest(tintval,tassimEnd);
        tintvalTemp=tintval(1:IndtAsF);
        Ntime_part=floor(length(tintvalTemp)/(Npartition-1));
        for I=1:Npartition-2
         tpartitionAss(I).time=tintvalTemp((I-1)*Ntime_part+1:I*Ntime_part+1);  
        end
        tpartitionAss(Npartition-1).time=tintvalTemp((Npartition-1-1)*Ntime_part+1:end); 
        tpartitionAss(Npartition).time=tintval(IndtAsF:end); 
    end
    
    
end
end