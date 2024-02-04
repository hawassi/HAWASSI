function InflDir=funG_findInfluxDirUserDefined(spatialinflux,ii,dom)

XYinfl=spatialinflux.line(ii).xy;
if strcmpi(spatialinflux.line(ii).Orientation,'Vertical')
    Xinfmin=min(XYinfl(:,1));Xinfmax=max(XYinfl(:,1));
    Xmid=min(XYinfl(:,1))+(Xinfmax-Xinfmin)/2;
    dist1= Xmid-dom.X(1);
    dist2= dom.X(end)-Xmid;
    if dist1<dist2
        InflDir={'UniR'};
    else
        InflDir={'UniL'};
    end
else
    Yinfmin=min(XYinfl(:,2));Yinfmax=max(XYinfl(:,2));
    Ymid=min(XYinfl(:,2))+(Yinfmax-Yinfmin)/2;
    dist1= Ymid-dom.Y(1);
    dist2= dom.Y(end)-Ymid;
    if dist1<dist2
        InflDir={'UniT'};
    else
        InflDir={'UniB'};
    end
end
