function ramp2d=funG_ramp2d(inputwave,dom,timesig,Spatinflux,Tp,I)
ramp=inputwave.ramp;
influxlineOr=Spatinflux.line(I).Orientation;
tapered=inputwave.tapered;
XYinfl=Spatinflux.line(I).xy;
if  ramp.check==1
    dt=timesig(2)-timesig(1);
    Lramp.b    = (round(ramp.lengthfact*Tp)/dt+1)*dt;
    Lramp.t    = Lramp.b;
    
    if strcmp(influxlineOr,'Vertical') 
        ymin=dom.Y(funC_closest(dom.Y,min(XYinfl(:,2))));
        ymax=dom.Y(funC_closest(dom.Y,max(XYinfl(:,2))));
        
        if tapered.check==1
            lengthInf=ymax-ymin;
            Lramp.r=(round(tapered.length*lengthInf/dom.dy)+1)*dom.dy;
            Lramp.l=Lramp.r;
        else
            Lramp.r    = 0; Lramp.l   = 0;
        end
        ramp2d= funC_tapered2d(dom.Y,timesig,ymin,ymax,timesig(1),timesig(end),Lramp);
    else
       
%         Lramp.b    = ramp.lengthfact*Tp;
%         Lramp.t    = ramp.lengthfact*Tp;
        xmin=dom.X(funC_closest(dom.X,min(XYinfl(:,1))));
        xmax=dom.X(funC_closest(dom.X,max(XYinfl(:,1))));
            
        if tapered.check==1
            lengthInf=xmax-xmin;
            Lramp.r=(round(tapered.length*lengthInf/dom.dx)+1)*dom.dx;
            Lramp.l=Lramp.r;
        else
            Lramp.r    = 0; Lramp.l    = 0;
        end
        ramp2d= funC_tapered2d(dom.X,timesig,xmin,xmax,timesig(1),timesig(end),Lramp);
    end
else
    ramp2d=1;
end

% [X,T]=meshgrid(dom.X,timesig);
% figure;
% surf(X,T,ramp2d,'edgecolor','none')
% view(2)