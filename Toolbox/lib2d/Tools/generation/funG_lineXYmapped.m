function XY=funG_lineXYmapped(dom,spatialinflux,ii)
if strcmp(spatialinflux.linetype(ii),'Straight')
    x1=cell2mat(spatialinflux.x1(ii));
    x2=cell2mat(spatialinflux.x2(ii));
    y1=cell2mat(spatialinflux.y1(ii));
    y2=cell2mat(spatialinflux.y2(ii));
    m=(y2-y1)/(x2-x1);
    if strcmpi(spatialinflux.line(ii).Orientation,'Vertical')
        Yline=dom.Y(funC_closest(dom.Y,y1):funC_closest(dom.Y,y2));
        Xline=dom.X(funC_closest(dom.X,(Yline-y1)/m+x1));
    else
        Xline=dom.X(funC_closest(dom.X,x1):funC_closest(dom.X,x2));
        Yline=dom.Y(funC_closest(dom.Y,m*(Xline-x1)+y1));
    end
else
    xc=cell2mat(spatialinflux.xc(ii));
    yc=cell2mat(spatialinflux.yc(ii));
    rad=cell2mat(spatialinflux.radius(ii));
    theta1=deg2rad(cell2mat(spatialinflux.theta1(ii)));
    theta2=deg2rad(cell2mat(spatialinflux.theta2(ii)));
    KwadTheta1=findG_kwadran_angle(rad2deg(theta1));
    
    if strcmpi(spatialinflux.line(ii).Orientation,'Horizontal')
        if KwadTheta1==1
            xmin=xc+rad*cos(theta2);
            xmax=xc+rad*cos(theta1);
            sign=1;
        elseif KwadTheta1==3
            xmin=xc+rad*cos(theta1);
            xmax=xc+rad*cos(theta2);
            sign=-1;
        end
        Xline=dom.X(funC_closest(dom.X,xmin):funC_closest(dom.X,xmax));
        Yline=dom.Y(funC_closest(dom.Y,sign*sqrt(rad^2-(Xline-xc).^2)+yc));
    else
        if KwadTheta1==2
            ymin=yc+rad*sin(theta2);
            ymax=yc+rad*sin(theta1);
            sign=-1;
        elseif KwadTheta1==4
            ymin=yc+rad*sin(theta1);
            ymax=yc+rad*sin(theta2);
            sign=1;
        end
        Yline=dom.Y(funC_closest(dom.Y,ymin):funC_closest(dom.Y,ymax));
        Xline=dom.X(funC_closest(dom.X,sign*sqrt(rad^2-(Yline-yc).^2)+xc));
    end
end
indYline=funC_closest(dom.Y,Yline);
indXline=funC_closest(dom.X,Xline);

linearindex = sub2ind(size(dom.XX),indYline ,indXline);

XY=[Xline' Yline' linearindex];

% figure;
% plot(dom.XX(linearindex),dom.YY(linearindex));