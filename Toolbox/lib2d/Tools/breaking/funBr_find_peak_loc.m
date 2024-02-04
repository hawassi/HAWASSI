function indexMax=funBr_find_peak_loc(param,eta,dom,CrestBreakprop,I)
indxMaxPrev=CrestBreakprop.indexNow(I);
theta=CrestBreakprop.DirProp.theta(I);
kwadran=CrestBreakprop.DirProp.kwadran(I);
XX=dom.XX;YY=dom.YY;
xp=XX(indxMaxPrev);yp=YY(indxMaxPrev);
XXc=XX-xp;YYc=YY-yp;
R=param.Rsearch;%%R=1/4 wave length;
Char=zeros(size(eta));
Char(XXc.^2+YYc.^2<=R^2)=1;

if kwadran==1
    alpha=theta-80*pi/180;
    beta=(80*pi./180+theta)-pi/2;
    XXr=XXc*cos(alpha)+YYc*sin(alpha);
    YYr=YYc*cos(beta)-XXc*sin(beta);
    
    Char(YYr<0)=0;Char(XXr<0)=0;
elseif kwadran==2
    alpha=80*pi/180-theta;
    beta =(pi/2-theta)-80*pi/180;
    XXr=XXc*cos(alpha)+YYc*sin(alpha);
    YYr=YYc*cos(beta)-XXc*sin(beta);
    
    
    Char(YYr<0)=0;Char(XXr>0)=0;
elseif kwadran==3
    alpha=theta-80*pi/180;
    beta=(80*pi./180+theta)-pi/2;
    
    XXr=XXc*cos(alpha)+YYc*sin(alpha);
    YYr=YYc*cos(beta)-XXc*sin(beta);
   
    Char(YYr>0)=0; Char(XXr>0)=0;
elseif kwadran==4
    alpha=80*pi/180-theta;
    beta =(pi/2-theta)-80*pi/180;
    
    XXr=XXc*cos(alpha)+YYc*sin(alpha);
    YYr=YYc*cos(beta)-XXc*sin(beta);
    
    Char(YYr>0)=0;Char(XXr<0)=0;
end

EtaChecked=eta.*Char.*param.char;
indexx=find(EtaChecked==max(max(EtaChecked)));
if length(indexx)>1
    indexMax=min(indexx);
else
    indexMax=indexx;
end















