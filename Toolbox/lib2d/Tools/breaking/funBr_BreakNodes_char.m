function [B, Break_nodes]=funBr_BreakNodes_char(param,U,U_star,U_F,eta,dom,CrestBreakprop,B,I)

indxMax=CrestBreakprop.indexNow(I);
theta=CrestBreakprop.DirProp.theta(I);
kwadran=CrestBreakprop.DirProp.kwadran(I);
XX=dom.XX;YY=dom.YY;
xp=XX(indxMax);yp=YY(indxMax);
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




Char(U<U_F)=0;
Char=Char.*param.char;
Char(eta<0)=0;
Char(B>0)=0;     %avoid same breaking nodes for different crest. 
Break_nodes=find(Char==1);
% figure;
% surf(XX,YY,Char,'edgecolor','none')
%   xlabel('x');ylabel('y');
% pause

for ii=1:length(Break_nodes)
    inB=Break_nodes(ii);
        if U(inB)>=2*U_star
            B(inB)=1;
        elseif  U(inB)<=U_star
            B(inB)=0;
        else
            B(inB)=(U(inB)/ U_star)-1;
        end
end
% B(U(Break_nodes)>2*U_star)=1;
% B(U(Break_nodes)<=2*U_star)=(U(U(Break_nodes)<=2*U_star)/ U_star )-1;
% B(U(Break_nodes)<=U_star)=0;
% B(U(Break_nodes)<U_F)=0;

 %[all(B(Break_nodes)==0) B(indxMax) U(indxMax) U_star]
% if  max( B(Break_nodes))==1
%    disp('1')
% elseif  max( B(Break_nodes))>1
%    disp('>>>>>>>>1')    
% end
 

% pause


