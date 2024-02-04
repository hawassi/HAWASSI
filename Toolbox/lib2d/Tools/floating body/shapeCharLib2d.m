function [chi]=shapeCharLib2d(dom,inputbody,I)
xyz0=inputbody.xyz0(I,:);
form=inputbody.form{I};


switch form
    case 'barge'
         LL =inputbody.length(I);
         WW =inputbody.width(I);
         yb=xyz0(2)-WW./2;yt=xyz0(2)+WW./2;         
         xl=xyz0(1)-LL./2;xr=xyz0(1)+LL./2;
         bb=yb-dom.Y(1);
         tt=dom.Y(end)-yt;
         HeavY=funC_Heav2d(dom.Y,bb,tt,dom.Nx);
         ll=xl-dom.X(1);
         rr=dom.X(end)-xr;
         HeavX=funC_Heav2d(dom.X,ll,rr,dom.Ny)';
         chi=(HeavX.*HeavY); 
    case 'vertical cylinder'
         RR =inputbody.radius(I);
         xc=xyz0(1);yc=xyz0(2);
         XX=dom.XX;YY=dom.YY;
         chi=zeros(size(dom.XX));
         chi((XX-xc).^2+(YY-yc).^2<=RR.^2)=1;
end
    
    