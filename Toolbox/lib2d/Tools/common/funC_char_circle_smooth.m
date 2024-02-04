function char=funC_char_circle_smooth(X,Y,xc,yc,r0,smoothfact,dir)
dy=Y(2)-Y(1);dx=X(2)-X(1);
Nx=length(X);Ny=length(Y);
char=zeros(Ny,Nx);
if strcmpi(dir,'South')
    xend=cos(5*pi/180)*r0;
    x1=xc-xend;
    x2=xc+xend;
    ind1=funC_closest(X,x1);
    ind2=funC_closest(X,x2);
    
    for ii=1:Nx
        if ii<ind1 
        ycirc=yc+sqrt(r0^2-(X(ind1)-xc)^2); 
        elseif ii>ind2
        ycirc=yc+sqrt(r0^2-(X(ind2)-xc)^2); 
        else
        ycirc=yc+sqrt(r0^2-(X(ii)-xc)^2);   
        end      
        char(:,ii)=funC_cf2d(Y,ycirc,smoothfact*dy,1)';
    end
 
   
elseif strcmpi(dir,'North')
  for ii=1:Nx
        if abs(X(ii))-xc>=r0 
        char(:,ii)=funC_cf2d(Y,yc,smoothfact*dy,1)';
        else
        ycirc=yc-sqrt(r0^2-(X(ii)-xc)^2);   
        char(:,ii)=funC_cf2d(Y,ycirc,smoothfact*dy,1)';
        end      
  end
elseif strcmpi(dir,'West')
   % char=funC_cf2d(x,xc,smoothfact*dx,1);
elseif strcmpi(dir,'East')
    %char=1-funC_cf2d(x,xc-smoothfact*dx,2*smoothfact*dx,1);
end

end