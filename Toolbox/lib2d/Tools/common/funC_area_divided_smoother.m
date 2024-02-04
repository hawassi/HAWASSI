function char = funC_area_divided_smoother(x,y,xc,yc,smoothfact,dir)
dy=y(2)-y(1);dx=x(2)-x(1);
if strcmpi(dir,'South') 
char=funC_cf2d(y,yc,2*smoothfact*dy,length(x))';
elseif strcmpi(dir,'North') 
char=1-funC_cf2d(y,yc-2*smoothfact*dy,2*smoothfact*dy,length(x))';    
elseif strcmpi(dir,'West') 
char=funC_cf2d(x,xc,2*smoothfact*dx,length(y));    
elseif strcmpi(dir,'East') 
char=1-funC_cf2d(x,xc-2*smoothfact*dx,2*smoothfact*dx,length(y));    
end

