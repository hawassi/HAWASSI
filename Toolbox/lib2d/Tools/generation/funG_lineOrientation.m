function Orientation=funG_lineOrientation(spatialinflux,ii)

if strcmp(spatialinflux.linetype(ii),'Straight')
x1=cell2mat(spatialinflux.x1(ii));
x2=cell2mat(spatialinflux.x2(ii));
y1=cell2mat(spatialinflux.y1(ii));
y2=cell2mat(spatialinflux.y2(ii));
if x2-x1>=y2-y1
   Orientation={'Horizontal'}; 
else
   Orientation={'Vertical'}; 
end

else
theta1=cell2mat(spatialinflux.theta1(ii));
if (theta1>=0 && theta1<90) || (theta1>=180 && theta1<270)
 Orientation={'Horizontal'};    
else
 Orientation={'Vertical'};        
end

end
end