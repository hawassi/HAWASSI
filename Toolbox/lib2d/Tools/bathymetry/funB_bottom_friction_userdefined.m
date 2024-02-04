function char=funB_bottom_friction_userdefined(dom,userdata)
data.x=userdata(:,1);    
data.y=userdata(:,2);
 
mappedindex=funC_mappedmatrixindexes(dom.X,dom.Y,data.x,data.y);
char=inpolygon(dom.XX,dom.YY,dom.XX(mappedindex),dom.YY(mappedindex));
