function Om2d = funOprt_Om2dUser(Kx,Ky,depth,funDisp,h) 
K=sqrt(Kx.^2+Ky.^2);
fun=['@(k,d)(',funDisp,')'];
DispUser=str2func(fun);
try
    Om2d=DispUser(K,depth);
catch
 [statusbarObj]=JavaFrame_handling();   
 statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
 statusbarTxt = statusbarObj.getComponent(0);
 statusbarTxt.setForeground(java.awt.Color.red);   
end