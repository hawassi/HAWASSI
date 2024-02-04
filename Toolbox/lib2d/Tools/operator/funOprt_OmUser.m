function w = funOprt_OmUser(k,d,funDisp,h) 
w = zeros(length(k),length(d));
fun=['@(k,d)(',funDisp,')'];
DispUser=str2func(fun);
try
for j=1:length(d)
    w(:,j) = DispUser(k,d(j)); 
end;
catch
 [statusbarObj]=JavaFrame_handling();   
 statusbarObj.setText('Wrong format of the dispersion. Use k and d as variables');
 statusbarTxt = statusbarObj.getComponent(0);
 statusbarTxt.setForeground(java.awt.Color.red);   
end