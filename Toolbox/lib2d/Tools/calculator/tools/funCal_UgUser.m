function u = UgUser(k,d,funUg)
J=length(d);
u = zeros(length(k),J);
fun=['@(k,d)(',funUg,')'];
VelUser=str2func(fun);
try
for j=1:J
    for i=1:length(k)
        if k(i)==0
            u(i,j)=sqrt(9.81*d(j));
        else
            u(i,j) = VelUser(k(i),d(j));
        end
    end
end
catch
 [statusbarObj]=JavaFrame_handling();   
 statusbarObj.setText('Wrong format of the group velocity. Use k and d as variables');
 statusbarTxt = statusbarObj.getComponent(0);
 statusbarTxt.setForeground(java.awt.Color.red);   
end