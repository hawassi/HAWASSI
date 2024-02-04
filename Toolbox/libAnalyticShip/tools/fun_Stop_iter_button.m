function fun_Stop_iter_button(jbStop)
jbh = handle(jbStop,'CallbackProperties');
set(jbh,'ActionPerformedCallback',{@stopbuttonfcn})
try
    drawnow limitrate;
catch
    pause(1e-100);%<2014b
end