function ID=eventLoopStop(jb)
global Idstop
jbh = handle(jb,'CallbackProperties');
set(jbh,'ActionPerformedCallback',{@stopbuttonfcn})
ID=Idstop;

try
drawnow limitrate;%
catch
pause(1e-100);%
end
