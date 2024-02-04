function ID=funGui_eventLoopStop(jb)
global Idstop
jbh = handle(jb,'CallbackProperties');
set(jbh,'ActionPerformedCallback',{@funGui_stopbuttonfcn})
ID=Idstop;
drawnow limitrate;%pause(1e-100);%

