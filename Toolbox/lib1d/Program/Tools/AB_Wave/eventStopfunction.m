function [value,isterminal,direction] = eventStopfunction(time,z,jb)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.

global Idstop

jbh = handle(jb,'CallbackProperties');
set(jbh,'ActionPerformedCallback',{@stopbuttonfcn})

try
drawnow limitrate; 
catch
pause(1e-100);%<2014b
end

if Idstop==1
    value      = 0;  % when value = 0, an event is triggered
    isterminal = 1; % terminate after the first event
    direction  = 0;  % get all the zeros
else
    value      = 1;
    isterminal = 0; % terminate after the first event
    direction  = 0;
end
