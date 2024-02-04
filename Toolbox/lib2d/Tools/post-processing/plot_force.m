g=9.81;
bb=0.2*bath.depth;
a=input.wave.Hs/2;
if strcmp(input.body.form,'barge');
swidth=input.body.width;
fact=(g*a*bb*swidth);
Tscale=1;
elseif strcmp(input.body.form,'vertical cylinder');
r=input.body.radius;
fact=r^2*a*g;
Tscale=input.wave.Tp;
end

if any(output.bodysavevar.time(2:end)==0)
IndEnd=find(output.bodysavevar.time(2:end)==0,1)-1;
else
IndEnd=length(time)-1; 
end

betaz=interp1(output.bodysavevar.time(1:IndEnd),output.bodysavevar.betaz(1:IndEnd),output.time,'spline');
betax=interp1(output.bodysavevar.time(1:IndEnd),output.bodysavevar.betax(1:IndEnd),output.time,'spline');
forcez=-gradient(betaz,output.time(2)-output.time(1))/fact;
forcex=-gradient(betax,output.time(2)-output.time(1))/fact;

figure;
subplot(2,1,1)
plot(output.time./Tscale,betaz,'r',output.time./Tscale,betax,'b');
title('Momentum');
xlabel('time');
plot_properties;
subplot(2,1,2)
plot(output.time./Tscale,forcez,'r',output.time./Tscale,forcex,'b');
plot_properties;
title('Force')
xlabel('time');