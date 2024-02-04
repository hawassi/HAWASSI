function axes_properties(axis,fact)
lineobj = findobj(axis,'type', 'line');
set(lineobj, 'linewidth', 1.8);
%set(lineobj, 'linestyle', '--');
textobj = findobj(axis,'type', 'text'); %buat legend
%set(textobj, 'fontunits', 'points');
set(textobj, 'Fontunits','normalized','fontsize', 0.1*fact);
set(textobj, 'fontweight', 'bold');
xl = get(axis, 'xlabel');
set(xl,'fontweight', 'bold');
set(xl,'Fontunits','normalized','fontsize', 0.03*fact);
yl = get(axis, 'ylabel');
set(yl,'fontweight', 'bold');
set(yl,'Fontunits','normalized','fontsize', 0.03*fact);
zl = get(axis, 'zlabel');
set(zl,'fontweight', 'bold');
set(zl,'Fontunits','normalized','fontsize', 0.03*fact);
set(axis, 'Fontunits','normalized','fontsize', 0.03*fact);
set(axis, 'Fontweight', 'bold');
tit = get(axis, 'title');
set(tit, 'Fontunits','normalized','fontsize', 0.03*fact);
set(tit, 'Fontweight', 'bold');