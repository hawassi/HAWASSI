function axes_properties_cb(axis,fact)
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
set(gca, 'Fontunits','normalized','fontsize', 0.03*fact);
set(gca, 'Fontweight', 'bold');
tit = get(axis, 'title');
set(tit, 'Fontunits','normalized','fontsize', 0.03*fact);
set(tit, 'Fontweight', 'bold');