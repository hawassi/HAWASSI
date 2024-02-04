%plot properties     
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', 1.6);
%set(lineobj, 'linestyle', '--');
textobj = findobj('type', 'text'); %buat legend
%set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 10);
set(textobj, 'fontweight', 'bold');
xl = get(gca, 'xlabel');
set(xl,'fontweight', 'bold');
set(xl,'fontsize', 10);
yl = get(gca, 'ylabel');
set(yl,'fontweight', 'bold');
set(yl,'fontsize', 10);
zl = get(gca, 'zlabel');
set(zl,'fontweight', 'bold');
set(zl,'fontsize', 10);
set(gca, 'FontSize', 10);
set(gca, 'Fontweight', 'bold');
tit = get(gca, 'title');
set(tit, 'FontSize', 10);
set(tit, 'Fontweight', 'bold');