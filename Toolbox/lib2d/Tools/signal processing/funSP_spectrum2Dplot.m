function funSP_spectrum2Dplot(spectprop,AxFig,nwpeak,Idspect,IDinput,prevSaveDir)
cla(AxFig);
halfomsig=spectprop.halfomsig;
mdir_rad=spectprop.mDir;
cos2pdf=spectprop.cos2pdf;
varDensSpect=spectprop.varDensSpect;
stdTheta=rad2deg(spectprop.StdTheta);
theta0  = linspace(-pi/2,pi/2, length(cos2pdf)); 

IndEndWW=funC_closest(halfomsig,nwpeak*0.5);
[theta,ww]=meshgrid(theta0,halfomsig(1:IndEndWW));
Spect1D=varDensSpect(1:IndEndWW);

if Idspect==3
Spect=Spect1D'*cos2pdf;
%convert to cartesian;
x=ww.*cos(theta+mdir_rad);
y=ww.*sin(theta+mdir_rad);
else
cla(AxFig);
Spect=Spect1D'*cos2pdf*(2*pi);
%convert to cartesian;
x=ww.*cos(theta+mdir_rad)./(2*pi);
y=ww.*sin(theta+mdir_rad)./(2*pi);    
end

% [x,y]=pol2cart(theta,ww);
axes(AxFig);
h = polar(x,y);
delete(h);                                                                  
hold on   
contour(AxFig,x,y,Spect,'linewidth',1);
title(AxFig,['2d Spectrum, st. dev.=',num2str(stdTheta,4),' deg'])
axes_properties(AxFig,1);
cb=colorbar('eastoutside');  
if Idspect==3
    ylabel(cb,'Var. dens. [m^2 s/rad^2]')
else
    ylabel(cb,'Var. dens. [m^2/(Hz rad)]')
end
axes_properties_cb(cb,1);

colormap(AxFig,'jet');
%axes_properties(hh,0.001);

% for saving figure
ax_old =  AxFig;
f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
ax_new = copyobj(ax_old,f_new);
set(f_new,'visible','off')
axes(ax_new)
cb=colorbar;
if Idspect==3
    ylabel(cb,'Var. dens. [m^2 s/rad^2]')
else
    ylabel(cb,'Var. dens. [m^2/(Hz rad)]')
end
axes_properties_cb(cb,1);
axes_properties(ax_new,1);
set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
set(ax_new,'outerposition',[0 0 1 1]);
saveas(ax_new,[prevSaveDir,'2DSpectrum_',num2str(IDinput),'.png']);
close(f_new)

% figure
% plot(rad2deg(theta0)-90,cos2pdf)
% 
% fig=figure;c1 = contourf(theta,ww,Spect);%,Nlevel); % calculate levels
% [z_level1,c1] = clevels(c1);                  % find contour levels
% h1            = polar(c1(1,:),c1(2,:),'b-');
