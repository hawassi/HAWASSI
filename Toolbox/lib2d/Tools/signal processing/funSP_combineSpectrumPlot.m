function funSP_combineSpectrumPlot(axes1,influx,timeSimul,IDspect,input,spatial,dom,omMax,prevSaveDir)

Nw=input.wave.N;
omTsimul=funC_freqspace(timeSimul.interval);
omS=omTsimul(1:floor(end/2));
IndEndW0=funC_closest(omS,omMax*0.5);
Ss_hat=0;
SpectS=0;

if IDspect==3 || IDspect==4
    
    theta0  = linspace(0,2*pi, 2*length(omS));
    [theta,ww]=meshgrid(theta0, omS(1:IndEndW0));
    if IDspect==3  
        x=ww.*cos(theta);
        y=ww.*sin(theta);
    else  
        x=ww.*cos(theta)./(2*pi);
        y=ww.*sin(theta)./(2*pi);
    end
%     assignin('base','x',x)
%     assignin('base','y',y)
    
    xx=reshape(x,[],1);
    yy=reshape(y,[],1);
    x1=xx;
    y1=yy;
end

for ii=1:Nw
    if IDspect==1||IDspect==2
        timesig=influx.wave(ii).time;
        omsig=funC_freqspace(timesig);
        
        ramp2d=funG_ramp2d(input.wave,dom,timesig,...
            spatial.influx,influx.par.T_p(ii),ii);
        Etaii=ramp2d.*influx.wave(ii).eta;
        
        inflXY=spatial.influx.line(ii).xy;
        if strcmpi(spatial.influx.line(ii).Orientation,'Vertical')
            Ymid=(max(inflXY(:,2))+min(inflXY(:,2)))/2;
            IndS=funC_closest(dom.Y,Ymid);
        else
            Xmid=(max(inflXY(:,1))+min(inflXY(:,1)))/2;
            IndS=funC_closest(dom.X,Xmid);
        end
        varSig=var(Etaii(:,IndS));
        S_hat=abs(fft(Etaii(:,IndS))).^2;
        S_hat=S_hat(1:floor(end/2));
        ww=omsig(1:floor(end/2));
        var0=trapz(ww,S_hat);
        S_hat=S_hat*varSig/var0;
        
        Ss_hat=Ss_hat+abs(interp1(ww,S_hat,omS,'spline'));
        
    else
        spectprop=influx.Spect(ii).prop;
        halfomsig=spectprop.halfomsig;
        mdir_rad=spectprop.mDir;
        cos2pdf=spectprop.cos2pdf;
        varDensSpect=spectprop.varDensSpect;
        stdTheta=rad2deg(spectprop.StdTheta);
        theta0  = linspace(-pi/2,pi/2, length(cos2pdf));
        
        IndEndWW=funC_closest(halfomsig,omMax*0.6);
        [theta,ww]=meshgrid(theta0,halfomsig(1:IndEndWW));     
        
        Spect1D=varDensSpect(1:IndEndWW);
        
        if IDspect==3
            Spect=Spect1D'*cos2pdf;
            %convert to cartesian;
            x=ww.*cos(theta+mdir_rad);
            y=ww.*sin(theta+mdir_rad);
        else
            Spect=Spect1D'*cos2pdf*(2*pi);
            %convert to cartesian;
            x=ww.*cos(theta+mdir_rad)./(2*pi);
            y=ww.*sin(theta+mdir_rad)./(2*pi);
        end
        
        xx=reshape(x,[],1);
        yy=reshape(y,[],1);
        SS=reshape(Spect,[],1);
        
            Fs=scatteredInterpolant(xx,yy,SS,'nearest');
            SpectS=SpectS+reshape(Fs(x1,y1),[],length(Spect(1,:)));
    end
end

if IDspect==1||IDspect==2
    if IDspect==1
        plot(axes1,omS(1:IndEndW0),Ss_hat(1:IndEndW0));
        xlim(axes1,[omS(1) omS(IndEndW0)]);
        xlabel(axes1,'\omega [rad/s]');ylabel(axes1,'Variance density [m^2 s/rad]')
    else
        plot(axes1,omS(1:IndEndW0)./(2*pi),2*pi*Ss_hat(1:IndEndW0));
        xlim(axes1,[omS(1)./(2*pi) omS(IndEndW0)./(2*pi)]);
        xlabel(axes1,'f [Hz]');ylabel(axes1,'Variance density [m^2/Hz]')
    end
    title(axes1,'Combined 1d Spectrum')
    axes_properties(axes1,1.5);
    
    
    % for saving figure
    ax_old = axes1;
    f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
    set(f_new,'Renderer','zbuffer');
    ax_new = copyobj(ax_old,f_new);
    set(f_new,'visible','off')
    colormap(ax_new,'jet');
    axes(ax_new)
    title(ax_new,'Combined 1d Spectrum')
    if IDspect==1
        xlabel(ax_new,'\omega [rad/s]','fontweight', 'bold')
        ylabel(ax_new,'Variance density [m^2 s/rad]','fontweight', 'bold')
    else
        xlabel(ax_new,'f [Hz]','fontweight', 'bold')
        ylabel(ax_new,'Variance density [m^2/Hz]','fontweight', 'bold')
    end
    axes_properties(ax_new,1);
    set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
    set(ax_new,'outerposition',[0 0 1 1]);
    saveas(ax_new,[prevSaveDir,'1DSpectrum_combined.png']);
    close(f_new)
elseif IDspect==3 || IDspect==4
    axes(axes1);
    x1=reshape(x1,[],length(SpectS(1,:)));
    y1=reshape(y1,[],length(SpectS(1,:)));
    h = polar(x1,y1);
    delete(h);
    hold on
   
    contour(axes1,x1,y1,SpectS,'linewidth',1);
    title(axes1,['Combined 2d Spectrum '])
    axes_properties(axes1,1);
    cb=colorbar('eastoutside');
    if IDspect==3
          ylabel(cb,'Var. dens. [m^2 s/rad^2]')
    else
          ylabel(cb,'Var. dens. [m^2/(Hz rad)]')
    end
    axes_properties_cb(cb,1);
    colormap(axes1,'jet');
%    axes_properties_cb(hh,0.001);
    
    
    % for saving figure
    ax_old = axes1;
    f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
    ax_new = copyobj(ax_old,f_new);
    set(f_new,'visible','off')
    colormap(ax_new,'jet');
    axes(ax_new)
    cb=colorbar('eastoutside');
   if IDspect==3
          ylabel(cb,'Var. dens. [m^2 s/rad^2]')
    else
          ylabel(cb,'Var. dens. [m^2/(Hz rad)]')
    end
    axes_properties_cb(cb,1);
    axes_properties(ax_new,1);
    set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
    set(ax_new,'outerposition',[0 0 1 1]);
    saveas(ax_new,[prevSaveDir,'CombinedSpectrum.png'],'png');
    close(f_new);
    
end