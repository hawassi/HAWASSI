% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

function callback_preview_wave(hObj,~,h)
set(h.monitorbox,'String','>>','foregroundcolor','k');

projname=get(h.project_edit_name,'string');
projdir=get(h.project_popup_projdir,'string');
prevSaveDir=[projdir,'\',projname,'\preview\'];

if ~isdir(prevSaveDir), mkdir(prevSaveDir); end


preprocdat=get(h.preview.wave_popup_var,'userdata');

%set(h.preview.wave_popup_var,'userdata',[]);%clear memory

if isempty(preprocdat)
    set(h.monitorbox,'String',['>>No data avilable, please do pre-processing'],'foregroundcolor','k');
    return;
end
delete(findall(h.panel.prev_wave,'type','axes'));

dom=preprocdat(1).dom;model=preprocdat(1).model;
par=preprocdat(1).par;influx=preprocdat(1).influx;
timeSimul=preprocdat(1).timeSimul;
g=par.g; ivp=preprocdat(1).ivp;bdyassim=preprocdat(1).bdyassim;
input=preprocdat(1).input;spatial=preprocdat(1).spatial;
if log2(dom.Nx)>9
    stepdx=2;
else
    stepdx=1;
end
if log2(dom.Ny)>9
    stepdy=2;
else
    stepdy=1;
end

Om_fun       = str2func(['funOprt_Om',model.dispersion(5:end)]);

meandepth=mean(mean(-dom.bathy.profile));

%Oprt_Up2d_fun= str2func(['Oprt_',model.phasevel]);
kx=linspace(0,2*max(dom.kx)/dom.cutfracx,100);%% max(dom.kx)is kmax/2
ky=linspace(0,2*max(dom.ky)/dom.cutfracy,100);
[Kx,Ky]=meshgrid(kx,ky);
KK=sqrt(Kx.^2+Ky.^2);
kk=linspace(0,max(max(KK)),100);

[xx,yy]=meshgrid(dom.X(1:stepdx:end),dom.Y(1:stepdy:end));


IDinput=get(hObj,'value');
IDspect=get(h.preview.wave_popup_spect,'value');

axes1= axes('parent', h.panel.prev_wave, 'Position', [0.07 0.12 0.38 0.65]);
axes2= axes('parent', h.panel.prev_wave, 'Position', [0.57 0.12 0.33 0.65]);

cla(axes1);cla(axes2);
if strcmpi(input.wave.option,'Yes')
    Nw=input.wave.N;
    
    if Nw>1
        Nend=Nw+1;
        flagcombineSpect=1;
    else
        Nend=Nw;
        flagcombineSpect=0;
    end
    
    if IDinput<=Nend
        set(h.preview.wave_popup_spect,'enable','on');
        if flagcombineSpect==0 || (flagcombineSpect==1&&IDinput<=Nw)
            timesig=influx.wave(IDinput).time;
            
            if floor(log2(timesig))>10
                stept=2+floor(log2(timesig))-10;
            else
                stept=1;
            end
            ramp2d=funG_ramp2d(input.wave,dom,timesig,...
                spatial.influx,influx.par.T_p(IDinput),IDinput);
            EtaIDinput=ramp2d.*influx.wave(IDinput).eta;
            inflXY=spatial.influx.line(IDinput).xy;
            if strcmp(spatial.influx.line(IDinput).Orientation,'Vertical')
                Ind1=funC_closest(dom.Y,min(inflXY(:,2)));
                Ind2=funC_closest(dom.Y,max(inflXY(:,2)));
                
                [Y,T]=meshgrid(dom.Y(Ind1:stepdy:Ind2),timesig(1:stept:end));
                
                surf(axes1,T,Y,EtaIDinput(1:stept:end,Ind1:stepdy:Ind2),'edgecolor','none');
                xlim(axes1,[min(min(T)) max(max(T))]);
                ylim(axes1,[dom.Y(Ind1) dom.Y(Ind2)]);
                xlabel(axes1,'time[s]');ylabel(axes1,'y[m]');
            else
                Ind1=funC_closest(dom.X,min(inflXY(:,1)));
                Ind2=funC_closest(dom.X,max(inflXY(:,1)));
                [T,X]=meshgrid(timesig(1:stept:end),dom.X(Ind1:stepdx:Ind2));
                surf(axes1,X,T,EtaIDinput(1:stept:end,Ind1:stepdx:Ind2)','edgecolor','none');
                xlim(axes1,[dom.X(Ind1) dom.X(Ind2)]);
                ylim(axes1,[min(min(T)) max(max(T))]);
                xlabel(axes1,'x[m]');ylabel(axes1,'time[s]');
            end
            title(axes1,['Wave-influx #', num2str(IDinput)]);
            axes_properties(axes1,1.5);
            view(axes1,2)
            axes(axes1)
            colormap(axes1,'jet');
            cb=colorbar;
            ylabel(cb,'elevation [m]')
            %            cbfreeze;cblabel('elevation [m]')
            %            freezeColors;
            axes_properties_cb(cb,1.5);
            
            % for saving figure
            ax_old = axes1;
            f_new = figure('position',[1.1 1.1 0.7 0.7]);
            set(f_new,'Renderer','zbuffer');
            ax_new = copyobj(ax_old,f_new);
            set(f_new,'visible','off')
            colormap(ax_new,'jet');
            axes(ax_new)
            cb=colorbar;
            ylabel(cb,'elevation [m]','fontweight', 'bold')
            axes_properties(ax_new,1);
            axes_properties_cb(cb,1.5);
            set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
            set(ax_new,'outerposition',[0 0 1 1]);
            saveas(ax_new,[prevSaveDir,'influxsignal_',num2str(IDinput),'.png']);
            close(f_new)
        end
        omMax= Om_fun(2*max(kk),meandepth,g,'');
        
        if flagcombineSpect==0 || (flagcombineSpect==1&&IDinput<=Nw)
            
            if   input.wave.type(IDinput).flag  ==  4
                if IDspect>2, IDspect=1;
                    set(h.monitorbox,'String',['>>2D spectrum is not available'],'foregroundcolor','k');
                    set(h.preview.wave_popup_spect,'value',1);
                end
            else
                if length(influx.Spect(IDinput).prop.cos2pdf)==1
                    IDspect=1;
                    set(h.monitorbox,'String',['>>2D spectrum is not available'],'foregroundcolor','k');
                    set(h.preview.wave_popup_spect,'value',1);
                end
            end
            
            if IDspect==3 || IDspect==4
                funSP_spectrum2Dplot(influx.Spect(IDinput).prop,axes2,omMax,IDspect,IDinput,prevSaveDir);
                axes_properties(axes2,1.5);
           end
            inflXY=spatial.influx.line(IDinput).xy;
            
            if IDspect==1 || IDspect==2
                omsig=funC_freqspace(timesig);
                IndEndW=funC_closest(omsig,omMax);
                if strcmp(spatial.influx.line(IDinput).Orientation,'Vertical')
                    Ymid=(max(inflXY(:,2))+min(inflXY(:,2)))/2;
                    IndS=funC_closest(dom.Y,Ymid);
                else
                    Xmid=(max(inflXY(:,1))+min(inflXY(:,1)))/2;
                    IndS=funC_closest(dom.X,Xmid);
                end
                
                [OmS,S_hat]=funSP_variance_density_spectrum1d(timesig,EtaIDinput(:,IndS));
              
                 if strcmpi(dom.wall.option,'yes')
                    if dom.wall.ReflCoef.FlagFreqDep==1
                     yyaxis(axes2,'left');
                    end
                 end
                if IDspect==1
                    plot(axes2,OmS,S_hat);
                    try
                        xlim(axes2,[OmS(1) OmS(IndEndW)]);
                    catch
                        
                    end
                    xlabel(axes2,'\omega [rad/s]');ylabel(axes2,'Variance density [m^2 s/rad]')
                else
                    plot(axes2,OmS./(2*pi),2*pi*S_hat);
                    try
                        xlim(axes2,[OmS(1)./(2*pi) OmS(IndEndW)./(2*pi)]);
                    catch
                        
                    end
                    xlabel(axes2,'f [Hz]');ylabel(axes2,'Variance density [m^2/Hz]')
                end
                
                if strcmpi(dom.wall.option,'yes')
                    if dom.wall.ReflCoef.FlagFreqDep==1
                     for i=1:dom.wall.N
                     if dom.wall.ReflCoef.FreqDep(i).flag==1
                         yyaxis(axes2,'right');
                         hold(axes2,'on');
                         if IDspect==1
                         plot(axes2,OmS(1:IndEndW),dom.wall.ReflCoef.fun_Om(i).coef(1:IndEndW));
                         else
                         plot(axes2,OmS(1:IndEndW)/(2*pi),dom.wall.ReflCoef.fun_Om(i).coef(1:IndEndW));
                         
                         end
                         ylabel(axes2,'Wall Refl. Coef')
                         hold(axes2,'off');
                     end
                     end
                    end
                end
                
                title(axes2,'1d Spectrum');
                axes_properties(axes2,1.5);
                
                ax_old = axes2;
                f_new = figure('position',[1.1 1.1 0.7 0.7]);
                set(f_new,'Renderer','zbuffer');
                
                try
                % for saving figure
                ax_new = copyobj(ax_old,f_new);
                set(f_new,'visible','off')
                axes(ax_new)
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
                saveas(ax_new,[prevSaveDir,'1DSpectrum_',num2str(IDinput),'.png']);
                catch
                    
                end
                 close(f_new);
            end
        else
            delete(axes1);
            funSP_combineSpectrumPlot(axes2,influx,timeSimul,IDspect,input,...
                spatial,dom,omMax,prevSaveDir);
        end
        
        
    else
        set(h.preview.wave_popup_spect,'enable','off');
        if input.bdyassim.option==0
            surf(axes1,xx,yy,ivp.eta(1:stepdy:end,1:stepdx:end),'edgecolor','none');
        else
            Nt=input.bdyassim.assimdata(1,3);
            Nx=input.bdyassim.assimdata(2,3);
            Ny=input.bdyassim.assimdata(3,3);
            bdyeta=input.bdyassim.assimdata(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
            Xdat=linspace(input.bdyassim.assimdata(2,1),input.bdyassim.assimdata(2,2),Nx);
            Ydat=linspace(input.bdyassim.assimdata(3,1),input.bdyassim.assimdata(3,2),Ny);
            EtaBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyeta,dom);
            meandepthAssim=mean(mean(-dom.bathy.profile));
            surf(axes1,xx,yy,EtaBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
        end
        xlim(axes1,[min(min(xx)) max(max(xx))]);
        ylim(axes1,[min(min(yy)) max(max(yy))]);
        xlabel(axes1,'x [m]');ylabel(axes1,'y[m]');
        title(axes1,'Initial condition of wave elevation');
        axes_properties(axes1,1.5);
        view(axes1,2)
        axes(axes1)
        colormap(axes1,'jet');
        cb=colorbar;
        ylabel(cb,'elevation [m]');
        axes_properties_cb(cb,1.5);
        %     freezeColors;cbfreeze;
        if input.bdyassim.option==0
            if model.phiForm ==1
                surf(axes2,xx,yy,ivp.phi(1:stepdy:end,1:stepdx:end),'edgecolor','none');
            else
                ivpUV=sqrt(ivp.u(1:stepdy:end,1:stepdx:end).^2+ivp.v(1:stepdy:end,1:stepdx:end).^2);
                surf(axes2,xx,yy,ivpUV,'edgecolor','none');
            end
        else
            if model.phiForm ==1
                if input.bdyassim.cb_phi==1
                    bdyphi=input.bdyassim.assimdata_phi(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                    PhiBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyphi,dom);
                else
                    propdir=bdyassim.propdir;
                    PhiBdy=funOprt_PhifromEta(g,dom,model,meandepthAssim,EtaBdy,propdir);
                end
                PhiBdy=PhiBdy.*bdyassim.charupdate.*dom.cfSA;
                
                surf(axes2,xx,yy,PhiBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
            else
                if input.bdyassim.cb_vel==1
                    bdyu=input.bdyassim.assimdata_u(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                    UBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyu,dom);
                    bdyv=input.bdyassim.assimdata_u(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                    VBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyv,dom);
                else
                    propdir=bdyassim.propdir;
                    PhiBdy=funOprt_PhifromEta(g,dom,model,meandepthAssim,EtaBdy,propdir);
                    gradphi=funOprt_grad2d(dom.Kx,dom.Ky,fft2(PhiBdy));
                    UBdy =gradphi.x; VBdy =gradphi.y;
                end
                UBdy=UBdy.*bdyassim.charupdate.*dom.cfSA;
                VBdy=VBdy.*bdyassim.charupdate.*dom.cfSA;
                UVBdy=sqrt(UBdy.^2+VBdy.^2);
                
                surf(axes2,xx,yy,UVBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
            end
        end
        xlim(axes2,[min(min(xx)) max(max(xx))]);
        ylim(axes2,[min(min(yy)) max(max(yy))]);
        xlabel(axes2,'x [m]');ylabel(axes2,'y[m]');
        if model.phiForm ==1
            title(axes2,'Initial condition of wave potential');
        else
            title(axes2,'Initial condition of absolute wave velocity');
        end
        
        axes_properties(axes2,1.5);
        view(axes2,2)
        axes(axes2)
        colormap(axes2,'jet');
        cb=colorbar;
        if  model.phiForm ==1
            ylabel(cb, 'potential [m^2/s]');
        else
            ylabel(cb, '|velocity| [m/s]');
        end
        axes_properties_cb(cb,1.5);
        %       freezeColors;cbfreeze;
        % for saving figure
        ax_old = axes1;
        f_new = figure('position',[1.1 1.1 0.7 0.7]);
        set(f_new,'Renderer','zbuffer');
        ax_new = copyobj(ax_old,f_new);
        title(ax_new,'Initial condition of wave elevation');
        xlabel(ax_new,'x [m]');ylabel(ax_new,'y[m]');
        set(f_new,'visible','off')
        colormap(ax_new,'jet');
        axes(ax_new)
        y=colorbar;
        ylabel(y,'elevation [m]','fontweight', 'bold')
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,[prevSaveDir,'IC_eta.png']);
        close(f_new)
        
        % for saving figure
        ax_old = axes2;
        f_new = figure('position',[1.1 1.1 0.7 0.7]);
        set(f_new,'Renderer','zbuffer');
        ax_new = copyobj(ax_old,f_new);
         if model.phiForm ==1
        title(ax_new,'Initial condition of wave potential');
         else
           title(ax_new,'Initial condition of absolute wave velocity');
         end
        xlabel(ax_new,'x [m]');ylabel(ax_new,'y[m]');
        set(f_new,'visible','off')
        colormap(ax_new,'jet');
        axes(ax_new)
        cb=colorbar;
        if  model.phiForm ==1
            ylabel(cb, 'potential [m^2/s]');
        else
            ylabel(cb, '|velocity| [m/s]');
        end
        axes_properties_cb(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
         if model.phiForm ==1
        saveas(ax_new,[prevSaveDir,'IC_phi.png']);
         else
         saveas(ax_new,[prevSaveDir,'IC_absvel.png']);     
         end
        close(f_new)
        
    end
    
else
    set(h.preview.wave_popup_spect,'enable','off');
    
    if input.bdyassim.option==0
        surf(axes1,xx,yy,ivp.eta(1:stepdy:end,1:stepdx:end),'edgecolor','none');
    else
        Nt=input.bdyassim.assimdata(1,3);
        Nx=input.bdyassim.assimdata(2,3);
        Ny=input.bdyassim.assimdata(3,3);
        Xdat=linspace(input.bdyassim.assimdata(2,1),input.bdyassim.assimdata(2,2),Nx);
        Ydat=linspace(input.bdyassim.assimdata(3,1),input.bdyassim.assimdata(3,2),Ny);
        bdyeta=input.bdyassim.assimdata(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
        EtaBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyeta,dom);
        meandepthAssim=mean(mean(-dom.bathy.profile));
        EtaBdy=EtaBdy.*bdyassim.charupdate.*dom.cfSA;
        surf(axes1,xx,yy,EtaBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
    end
    %axis(axes1,'equal')
    xlabel(axes1,'x [m]');ylabel(axes1,'y[m]');
    xlim(axes1,[min(min(xx)) max(max(xx))]);
    ylim(axes1,[min(min(yy)) max(max(yy))]);
    title(axes1,'Initial condition of wave elevation');
    axes_properties(axes1,1.5);
    view(axes1,2)
    axes(axes1)
    colormap(axes1,'jet');cb=colorbar;
    ylabel(cb, 'elevation [m]');
    axes_properties_cb(cb,1.5)
    
    
    
   
        
    %   freezeColors;cbfreeze;
    if input.bdyassim.option==0
        if model.phiForm ==1
        surf(axes2,xx,yy,ivp.phi(1:stepdy:end,1:stepdx:end),'edgecolor','none');
        else
        ivpUV=sqrt(ivp.u(1:stepdy:end,1:stepdx:end).^2+ivp.v(1:stepdy:end,1:stepdx:end).^2);
        surf(axes2,xx,yy,ivpUV,'edgecolor','none');
        end
    else
        
        if model.phiForm ==1
            if input.bdyassim.cb_phi==1
                bdyphi=input.bdyassim.assimdata_phi(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                PhiBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyphi,dom);
            else
                propdir=bdyassim.propdir;
                PhiBdy=funOprt_PhifromEta(g,dom,model,meandepthAssim,EtaBdy,propdir);
            end
            PhiBdy=PhiBdy.*bdyassim.charupdate.*dom.cfSA;
            
             surf(axes2,xx,yy,PhiBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
        else
            if input.bdyassim.cb_vel==1
                bdyu=input.bdyassim.assimdata_u(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                UBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyu,dom);
                bdyv=input.bdyassim.assimdata_u(4+(floor(Nt/2)-1)*Ny:3+floor(Nt/2)*Ny,:);
                VBdy=funBAssim_interpEtatoSimDom(Xdat,Ydat,bdyv,dom);
            else
                propdir=bdyassim.propdir;
                PhiBdy=funOprt_PhifromEta(g,dom,model,meandepthAssim,EtaBdy,propdir);
                gradphi=funOprt_grad2d(dom.Kx,dom.Ky,fft2(PhiBdy));
                UBdy =gradphi.x; VBdy =gradphi.y;
            end
            UBdy=UBdy.*bdyassim.charupdate.*dom.cfSA;
            VBdy=VBdy.*bdyassim.charupdate.*dom.cfSA;
            UVBdy=sqrt(UBdy.^2+VBdy.^2);
            
             surf(axes2,xx,yy,UVBdy(1:stepdy:end,1:stepdx:end),'edgecolor','none');
        end
       
    end
    xlabel(axes2,'x [m]');ylabel(axes2,'y[m]');
    % axis(axes2,'equal')
    xlim(axes2,[min(min(xx)) max(max(xx))]);
    ylim(axes2,[min(min(yy)) max(max(yy))]);
     if model.phiForm ==1
    title(axes2,'Initial condition of wave potential');
     else
    title(axes2,'Initial condition of absolute wave velocity');
     end
    axes_properties(axes2,1.5);
    view(axes2,2)
    axes(axes2)
    colormap(axes2,'jet'); 
    cb=colorbar;
    if  model.phiForm ==1
         ylabel(cb, 'potential [m^2/s]');
    else
         ylabel(cb, '|velocity| [m/s]'); 
    end
    axes_properties_cb(cb,1.5)
    
    % for saving figure
    ax_old = axes1;
    f_new = figure('position',[1.1 1.1 0.7 0.7]);
    set(f_new,'Renderer','zbuffer');
    ax_new = copyobj(ax_old,f_new);
    title(ax_new,'Initial condition of wave elevation');
    xlabel(ax_new,'x [m]');ylabel(ax_new,'y[m]');
    set(f_new,'visible','off')
    colormap(ax_new,'jet');
    axes(ax_new)
    y=colorbar;
    ylabel(y,'elevation [m]','fontweight', 'bold')
    axes_properties(ax_new,1);
    axes_properties_cb(y,1);
    set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
    set(ax_new,'outerposition',[0 0 1 1]);
    saveas(ax_new,[prevSaveDir,'IC_eta.png']);
    close(f_new)
    
    % for saving figure
    ax_old = axes2;
    f_new = figure('position',[1.1 1.1 0.7 0.7]);
    set(f_new,'Renderer','zbuffer');
    ax_new = copyobj(ax_old,f_new);
    if model.phiForm ==1
    title(ax_new,'Initial condition of wave potential');
    else
    title(ax_new,'Initial condition of absolute wave velocity');
    end
    xlabel(ax_new,'x [m]');ylabel(ax_new,'y[m]');
    set(f_new,'visible','off')
    colormap(ax_new,'jet');
    axes(ax_new)
    y=colorbar;
    if model.phiForm ==1
    ylabel(y,'potential [m^2/s]','fontweight', 'bold')
    else
    ylabel(y,'velocity [m/s]','fontweight', 'bold')
    end
    axes_properties(ax_new,1);
    axes_properties_cb(y,1);
    set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
    set(ax_new,'outerposition',[0 0 1 1]);
    if model.phiForm ==1
    saveas(ax_new,[prevSaveDir,'IC_phi.png']);
    else
    saveas(ax_new,[prevSaveDir,'IC_absvel.png']);
    end
    close(f_new)
end


end