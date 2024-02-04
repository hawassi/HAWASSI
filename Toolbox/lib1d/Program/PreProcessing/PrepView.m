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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Preview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutfrac = spatial.cutfrac;
kcut_Simul=k(floor(length(k)/cutfrac));
OmKcut_simul=OmExact(kcut_Simul,depthinflux);

nFreqP=model.nonlinear;      %
OmExct=str2func('OmExact');
influx.k_cut=invOm(nFreqP*nupeak,depthinflux,OmExct,omAdd);

cfSA  =par.cfSA;
bathy =par.bathy; depth=par.depth;

maxbath =abs(max(-bathy));

Xinflux =par.Xinflux;
dk      =par.dk;

if input.type~=1
    hfig=figure('Name','Prepared', 'Position',[10,10,2,600]);
    set(hfig,'Renderer','zbuffer'); %due to graphics driver
    if IVP.type==5
        hf=subplot(6,1,1);
    elseif IVP.type~=5 && IVP.type~=1
        hf=subplot(5,1,1);
    else
        hf=subplot(4,1,1);
    end
    if model.nonlinear==1
        if strcmp(bath.wall.check, 'Yes')
            
            if bath.friction.check==1
                plot(x,spatgamX,'b',x,cfSA,'r--',x,bathy/maxbath,'k',...
                    Xinflux,0,'*b',x(par.friction.index),par.friction.Cf,'+y',Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
                leg=legend('influx function', 'simulation area','bathy/maxbathy','influx position','friction coef.','wall position');
            else
                plot(x,spatgamX,'b',x,cfSA,'r--',x,bathy/maxbath,'k',...
                    Xinflux,0,'*b',Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
                leg=legend('influx function', 'simulation area','bathy/maxbathy','influx position','wall position');
            end
        else
            if bath.friction.check==1
                plot(x,spatgamX,'b',x,cfSA,'r--',x,bathy/maxbath,'k',...
                    Xinflux,0,'b*',x(par.friction.index),par.friction.Cf,'+y', x,0*x,'b:');
                leg=legend('influx function', 'simulation area','bathy/maxbathy','influx position','friction coef.');
            else
                plot(x,spatgamX,'b',x,cfSA,'r--',x,bathy/maxbath,'k',...
                    Xinflux,0,'b*', x,0*x,'b:');
                leg=legend('influx function', 'simulation area','bathy/maxbathy','influx position');
            end
        end
    else
        if strcmp(bath.wall.check, 'Yes')
            if bath.friction.check==1
                plot(x,spatgamX,'b',x,cfSA,'r--',x,ChiAdj,'g',x,bathy/maxbath,'k',...
                    Xinflux,0,'*b',x(par.friction.index),par.friction.Cf,'+y', Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
                leg=legend('influx function', 'simulation area','adjustment zone','bathy/maxbathy','influx position','friction coef.','wall position');
            else
                plot(x,spatgamX,'b',x,cfSA,'r--',x,ChiAdj,'g',x,bathy/maxbath,'k',...
                    Xinflux,0,'*b', Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
                leg=legend('influx function', 'simulation area','adjustment zone','bathy/maxbathy','influx position','wall position');
            end
            
        else
            if bath.friction.check==1
                plot(x,spatgamX,'b',x,cfSA,'r--',x,ChiAdj,'g',x,bathy/maxbath,'k',...
                    Xinflux,0,'b*',x(par.friction.index),par.friction.Cf,'+y', x,0*x,'b:');
                leg=legend('influx function', 'simulation area','adjustment zone','bathy/maxbathy','influx position','friction coef.');
            else
                plot(x,spatgamX,'b',x,cfSA,'r--',x,ChiAdj,'g',x,bathy/maxbath,'k',...
                    Xinflux,0,'b*', x,0*x,'b:');
                leg=legend('influx function', 'simulation area','adjustment zone','bathy/maxbathy','influx position');
            end
            
        end
    end
    set(leg,'Location','SouthEast','visible','off')
    jbh = handle(hf);
    set(jbh,'buttondownfcn',{@showlegend,leg})
    
    xlabel('x [m]');  xlim([x(1) x(end)]);  ylim([min(bathy/maxbath) 2])
    title('Xspace','fontweight','normal');
    
    
    if IVP.type~=1
        if IVP.type==5
            subplot(6,1,2)
        else
            subplot(5,1,2)
        end
        plot(x,InitVal.eta.*cfSA','r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
        xlim([x(1) x(end)])
        xlabel('x [m]')
        ylabel('\eta [m]')
        title('Initial wave-elevation condition','fontweight','normal')
        if IVP.type==5
            subplot(6,1,3)
            plot(x,InitVal.u.*cfSA','r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
            xlim([x(1) x(end)])
            xlabel('x [m]')
            ylabel('u [m/s]')
            title('Initial tangential velocity condition','fontweight','normal')
        end
    end
    
    if IVP.type==5
        hf2=subplot(6,1,4);
    elseif IVP.type~=5 && IVP.type~=1
        hf2=subplot(5,1,3);
    else
        hf2=subplot(4,1,2);
    end
    kk = k(1:floor(length(k)/2));
    
    if nupeak<=0
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,0,'*b')
    else
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,(0:nupeak:(nFreqP)*nupeak),'*b')
    end
    
    
    if ~strcmp(model.dispersion,'OmExact')
        leg2=legend('Exact',model.dispersion(3:end));
        set(leg2,'Location','SouthEast','visible','off')
        jbh = handle(hf2);
        set(jbh,'buttondownfcn',{@showlegend,leg2})
    end
    
    axis([0 kcut_Simul 0 OmKcut_simul]);
    
    grid on;
    title(['Dispersion @ depth = ',num2str(depth),' [m]'],'fontweight','normal')
    xlabel('k [1/m]');
    ylabel('\omega [rad/s]');
    
    if IVP.type==5
        subplot(6,1,5)
    elseif IVP.type~=5 && IVP.type~=1
        subplot(5,1,4)
    else
        subplot(4,1,3)
    end
    %     subplot(2,1,1)
    %             if length(timesig)<1005, dp=1; ...
    %             elseif length(timesig)>1005 && length(timesig)<5005,dp=2;...
    %             else   dp=4;
    %             end;
    dp=1;
   
    plot(timesig(1:dp:end),Orginsig(1:dp:end),'b')  % Orginsig is original signal  %insig updated signal
    xlim([t_init t_end])
    xlabel('Time [s]')
    ylabel('\eta [m]')
    title('Influx signal','fontweight','normal')
    
    %         subplot(2,1,2)
    if IVP.type==5
        subplot(6,1,6)
    elseif IVP.type~=5 && IVP.type~=1
        subplot(5,1,5)
    else
        subplot(4,1,4);
    end
    if strcmp(bath.wall.check, 'Yes') && par.wall.type==2
        Omk=par.wall.Om;
        refl_coef=par.wall.refl_Coef;
        rhoOm=par.wall.Rho_Om;
        indOmend=floor(length(Omk)/2);
        refl_coefN=interp1(Omk(1:indOmend),refl_coef(1:indOmend),halfomsig);
        refl_coefN(isnan(refl_coefN))=refl_coef(indOmend);
        rhoOmN=interp1(Omk(1:indOmend),rhoOm(1:indOmend),halfomsig);
        rhoOmN(isnan(refl_coefN))=rhoOm(indOmend);
        if indOmend<length(halfomsig);
            refl_coefN(indOmend+1:length(halfomsig))=refl_coef(indOmend);
            rhoOmN(indOmend+1:length(halfomsig))=rhoOm(indOmend);
        end
        [AX,Hlines]=plotyy(halfomsig,abs(Orginsighat(1:floor(length(omsig)/2))),halfomsig,refl_coefN,'plot');
        set(AX(1),'XLim',[0 OmKcut_simul])
        set(AX(2),'XLim',[0 OmKcut_simul])
        hold(AX(2),'on');
        plot(AX(2),halfomsig,rhoOmN)
        hold(AX(2),'off')
        ylabel(AX(1),'Amplitude Sp.');
        %ylabel(AX(2),'$R(\omega)$','interpreter','latex');
        
        ylabel(AX(2),'$R(\omega)$ and $\rho(\omega)$','interpreter','latex');
    else
        
        plot(halfomsig,orginflShat(1:floor(length(omsig)/2)))
        ylabel('Variance density Sp.')
        xlim([0,OmKcut_simul])
    end
    
    grid on;
    xlabel('\omega [rad/s]')
 
    
elseif IVP.type~=1  %Initial Value Problem
    
    hfig=figure('Name','Prepared', 'Position',[10,10,2,600]);
    set(hfig,'Renderer','zbuffer'); %due to graphics driver
    
    hf=subplot(5,1,1);
    if strcmp(bath.wall.check, 'Yes')
        if model.nonlinear>1
        plot(x,cfSA,'r--',x,bathy/maxbath,'k',x,ChiAdj,'g-.',...
            Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
        else
        plot(x,cfSA,'r--',x,bathy/maxbath,'k',...
            Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
        end
        leg=legend('simulation area','bathy/maxbathy','wall position','Location','SouthEast');
    else
        plot(x,cfSA,'r--',x,bathy/maxbath,'k', x,0*x,'b:');
        leg=legend('simulation area','bathy/maxbathy','Location','SouthEast');
    end
    set(leg,'visible','off')
    jbh = handle(hf);
    set(jbh,'buttondownfcn',{@showlegend,leg})
    
    xlim([x(1) x(end)]); ylim([min(bathy/maxbath) 2])
    title('Xspace','fontweight','normal');
    
    
    subplot(5,1,2)
    kk = k(1:length(k)/2);
    
    if nupeak<=0
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,0,'*b')
    else
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,(0:nupeak:(nFreqP)*nupeak),'*b')
    end
    
    axis([0 kcut_Simul 0 OmKcut_simul]);
    
    
    grid on;
    title(['Dispersion @ depth = ',num2str(depth),' [m]'],'fontweight','normal')
    xlabel('k [1/m]');
    ylabel('\omega [rad/s]');
    
    subplot(5,1,3)
    plot(x,InitVal.eta.*cfSA','r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
    xlim([x(1) x(end)])
    xlabel('x [m]')
    ylabel('\eta [m]')
    title('Initial wave-elevation','fontweight','normal')
    subplot(5,1,4)
    Spect_k=abs(fft(InitVal.eta.*cfSA'));
    Nk=floor(length(k)/2);
    plot(k(1:Nk),Spect_k(1:Nk),'r',k(1:Nk),0.*k(1:Nk),':b')  % Orginsig is original signal  %insig updated signal
    xlim([0 kcut_Simul])
    xlabel('k [1/m]')
    ylabel('Amplitude Sp.')
    title('Spectrum of the intial wave profile','fontweight','normal')
    subplot(5,1,5)
    plot(x,InitVal.u.*cfSA','r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
    xlim([x(1) x(end)])
    xlabel('x [m]')
    ylabel('u [m/s]')
    title('Initial tangential velocity','fontweight','normal')
    
 elseif input.assim.check==1  %Assimilation
    
    hfig=figure('Name','Prepared', 'Position',[10,10,2,600]);
    set(hfig,'Renderer','zbuffer'); %due to graphics driver
    hf=subplot(5,1,1);
    if strcmp(bath.wall.check, 'Yes')
        if model.nonlinear>1
        plot(x,cfSA,'r--',x,bdyassim.charupdate,'-b',x,bathy/maxbath,'k',x,ChiAdj,'g-.',...
            Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
        leg=legend('simulation area','Assimilation area','bathy/maxbathy','NonlinAdj.','wall position','Location','SouthEast');
        else
        plot(x,cfSA,'r--',x,bdyassim.charupdate,'-b',x,bathy/maxbath,'k',...
            Xwall,linspace(-1,3/4,4),'oc',x,0*x,'b:');
        leg=legend('simulation area','Assimilation area','bathy/maxbathy','wall position','Location','SouthEast');
        end
        
    else
        if model.nonlinear>1
        plot(x,cfSA,'r--',x,bdyassim.charupdate,'-b',x,bathy/maxbath,'k',x,ChiAdj,'g-.', x,0*x,'b:');
        leg=legend('simulation area','Assimilation area','bathy/maxbathy','NonlinAdj.','Location','SouthEast');    
        else
        plot(x,cfSA,'r--',x,bdyassim.charupdate,'-b',x,bathy/maxbath,'k', x,0*x,'b:');
        leg=legend('simulation area','Assimilation area','bathy/maxbathy','Location','SouthEast');
        end
    end
    set(leg,'visible','off')
    jbh = handle(hf);
    set(jbh,'buttondownfcn',{@showlegend,leg})
    
    xlim([x(1) x(end)]); ylim([min(bathy/maxbath) 2])
    title('Xspace','fontweight','normal');
    
    
    subplot(5,1,2)
    kk = k(1:length(k)/2);
    
    if nupeak<=0
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,0,'*b')
    else
        plot(kk,OmExact(kk,depth,omAdd),'b',kk,Om(kk,depth,omAdd),'r--',...
            influx.k_cut,(0:nupeak:(nFreqP)*nupeak),'*b')
    end
    
    axis([0 kcut_Simul 0 OmKcut_simul]);
    
    
    grid on;
    title(['Dispersion @ depth = ',num2str(depth),' [m]'],'fontweight','normal')
    xlabel('k [1/m]');
    ylabel('\omega [rad/s]');
    
    subplot(5,1,3)
    plot(x,bdyeta.*bdyassim.charupdate,'r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
    xlim([x(1) x(end)])
    xlabel('x [m]')
    ylabel('\eta [m]')
    title('Initial wave-elevation','fontweight','normal')
    subplot(5,1,4)
    Spect_k=abs(fft(bdyeta.*bdyassim.charupdate));
    Nk=floor(length(k)/2);
    plot(k(1:Nk),Spect_k(1:Nk),'r',k(1:Nk),0.*k(1:Nk),':b')  % Orginsig is original signal  %insig updated signal
    xlim([0 kcut_Simul])
    xlabel('k [1/m]')
    ylabel('Amplitude Sp.')
    title('Spectrum of the intial wave profile','fontweight','normal')
    subplot(5,1,5)
    plot(x,bdyeta.*bdyassim.charupdate,'r',x,0*x,'b:')  % Orginsig is original signal  %insig updated signal
    xlim([x(1) x(end)])
    xlabel('x [m]')
    ylabel('u [m/s]')
    title('Initial tangential velocity NOT CORRECT YET','fontweight','normal')
       
end

try
    saveas(hfig,[Proj.Dir,'Preview.fig'])
catch
end

if shippar.check==1
    if strcmp(bath.type,'F')==1
        hfigS=figure('Name','Ship Set-up', 'Position',[170,10,2,600]);
        set(hfigS,'Renderer','zbuffer'); %due to graphics driver
    else
        hfigS=figure('Name','Ship Set-up', 'Position',[330,10,2,600]);
        set(hfigS,'Renderer','zbuffer'); %due to graphics driver
    end
    
    if any(strcmp(shippar.data(:,2),'Fixed'))  % for moving ship
        RadPot=zeros(size(x));
    elseif any(strcmp(shippar.data(:,2),'Surge'))
        RadPot=shippar.rad.psi.surge;
    elseif any(strcmp(shippar.data(:,2),'Heave'))
        RadPot=shippar.rad.psi.heave;
    elseif any(strcmp(shippar.data(:,2),'Pitch'))
        RadPot=shippar.rad.psi.pitch;
    end
    subplot(2,1,1)
    plot(x,shippar.form.shapeXcZc(end,:),'b',x,bathy,'k--',x,shippar.form.cfSA_ship,'r');
    title(['Initial Ship profile above bathymetry'])
    legend('Profile','Bathy')
    ylim([min(bathy) 0.1])
    xlabel('x [m]');ylabel('profile')
    plot_properties;
    
    subplot(2,1,2)
    if any(strcmp(shippar.data(:,2),'Free'))
    plot(x,shippar.rad.psi.surge,'r',x,shippar.rad.psi.heave,'g',x,shippar.rad.psi.pitch,'b',x,shippar.init.phi0,'k--');
     legend('Norm. Rad. pot (Surge)','Norm. Rad. pot (Heav)','Norm. Rad. pot (Pitch)','Initial wave potential')
    else
     plot(x,RadPot,'b',x,shippar.init.phi0,'k--');
      legend('Normalized Radiation Pot.','Initial wave potential')
    end
   
    xlabel('x [m]'), ylabel('potential')
    plot_properties;
    
    try
        saveas(hfigS,[Proj.Dir,'PreviewS.fig'])
    catch
    end
    
end

