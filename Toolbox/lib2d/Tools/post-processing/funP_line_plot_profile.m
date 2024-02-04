function funP_line_plot_profile(tsnap,spatsnap,simuldata,setting,axesfig,h)
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
bathy=simuldata.dom.bathy.profile;

if setting.coarse.check==0
stepS=1;
else
stepS=setting.coarse.param;  
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;
if spatsnap.y_check==1
if setting.spatlim.check==1   
spati=setting.spatlim.param(1);
spatf=setting.spatlim.param(2);    
indspat1=funC_closest(X,spati);indspat2=funC_closest(X,spatf);
else
indspat1=1; indspat2=length(X);   
end
spat1=X;spat2=Y;spat2snap=spatsnap.y_snap;
else
if setting.spatlim.check==1   
spati=setting.spatlim.param(1);
spatf=setting.spatlim.param(2);    
indspat1=funC_closest(Y,spati);indspat2=funC_closest(Y,spatf);
else
indspat1=1; indspat2=length(Y);       
end
spat1=Y;spat2=X;spat2snap=spatsnap.x_snap;
end

if setting.spatlim.check==1
xlimm=[spat1(indspat1) spat1(indspat2)];
else
xlimm=[spat1(1) spat1(end)];   
end

if setting.MTAtlim.check==1
indt1=funC_closest(T,setting.MTAtlim.param(1));
indt2=funC_closest(T,setting.MTAtlim.param(2));
else
indt1=1;
indt2=length(T);  
end

indt=funC_closest(T,tsnap);
indspatsnap=funC_closest(spat2,spat2snap);

if varId==1
varb=output.eta;   
elseif varId==2
 if simuldata.model.phiForm==1   
 varb=output.phi;  
 else
 varb=[];    
 end
elseif varId==3
 if simuldata.model.phiForm==1   
 varb=[];  
 else
 varb=output.u;    
 end
elseif varId==4
 if simuldata.model.phiForm==1   
 varb=[];  
 else
 varb=output.v;    
 end 
end

if ~isempty(varb)
    if setting.MTCcheck==1 || setting.MTTcheck==1
        if spatsnap.y_check==1
            varbnow=squeeze(varb(indt1:indt2,indspatsnap,:));
        else
            varbnow=squeeze(varb(indt1:indt2,:,indspatsnap));
        end
        MTC=max(varbnow,[],1);
        MTT=min(varbnow,[],1);
        ATC=zeros(size(MTC));ATT=zeros(size(MTC));
        Nxm=length(MTC(1,:));
        Tnow=T(indt1:indt2);
        [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
        set(jProgressBar,'Maximum',Nxm, 'Value',0);
        jProgressBar.setStringPainted( true );
        statusbarObj.setText('Calculating Average Temporal Amplitude ...');
        [jbStop]=funGui_Java_stopbutton(statusbarObj);
        jProgressBar.setStringPainted( true );
        tic;
        iter=1;
        
        for ii=1:Nxm
            IDstop=funGui_eventLoopStop(jbStop);
            if IDstop==1, break;end
            varbnow_t=squeeze(varbnow(:,ii));
            [ATC(ii),ATT(ii)]=funSP_find_average_temporal_crest_and_trough(Tnow,varbnow_t);
       
            if mod(iter,floor(0.1*Nxm))==0 || iter==floor(0.01*Nxm)
                set(jProgressBar,'Maximum',Nxm, 'Value',iter);
                ETA=funGui_remain_time(iter,Nxm);
                statusbarObj.setText(['time remaining=', num2str(ETA)]);
            end
            iter=iter+1;
       end
       set(jProgressBar,'Maximum',Nxm, 'Value',Nxm);
       jProgressBar.setVisible(0);
       jbStop.setVisible(0);
       statusbarObj.setText('done.');

    else
        MTC=[];MTT=[];ATC=[];ATT=[];
    end
else
    MTC=[];MTT=[];ATC=[];ATT=[];
end

cla(axesfig);
strc={'r','--g','-.b',':c','*k'};
indstr=1;
hold(axesfig,'off');


for ii=1:length(indt)
if varId==1
varb=squeeze(output.eta(indt(ii),:,:));
varname='elevation';
varlabel=[varname, '[m]'];
elseif varId==2
 if simuldata.model.phiForm==1
 varb=squeeze(output.phi(indt(ii),:,:));
 else
  Velx=squeeze(output.u(indt(ii),:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
  Vely=squeeze(output.v(indt(ii),:,:));
  varb=intgrad2(Velx,Vely,X(2)-X(1),Y(2)-Y(1));
 end
varname='potential'; 
varlabel=[varname, '[m^2/s]'];
elseif varId==3
    if simuldata.model.phiForm==1
        varb=funC_ifft2(1i.*dom.*Kx.*fft2(squeeze(output.phi(indt(ii),:,:))));
    else
        Velx=squeeze(output.u(indt(ii),:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        varb=Velx;
    end
    varname='Velocity (in x dir.)';
    varlabel=[varname, '[m/s]'];
elseif varId==4
    if simuldata.model.phiForm==1
        varb=funC_ifft2(1i.*dom.*Ky.*fft2(squeeze(output.phi(indt(ii),:,:))));
    else
        Vely=squeeze(output.v(indt(ii),:,:));%vel.x(indy1:stepS:indy2,indx1:stepS:indx2);
        varb=Vely;
    end
    varname='Velocity (in y dir.)';
    varlabel=[varname, '[m/s]'];
        
end
    
if input.wall.option==1
    varb(dom.wall.char<1)=NaN;
end
if spatsnap.y_check==1
var_xy=varb(indspatsnap,1:stepS:end);
bathy_xy=bathy(indspatsnap,1:stepS:end);
xlab='x[m]';labb=['_at_y=',num2str(spat2snap),'_'];
else
var_xy=varb(1:stepS:end,indspatsnap);    
bathy_xy=bathy(1:stepS:end,indspatsnap);
xlab='y[m]';labb=['_at_x=',num2str(spat2snap),'_'];
end
if indstr==6,indstr=1;end
if setting.MTCcheck==1 && setting.MTTcheck==1
plot(axesfig,spat1(1:stepS:end),var_xy,strc{indstr},...
    spat1(1:stepS:end),MTC,'--k',spat1(1:stepS:end),MTT,'-.c',...
    spat1(1:stepS:end),ATC,'--m',spat1(1:stepS:end),ATT,'-.y');
varlegend={varname,'MTC','MTT','ATC','ATT'};
elseif setting.MTCcheck==1 && setting.MTTcheck==0
plot(axesfig,spat1(1:stepS:end),var_xy,strc{indstr},...
    spat1(1:stepS:end),MTC,'--k',spat1(1:stepS:end),ATC,'--m');
varlegend={varname,'MTC','ATC'};
elseif setting.MTCcheck==0 && setting.MTTcheck==1
plot(axesfig,spat1(1:stepS:end),var_xy,strc{indstr},...
    spat1(1:stepS:end),MTT,'-.c',spat1(1:stepS:end),ATT,'-.y');
varlegend={varname,'MTT','ATT'};
else
plot(axesfig,spat1(1:stepS:end),var_xy,strc{indstr});
varlegend={varname};
end
hold(axesfig,'on');
indstr=indstr+1;
end

if setting.bathy.cb==1
    hold(axesfig,'on');
    Bscale=setting.bathy.scale;
    plot(axesfig,spat1(1:stepS:end),bathy_xy*Bscale,'-k')
    hold(axesfig,'off');
    varlegend{end+1}='Bathy';
end
legend(axesfig,varlegend)

if setting.zlim.check==1
    zlimm=setting.zlim.param';
    ylim(axesfig,zlimm);
end

xlabel(axesfig,xlab);ylabel(axesfig,varlabel);

xlim(axesfig,xlimm);
axes_properties(axesfig,1);

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
         ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_profile_',varname,labb,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_profile_',varname,labb,CheckId.SaveFig_type]))
        close;
    end
end


end