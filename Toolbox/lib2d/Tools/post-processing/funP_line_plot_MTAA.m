function funP_line_plot_MTAA(simuldata,setting,axesfig,h)
output=simuldata.output;
varId=setting.var;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;

if varId==1
    varP=output.eta;
    varname='elevation';
    varlabel=[varname, ' [m]'];
elseif varId==2
    if simuldata.model.phiForm==1
        varP=output.phi;
        varname='potential';
        varlabel=[varname, ' [m^2/s]'];
    else
        varP=output.u;
        varname='velocity (in x dir)';
        varlabel=[varname, ' [m/s]'];
    end
else
    varP=output.v;
    varname='velocity (in y dir)';
    varlabel=[varname, ' [m/s]'];
end

if setting.coarse.check==0
    stepT=1;
else
    stepT=setting.coarse.param;
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;

if setting.xlim.check==1
    indx1=funC_closest(X,setting.xlim.param(1));
    indx2=funC_closest(X,setting.xlim.param(2));
else
    indx1=1;indx2=length(X);
end

if setting.ylim.check==1
    indy1=funC_closest(Y,setting.ylim.param(1));
    indy2=funC_closest(Y,setting.ylim.param(2));
else
    indy1=1;indy2=length(Y);
end

if setting.tlim.check==1
    indt1=funC_closest(T,setting.tlim.param(1));
    indt2=funC_closest(T,setting.tlim.param(2));
else
    indt1=1;indt2=length(T);
end

if setting.plotopt==1
    Xlabel='time[s]';
    Ylabel='MAA';
    varplot='Area';
    savevar='MAAvsT';
    IterVar=T;
    inditer1=indt1;inditer2=indt2;
elseif setting.plotopt==2
    IterVar=X;inditer1=indx1;inditer2=indx2;
    Xlabel='x[m]';
    Ylabel='MTA';
    varplot='Temporal';
    savevar='MTAvsX';
elseif setting.plotopt==3
    IterVar=Y;inditer1=indy1;inditer2=indy2;
    Xlabel='y[m]';
    Ylabel='MTA';
    varplot='Temporal';
    savevar='MTAvsY';
elseif setting.plotopt==4
    Xlabel='x[m]';
    Ylabel='y[m]';
    Zlabel='t[s]';
    savevar='MACrest_position_';
    IterVar=T;
    inditer1=indt1;inditer2=indt2;
elseif setting.plotopt==5
    Xlabel='x[m]';
    Ylabel='y[m]';
    Zlabel='t[s]';
    savevar='MATrough_position_';
    IterVar=T;
    inditer1=indt1;inditer2=indt2;
elseif setting.plotopt==6
    Xlabel='time[s]';
    Ylabel='AAA';
    varplot='Area';
    savevar='AAAvsT';
    IterVar=T;
    inditer1=indt1;inditer2=indt2;
elseif setting.plotopt==7
    IterVar=X;inditer1=indx1;inditer2=indx2;
    Xlabel='x[m]';
    Ylabel='ATA';
    varplot='Temporal';
    savevar='ATAvsX';
elseif setting.plotopt==8
    IterVar=Y;inditer1=indy1;inditer2=indy2;
    Xlabel='y[m]';
    Ylabel='ATA';
    varplot='Temporal';
    savevar='ATAvsY';
end

if setting.plotopt>=1 && setting.plotopt<=3
    
    if varId==1
        Ylabel=[Ylabel,' [m]'];
    else
        if simuldata.model.phiForm==1
            Ylabel=[Ylabel,' [m^2/s]'];
        else
            Ylabel=[Ylabel,' [m/s]'];
        end
    end
end
xlimm=[IterVar(inditer1) IterVar(inditer2)];
%dts=IterVar(2)-IterVar(1);
Niter=length(IterVar);%floor((IterVar(inditer2)-IterVar(inditer1))/dts/stepT);


cla(axesfig);


[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
set(jProgressBar,'Maximum',Niter, 'Value',0);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Calculating...');
[jbStop]=funGui_Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
tic;
iter=1;

MTAC=zeros(Niter,1)*NaN;
MTAT=zeros(Niter,1)*NaN;
IterVarnow=zeros(Niter,1);
ATAC=zeros(Niter,1)*NaN;
ATAT=zeros(Niter,1)*NaN;
if setting.plotopt==4
    Xmax=zeros(Niter,1);
    Ymax=zeros(Niter,1);
    MaxVar=zeros(Niter,1);
end
for ii=1:Niter
    
    IDstop=funGui_eventLoopStop(jbStop);
    if IDstop==1, break;end
    
    inditer_i=inditer1+stepT*(ii-1);
    
    if setting.plotopt==1 || setting.plotopt==4 || setting.plotopt==5
        if inditer_i>indt2, break; end
        varTi=squeeze(varP(inditer_i,indy1:indy2,indx1:indx2));
    elseif  setting.plotopt==2
        if inditer_i>indx2, break; end
        varTi=squeeze(varP(indt1:indt2,indy1:indy2,inditer_i));
    elseif  setting.plotopt==3
        if inditer_i>indy2, break; end
        varTi=squeeze(varP(indt1:indt2,inditer_i,indx1:indx2));
    elseif  setting.plotopt==6
        if inditer_i>indt2, break; end
        varTi=squeeze(varP(inditer_i,indy1:indy2,indx1:indx2));
        Nyi=length(varTi(:,1)); Nxi=length(varTi(1,:));
        if Nyi<Nxi
          ATAC_y=zeros(Nyi,1);
          ATAT_y=zeros(Nyi,1);
          xnow=X(indx1:indx2);
          for jj=1:Nyi
              varbnow_x=squeeze(varTi(jj,:));
              [ATAC_y(jj),ATAT_y(jj)]=funSP_find_average_temporal_crest_and_trough(xnow.',varbnow_x.');
          end
          ATAC(ii)=mean(ATAC_y); ATAT(ii)=mean(ATAT_y);
        else
            ATAC_x=zeros(Nxi,1);
            ATAT_x=zeros(Nxi,1);
            ynow=Y(indy1:indy2);
            for jj=1:Nxi
                varbnow_y=squeeze(varTi(:,jj));
                [ATAC_x(jj),ATAT_x(jj)]=funSP_find_average_temporal_crest_and_trough(ynow.',varbnow_y);
            end
            ATAC(ii)=mean(ATAC_x); ATAT(ii)=mean(ATAT_x);
       end
    elseif  setting.plotopt==7
        if inditer_i>indx2, break; end
        varTi=squeeze(varP(indt1:indt2,indy1:indy2,inditer_i));
        Nyi=length(varTi(1,:));
        ATAC_y=zeros(Nyi,1);
        ATAT_y=zeros(Nyi,1);
        tnow=T(indt1:indt2);
        for jj=1:Nyi
            varbnow_t=squeeze(varTi(:,jj));
            [ATAC_y(jj),ATAT_y(jj)]=funSP_find_average_temporal_crest_and_trough(tnow,varbnow_t);
        end
        ATAC(ii)=mean(ATAC_y); ATAT(ii)=mean(ATAT_y);
    elseif setting.plotopt==8
        if inditer_i>indy2, break; end
        varTi=squeeze(varP(indt1:indt2,inditer_i,indx1:indx2));
        Nxi=length(varTi(1,:));
        ATAC_x=zeros(Nxi,1);
        ATAT_x=zeros(Nxi,1);
        tnow=T(indt1:indt2);
        for jj=1:Nxi
            varbnow_t=squeeze(varTi(:,jj));
            [ATAC_x(jj),ATAT_x(jj)]=funSP_find_average_temporal_crest_and_trough(tnow,varbnow_t);
            
        end
        ATAC(ii)=mean(ATAC_x); ATAT(ii)=mean(ATAT_x);
    end
  
    
    if  setting.plotopt<6
        if setting.plotopt==4  || setting.plotopt==5
            if setting.plotopt==4
                [maxA,ind] = max(varTi(:));
            else
                [maxA,ind] = min(varTi(:));
            end
            [m,n] = ind2sub(size(varTi),ind);
            Xmax(ii)=X(n);
            Ymax(ii)=Y(m);
            MaxVar(ii)=maxA;
        else
            MTAC(ii)=max(max(varTi));
            MTAT(ii)=min(min(varTi));
        end
    end
    IterVarnow(ii)=IterVar(inditer_i);
    
    if mod(iter,floor(0.1*Niter))==0 || iter==floor(0.01*Niter)
        set(jProgressBar,'Maximum',Niter, 'Value',iter);
        ETA=funGui_remain_time(iter,Niter);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
    end
    iter=iter+1;
end

set(jProgressBar,'Maximum',Niter, 'Value',Niter);
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');
if setting.plotopt==4 || setting.plotopt==5
    %plot3(axesfig,Xmax,Ymax,Tnow,'or');
    scatter3(axesfig, Xmax(:), Ymax(:), IterVarnow(:), [],MaxVar(:), 'filled' )
    view([-40 10]);
    colormap(axesfig,'jet');
    cb=colorbar;
    axes_properties_cb(cb,1);
    cb.Label.String=varlabel;
    axes_properties(axesfig,1);
    axes(axesfig)
    axes(axesfig)
    %cbfreeze;cblabel(varname);
    %freezeColors;
    zlabel(axesfig,Zlabel);
else
    if  setting.plotopt<6
        if IDstop==0
            plot(axesfig,IterVarnow,MTAC,'r',IterVarnow,MTAT,'--b');
        else
            plot(axesfig,IterVarnow(1:inditer_i-1),MTAC(1:inditer_i-1),'r',IterVarnow(1:inditer_i-1),MTAT(1:inditer_i-1),'--b');
        end
          title(axesfig,['Maximal ', varplot ,' Amplitude of ', varname])
    else
        if IDstop==0
            plot(axesfig,IterVarnow,ATAC,'r',IterVarnow,ATAT,'--b');
        else
            plot(axesfig,IterVarnow(1:inditer_i-1),ATAC(1:inditer_i-1),'r',IterVarnow(1:inditer_i-1),ATAT(1:inditer_i-1),'--b');
        end
          title(axesfig,['Average ', varplot ,' Amplitude of ', varname])
    end
    legend(axesfig,[Ylabel(1:2),'C'],[Ylabel(1:2),'T']);
end

ylabel(axesfig,Ylabel);
xlabel(axesfig,Xlabel);

if  setting.plotopt<6
maxvar=max(MTAC);
minvar=min(MTAT);
else
maxvar=max(ATAC);
minvar=min(ATAT);    
end

if setting.zlim.check==0
    zlimm=[minvar maxvar];
    if minvar==maxvar
        zlimm=[minvar-0.01 maxvar+0.01];
    end
else
    zlimm=setting.zlim.param';
end

axes_properties(axesfig,1);

if setting.plotopt==4 || setting.plotopt==5
    xlim(axesfig,[X(indx1) X(indx2)]);ylim(axesfig,[Y(indy1) Y(indy2)]);
    zlim(axesfig,[T(indt1) T(indt2)]);
else
    xlim(axesfig,xlimm);ylim(axesfig,zlimm);
end
axes_properties(axesfig,1);


if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        if setting.plotopt==4 || setting.plotopt==5
            colormap(ax_new,'jet');
            cb=colorbar;
            ylabel(cb,varname,'fontweight', 'bold')
        end
        axes_properties(ax_new,1);
        axes_properties(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',savevar,'_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        set(f_new,'visible','off')
        axes(ax_new)
        if setting.plotopt==4 || setting.plotopt==5
            colormap(ax_new,'jet');
            cb=colorbar;
            ylabel(cb,varname,'fontweight', 'bold')
        end
        axes_properties(ax_new,1);
        axes_properties(cb,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',savevar,'_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end