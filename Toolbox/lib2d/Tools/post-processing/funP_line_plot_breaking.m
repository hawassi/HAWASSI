function funP_line_plot_breaking(simuldata,setting,axesfig,h)

varId=setting.var;
output=simuldata.output;
input=simuldata.input;
par=simuldata.par;
dom=simuldata.dom;
model=simuldata.model;
break_crest=output.break_crest;
break_speed=output.break_speed;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
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




cla(axesfig);
if varId==1
    if setting.tlim.check==1
        indt1=funC_closest(T,setting.tlim.param(1));
        indt2=funC_closest(T,setting.tlim.param(2));
    else
        indt1=funC_closest(T,T(1));
        indt2=funC_closest(T,T(end));
    end

    Nrow=length(break_crest);
    iter=1;
    for jj=1:Nrow
      Nx=length(break_crest(jj).X);
      tbreak(iter:iter+Nx-1)=break_crest(jj).t;
      Xbreak(iter:iter+Nx-1)=break_crest(jj).X;
      Ybreak(iter:iter+Nx-1)=break_crest(jj).Y;
      iter=iter+Nx+1;
    end
    plot3(axesfig,Xbreak,Ybreak,tbreak,'ob');
    xlabel(axesfig,'x [m]');
    ylabel(axesfig,'y [m]');
    zlabel(axesfig,'t [s]');
    xlimm=[X(indx1) X(indx2)];   
    ylimm=[Y(indy1) Y(indy2)];
    zlimm=[T(indt1) T(indt2)];
    xlim(axesfig,xlimm);ylim(axesfig,ylimm);
    zlim(axesfig,zlimm);
    axes_properties(axesfig,1);
    varname='position';
else
    tt=break_speed(:,1);
    Ccrest=break_speed(:,2);
    Ucrest=break_speed(:,3);  
    
    
    if setting.tlim.check==1
        indt1=funC_closest(T,setting.tlim.param(1));
        indt2=funC_closest(T,setting.tlim.param(2));
    else
        indt1=funC_closest(T,tt(1));
        indt2=funC_closest(T,tt(end));
    end
    
    
    
    
    plot(axesfig,tt(1:stepT:end),Ucrest(1:stepT:end),'*r',tt(1:stepT:end),Ccrest(1:stepT:end),'db');
    xlabel(axesfig,'t [s]');
    ylabel(axesfig,'Speed [m/s]');
    legend(axesfig,'U','C')
    xlimm=[T(indt1) T(indt2)];   
    xlim(axesfig,xlimm);
    
    if setting.zlim.check==1
        ylim(axesfig,[setting.zlim.param(1) setting.zlim.param(2)]);
    
    else
         ylim(axesfig,[min(min(Ucrest),min(Ccrest)) max(max(Ucrest),max(Ccrest))]);

    end
    axes_properties(axesfig,1);
    varname='kinematics'; 
end

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_breaking_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_breaking_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end