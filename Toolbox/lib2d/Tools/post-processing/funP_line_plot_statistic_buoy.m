function funP_line_plot_statistic_buoy(spatsnap,simuldata,setting,axesfig,h)
hfig=h.fig;
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
X=output.X;Y=output.Y;
T=output.time;
if setting.coarse.check==0
stepS=1;
else
stepS=setting.coarse.param;  
end
CheckId.SaveFig=setting.savefig.check;
CheckId.SaveFig_type=setting.savefig.format;



if setting.tlim.check==1
indt1=funC_closest(T,setting.tlim.param(1));
indt2=funC_closest(T,setting.tlim.param(2));
else
indt1=1;indt2=length(T);   
end

tt=output.time(indt1:indt2);dt=tt(2)-tt(1);
Nt=length(tt);Ndat=length(spatsnap.buoy_x_snap);

indX=funC_closest(X,spatsnap.buoy_x_snap);
indY=funC_closest(Y,spatsnap.buoy_y_snap);
Eta=zeros(Nt*Ndat,1);
iter=1;
for ii=1:Ndat
    Eta(iter:iter+Nt-1)=squeeze(output.eta(indt1:indt2,indY(ii),indX(ii)));
    iter=iter+Nt;
end
T=linspace(tt(1),length(Eta)*dt,length(Eta)).';
Eta=Eta(1:stepS:end);
if varId==1
var_xy=4*sqrt(var(Eta));
elseif varId==2
var_xy=max(Eta);
elseif varId==3
var_xy=min(Eta);
elseif varId==4
[var_xy,~]=funSP_find_average_temporal_crest_and_trough(T,Eta);
elseif varId==5
[~,var_xy]=funSP_find_average_temporal_crest_and_trough(T,Eta);
elseif varId==6
var_xy=mean(Eta);
elseif varId==7
var_xy=funSP_find_peak_period(T,Eta);
elseif varId==8
var_xy=funSP_find_skewness(Eta);
elseif varId==9
var_xy=-funSP_find_asymmetry(Eta);
elseif varId==10
var_xy=kurtosis(Eta-mean(Eta));
end

if varId<11    
set(h.pp_line_statistic_show_value_edit,'string',var_xy,'fontweight','bold')
else
    %Analyze simulated wave
    
    [waveheights_AB_eta,T1_AB_eta] = WaveExtremeAnalysis(T,Eta');
    %Determine probabilities of exceedance for measured wave heights
    [Probability_AB_eta] = WaveExtremeProbability(waveheights_AB_eta.CrestTrough,waveheights_AB_eta.TroughCrest,waveheights_AB_eta.Crest,waveheights_AB_eta.Trough);
    
    [waveheights_AB_mineta,T1_AB_mineta] = WaveExtremeAnalysis(T,-Eta');
    %Determine probabilities of exceedance for measured wave heights
    [Probability_AB_mineta] = WaveExtremeProbability(waveheights_AB_mineta.CrestTrough,waveheights_AB_mineta.TroughCrest,waveheights_AB_mineta.Crest,waveheights_AB_mineta.Trough);
   
    
    %Plot distributions on log scale
    if setting.Hsref.check==0
    Hs_eta=waveheights_AB_eta.SignificantHm0;
   % Hs_mineta=waveheights_AB_mineta.SignificantHm0;
    else
    Hs_eta=setting.Hsref.param;  
    end
    Hs_mineta=Hs_eta; 
   

    semilogy(axesfig,waveheights_AB_eta.Crest/Hs_eta,Probability_AB_eta.Crest,'or',...
        waveheights_AB_mineta.Crest/Hs_mineta,Probability_AB_mineta.Crest,'*b');

    
    if setting.savedat==1
        varrcrest=[(waveheights_AB_eta.Crest/Hs_eta)' Probability_AB_eta.Crest'];
        varrtrough=[(waveheights_AB_mineta.Crest/Hs_mineta)' Probability_AB_mineta.Crest'];
        save([Proj.workdir,'PP_line_statistic_crest_exceedance','.txt'],'varrcrest','-ascii');
        save([Proj.workdir,'PP_line_statistic_trough_exceedance','.txt'],'varrtrough','-ascii');
    end
    
    xlabel(axesfig,'|\eta|/H_{s}ref. [-]')
    ylabel(axesfig,'Probability of exceedance [-]')
    legend(axesfig,'Crest exc.','Trough exc.');
   


   % title(['wave probe at X= ' num2str(data223002F.waveprobes(i))])
   if setting.spatlim.check==1
       spati=setting.spatlim.param(1);
       spatf=setting.spatlim.param(2);
       xlimm=[spati spatf];
       xlim(axesfig,xlimm);
   end
   
   if setting.zlim.check==1
      zlimm=[setting.zlim.param(1) setting.zlim.param(2)];
      ylim(axesfig,zlimm);
   end
   
  

axes(axesfig);
axes_properties(axesfig,1);

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        legend(ax_new,'Crest exc.','Trough exc.');
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_crest_trough_exceedance',CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        legend(ax_new,'Crest exc.','Trough exc.');
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_statistic_crest_trough_exceedance',CheckId.SaveFig_type]))
        close;
    end
end
end
end