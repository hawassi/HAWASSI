function funP_anim_line(spatsnap,simuldata,setting,axesfig,h)
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
dom=simuldata.dom;
Proj=simuldata.Proj;
bathy=simuldata.dom.bathy.profile;

X=output.X;Y=output.Y;
T=output.time;
if setting.coarse.check==0
stepT=1;
else
stepT=setting.coarse.param;  
end
stepS=1;

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

if setting.tlim.check==1
    indt1=funC_closest(T,setting.tlim.param(1));
    indt2=funC_closest(T,setting.tlim.param(2));
else
    indt1=1;indt2=length(T);
end
Tnow=T(indt1:stepT:indt2);
    
if setting.saveanim.check==1
  GIF_delaytime=setting.gifset.param(1);
  GIF_loopcount=setting.gifset.param(2);
saveanim=1;
else
saveanim=0;
end

indspatsnap=funC_closest(spat2,spat2snap);

cla(axesfig);
axes(axesfig)
hold(axesfig,'off');
 if saveanim==1
%  ax = axesfig;
%  set(ax,'units','pixels');
%  pos = get(ax,'Position');
%  ti = get(ax,'TightInset');
%  rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
 
 ax = axesfig;
 set(ax,'units','pixels');
 pos =get(ax,'Position');
 ti = get(ax,'TightInset');
 marg=ti(1);
 rect=[-marg, -marg, pos(3)+1.5*marg, pos(4)+1.5*marg];
 set(ax,'units','normalized');
 end
 if varId==1
     if spatsnap.y_check==1
         varS=squeeze(output.eta(indt1:stepT:indt2,indspatsnap,:));
         xlab='x[m]';labb=[', @y=',num2str(spat2snap),'[m], '];
         
     else
         varS=squeeze(output.eta(indt1:stepT:indt2,:,indspatsnap));
         xlab='y[m]';labb=[', @x=',num2str(spat2snap),'[m], '];
     end
     varname='elevation'; varlabel=[varname,' [m]'];
 elseif varId==2
     if simuldata.model.phiForm==1
         if spatsnap.y_check==1
             varS=squeeze(output.phi(indt1:stepT:indt2,indspatsnap,:));
             xlab='x[m]';labb=[', @y=',num2str(spat2snap),'[m], '];
             
         else
             varS=squeeze(output.phi(indt1:stepT:indt2,:,indspatsnap));
             xlab='y[m]';labb=[', @x=',num2str(spat2snap),'[m], '];
         end
         varname='potential'; varlabel=[varname,' [m^2/s]'];
     else
         if spatsnap.y_check==1
             varS=squeeze(output.u(indt1:stepT:indt2,indspatsnap,:));
             xlab='x[m]';labb=[', @y=',num2str(spat2snap),'[m], '];
             
         else
             varS=squeeze(output.u(indt1:stepT:indt2,:,indspatsnap));
             xlab='y[m]';labb=[', @x=',num2str(spat2snap),'[m], '];
         end
         varname='velocity (in x dir.)'; varlabel=[varname,' [m/s]']; 
     end
 elseif varId==3
     if spatsnap.y_check==1
         varS=squeeze(output.v(indt1:stepT:indt2,indspatsnap,:));
         xlab='x[m]';labb=[', @y=',num2str(spat2snap),'[m], '];
         
     else
         varS=squeeze(output.v(indt1:stepT:indt2,:,indspatsnap));
         xlab='y[m]';labb=[', @x=',num2str(spat2snap),'[m], '];
     end
     varname='velocity (in y dir.)'; varlabel=[varname,' [m/s]'];
 end
 maxvar=max(max(varS)); minvar=min(min(varS));  
if setting.MTCcheck==1 || setting.MTTcheck==1 
 MTC=max(varS,[],1);
 MTT=min(varS,[],1);
 
 ATC=zeros(size(MTC));ATT=zeros(size(MTC));
 Nxm=length(MTC(1,:));
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
     varbnow_t=squeeze(varS(:,ii));
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
 

IdLegend=0;
for ii=1:length(Tnow)
indt=funC_closest(T,Tnow(ii));
  
if varId==1
varb=squeeze(output.eta(indt,:,:));
elseif varId==2
    if simuldata.model.phiForm==1
        varb=squeeze(output.phi(indt,:,:));
    else
        varb=squeeze(output.u(indt,:,:));
    end
else
    varb=squeeze(output.v(indt,:,:));
end
    
if input.wall.option==1
    varb(dom.wall.char<1)=NaN;
end
if spatsnap.y_check==1
var_xy=varb(indspatsnap,1:stepS:end);
bathy_xy=bathy(indspatsnap,1:stepS:end);
else
var_xy=varb(1:stepS:end,indspatsnap); 
bathy_xy=bathy(1:stepS:end,indspatsnap);
end
if setting.MTCcheck==1 && setting.MTTcheck==1
plot(axesfig,spat1(1:stepS:end),var_xy,'r',...
    spat1(1:stepS:end),MTC,'--k',spat1(1:stepS:end),MTT,'-.c',...
    spat1(1:stepS:end),ATC,'--m',spat1(1:stepS:end),ATT,'-.y');
IdLegend=1;
varlegend={varname,'MTC','MTT','ATC','ATT'};
elseif setting.MTCcheck==1 && setting.MTTcheck==0
plot(axesfig,spat1(1:stepS:end),var_xy,'r',...
    spat1(1:stepS:end),MTC,'--k', spat1(1:stepS:end),ATC,'--m');
varlegend={varname,'MTC','ATC'};
IdLegend=1;
elseif setting.MTCcheck==0 && setting.MTTcheck==1
plot(axesfig,spat1(1:stepS:end),var_xy,'r',...
    spat1(1:stepS:end),MTT,'-.c',spat1(1:stepS:end),ATT,'-.y');
varlegend={varname,'MTT','ATT'};
IdLegend=1;
else
plot(axesfig,spat1(1:stepS:end),var_xy,'r'); 
varlegend={varname};
IdLegend=0;
end

if setting.bathy.cb==1
    hold(axesfig,'on');
    Bscale=setting.bathy.scale;
    plot(axesfig,spat1(1:stepS:end),bathy_xy*Bscale,'-k')
    hold(axesfig,'off');
    varlegend{end+1}='Bathy';
end
if IdLegend==1
legend(axesfig,varlegend);
end

stringtitle=['Evolution of ',varname,labb,' @ time: ',num2str(roundn(T(indt),-2)), ' [s]'];
title(axesfig,stringtitle);
if setting.zlim.check==1
    zlimm=setting.zlim.param';
    ylim(axesfig,zlimm);
else
    if setting.bathy.cb==1
    ylim(axesfig,[min(minvar,min(bathy_xy*Bscale)) max(maxvar,max(bathy_xy*Bscale))]);   
    else
    ylim(axesfig,[minvar maxvar]);
    end
end
xlim(axesfig,xlimm);
xlabel(axesfig,xlab);ylabel(axesfig,varlabel);
axes_properties(axesfig,1);
drawnow limitrate
if saveanim==1
    drawnow limitrate;
    frame = getframe(ax,rect);
    pause(GIF_delaytime)
%     pause(GIF_delaytime)
    
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename=[Proj.workdir,'PP_anim_line_',varname,'.gif'];
    if ii == 1
        imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
    end
else
    drawnow limitrate;
end

if get(h.pp_anim_line_pause,'userdata')==1
    uiwait(gcf);
end
if get(h.pp_anim_line_stop,'userdata')==1
    set(h.pp_anim_line_stop,'userdata',0);
    break;
end
end


end