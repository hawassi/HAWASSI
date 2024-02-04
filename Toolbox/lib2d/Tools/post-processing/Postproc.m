%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% HAWASSI-AB 2D                                         %%%%%%%%%%
%%%%%%%%%% Hamiltonian Wave-Ship-Structure Interaction           %%%%%%%%%%
%%%%%%%%%% Copyright (c): LabMath-Indonesia                      %%%%%%%%%%
%%%%%%%%%% version: 5 July 2016                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Post-Processing                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExecutePP=7; %1=2D profile at t=t1 
             %2=time signal at (x,y)=(x1,y1) v
             %3=1D  profile at t=t1 y=y1; 
             %4=1D  profile at t=t1 x=x1; 
             %5=2D  MTC;
             %6=2D  MTT;
             %7=2D  Hs;
             %8=Momentum and Energy
             %9=validation with measurement
             %10=2D animation
             %11=1D animatio at y=yi;
             %12= HS in Buoy

varproc=1;



if varproc==1
    varb=output.eta;
    varname='wave elevation';
else
    varb=output.phi;
    varname='wave potential';
end
if input.wall.option==1
    varb(:,dom.wall.char<1)=NaN;
end

X=output.X;Y=output.Y;T=output.time;
[XX,YY]=meshgrid(X,Y);
if ExecutePP==1
    tsnap=60;
    viewplot=2;%[-10 10];
    colmap=2;
    xi=0;xf=20;
    yi=Y(1);yf=30;
    indx1=funC_closest(X,xi);indx2=funC_closest(X,xf);
    indy1=funC_closest(Y,yi);indy2=funC_closest(Y,yf);
    CheckId.SaveFig=1;stepS=1;
    
    
    [XX,YY]=meshgrid(X(1:stepS:end),Y(1:stepS:end));
    CheckId.SaveFig_type='.png';
    xlimm=[X(indx1) X(indx2)];ylimm=[Y(indy1) Y(indy2)];
    indt=funC_closest(T,tsnap);
    
    if input.wall.option==1
        varb(:,dom.wall.char<1)=NaN;
    end
    
    var_xy=squeeze(varb(indt,1:stepS:end,1:stepS:end));
     maxvar=max(max(max(var_xy))); climm=[-maxvar maxvar];
    hf1=figure;
    set(hf1, 'units','normalized','Position', [0.2, 0.2, 0.6,0.6]);
    set(hf1,'Renderer','zbuffer'); %due to graphics driver
    surf(XX,YY,var_xy,'edgecolor','none');
    xlabel('x[m]');ylabel('y[m]');zlabel('elevation [m]');
    grid off;
    xlim(xlimm);ylim(ylimm);caxis(climm);
    view(viewplot);
       if colmap==1
            colormap('jet');
        else
            colormap('winter');
        end
    y=colorbar;
    ylabel(y,'elevation [m]','fontweight', 'bold')
    plot_properties;
    
    if CheckId.SaveFig==1
        if strcmp(CheckId.SaveFig_type,'.eps')
            saveas(gcf,[Proj.Dir,'PP1_',varname,'_@t=',num2str(tsnap),CheckId.SaveFig_type],'epsc')
        else
            saveas(gcf,[Proj.Dir,'PP1_',varname,'_@t=',num2str(tsnap),CheckId.SaveFig_type])
        end
    end
    
    
elseif ExecutePP==2
    %     Tpp=5;DD=2;ampl=0.1;
    %     xobs=130;
    %     wpp=2*pi/Tpp;
    %     kpp=invOmExact(wpp,DD);
    %     lambda_pp=2*pi./kpp;
    %     Sol=ampl.*(cos(kpp.*xobs-wpp.*time));
    %     indxobs=funC_closest(dom.X,xobs);
    %     etaSimul=var(:,end/2,indxobs);
    %     figure;
    %     plot(time,Sol,'b',time,etaSimul,'r');
    %     xlabel('time');ylabel('\eta');
    
    xobs=80;
    yobs=19.42;
    indxobs=funC_closest(dom.X,xobs);
    indyobs=funC_closest(dom.Y,yobs);
    
    etaSimul=squeeze(varb(:,indyobs,indxobs));
    figure;
    plot(output.time,etaSimul,'r');
    xlabel('time');ylabel('\eta');
    plot_properties;
    
elseif ExecutePP==3
    tobs=8;
    yobs=15;
    indtobs=funC_closest(output.time,tobs);
    indyobs=funC_closest(dom.Y,yobs);
    
    etaSimul=squeeze(varb(indtobs,indyobs,:));
    
    figure;
    plot(dom.X,etaSimul,'r');
    xlabel('x [m]');ylabel('\eta');
    plot_properties;
    
elseif ExecutePP==4
    tobs=20;
    xobs=0;
    
    indtobs=funC_closest(output.time,tobs);
    indxobs=funC_closest(dom.X,xobs);
    
    etaSimul=squeeze(varb(indtobs,:,indxobs));
    
    figure;
    plot(dom.Y,etaSimul,'r');
    xlabel('y [m]');ylabel('\eta');
    plot_properties;
    
    
elseif ExecutePP==5
    viewplot=[30,60];%viewplot=2;
    xlimm=[X(1) X(end)];ylimm=[Y(1) Y(end)];climm=[-0.2 0.2];
    var_xy=squeeze(max(varb));
    
    hf5=figure;
    set(0,'CurrentFigure',hf5);
    set(gcf,'Renderer','zbuffer'); %due to graphics driver
    surf(XX,YY,var_xy,'edgecolor','none');
    xlabel('x[m]');ylabel('y[m]');
    xlim(xlimm);ylim(ylimm);%caxis(climm);
    view(viewplot);
    colorbar;
    plot_properties;
    
elseif ExecutePP==6
    viewplot=[30,60];%viewplot=2;
    xlimm=[X(1) X(end)];ylimm=[Y(1) Y(end)];climm=[-0.2 0.2];
    var_xy=squeeze(min(varb));
    
    hf6=figure;
    set(0,'CurrentFigure',hf6);
    set(gcf,'Renderer','zbuffer'); %due to graphics driver
    surf(XX,YY,var_xy,'edgecolor','none');
    xlabel('x[m]');ylabel('y[m]');
    xlim(xlimm);ylim(ylimm);%caxis(climm);
    view(viewplot);
    colorbar;
    plot_properties;
    
elseif ExecutePP==7
    viewplot=2;%[30,60];%viewplot=2;
    xlimm=[X(1) X(end)];
    ylimm=[Y(1) Y(end)];
    indti=funC_closest(output.time,T(1));
    indtf=funC_closest(output.time,800);
    stept=2;
    
    var_xy=squeeze(4*sqrt(var(varb(indti:stept:indtf,:,:))));
    climm=[0 0.15];%max(max(var_xy))];
    
    CheckId.SaveFig=1;
    CheckId.SaveFig_type='.png';
    
    hf5=figure;
    set(0,'CurrentFigure',hf5);
    set(gcf,'Renderer','zbuffer'); %due to graphics driver
    surf(XX,YY,var_xy,'edgecolor','none');
    xlabel('x[m]');ylabel('y[m]');
    xlim(xlimm);ylim(ylimm);caxis(climm);
    view(viewplot);
    cb=colorbar;
    ylabel(cb,'Significant wave height [m]','fontweight','bold')
    plot_properties;
    
    if CheckId.SaveFig==1
        if strcmp(CheckId.SaveFig_type,'.eps')
            saveas(gcf,[Proj.Dir,'PP1','-Hs',CheckId.SaveFig_type],'epsc')
        else
            saveas(gcf,[Proj.Dir,'PP1','-Hs',CheckId.SaveFig_type])
        end
    end
    
elseif ExecutePP==8  %%%%%%%%%%% Energy dan momentum
    viewplot=[30,60];%viewplot=2;
    ti=0;tf=T(end);stept=2;
    xi=X(1);xf=X(end);stepx=1;
    yi=Y(1);yf=Y(end);stepy=1;
    
    xlimm=[ti tf];ylimm=[Y(1) Y(end)];
    
    
    indti=funC_closest(T,ti);indtf=funC_closest(T,tf);
    dts=output.time(2)-output.time(1);
    Ntt=floor((tf-ti)/dts/stept);
    indxi=funC_closest(X,xi);indxf=funC_closest(X,xf);
    indyi=funC_closest(Y,yi);indyf=funC_closest(Y,yf);
    
    OperatorSetup;
    dxdy=(dom.X(2)-dom.X(1))*(dom.Y(2)-dom.Y(1));
    Ep=zeros(Ntt,1);
    Ek=zeros(Ntt,1);
    M=zeros(Ntt,1);
    timenow=zeros(Ntt,1);
    xx=dom.X(indxi:stepx:indxf);
    yy=dom.Y(indyi:stepy:indyf);
    KK=sqrt(dom.Kx.^2+dom.Ky.^2);
    for i=1:Ntt
        indtti=indti+stept*(i-1);
        if indtti>=indtf, break; end;
        eta=squeeze(output.eta(indtti,:,:));
        phi=squeeze(output.phi(indtti,:,:));
        phi_hat=fft2(phi);
        if strcmp(model.evol,'HS')
            Lphi_hat=Oprt.L2d.*phi_hat;
            Lphi          = funC_ifft2(Lphi_hat);
            K2=0.5.*(phi.*Lphi);
            P=0.5.*par.g.*eta.^2;
            if model.nonlinear==1
                K=K2;
            else
                LetaLphi_hat  = Oprt.L2d.*fft2(eta.*Lphi);
                gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
                gradphi2      = funOprt_innerproduct(gradphi,gradphi);
                K3=eta.*gradphi2-phi.*funC_ifft2(LetaLphi_hat);
                if model.nonlinear==2
                    K=K2+K3;
                else
                    LetaLphi        =funC_ifft2(LetaLphi_hat);
                    LetaLetaLphi_hat=Oprt.L2d.*fft2(eta.*LetaLphi);
                    eta2Lphi_hat=fft2(eta.^2.*Lphi);
                    gradeta2Lphi_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,eta2Lphi_hat);
                    DivGradeta2Lphi_hat=funOprt_div2d(dom.Kx,dom.Ky,gradeta2Lphi_hat);
                    Gradphi_hat    =funOprt_grad2d_hat(dom.Kx,dom.Ky,phi_hat);
                    DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,Gradphi_hat);
                    DivGradphi     =funC_ifft2(DivGradphi_hat);
                    Leta2DivGradphi_hat=Oprt.L2d.*fft2(eta.^2.*DivGradphi);
                    delphiH3_hat  =LetaLetaLphi_hat+0.5.*(DivGradeta2Lphi_hat+Leta2DivGradphi_hat);
                    K4            =phi.*funC_ifft2(delphiH3_hat);
                    if model.nonlinear==3
                        K=K2+K3+K4;
                    end
                end
            end
            Ep(i)=trapz(yy,trapz(xx,P(indyi:stepy:indyf,indxi:stepx:indxf),2));
            Ek(i)=trapz(yy,trapz(xx,K(indyi:stepy:indyf,indxi:stepx:indxf),2));
            timenow(i)=output.time(indtti);
            Mom=(bath.depth+eta).*funC_ifft2(1i.*KK.*phi_hat);
            M(i)=trapz(yy,trapz(xx,Mom(indyi:stepy:indyf,indxi:stepx:indxf),2));
        end
    end
    %%
    if Ep(end)==0
        IndEnd=find(Ep(2:end)==0,1,'first')-1;
    else
        IndEnd=Ntt;
    end
    figure
    subplot(2,1,1)
    plot(timenow(1:IndEnd),Ep(1:IndEnd)+Ek(1:IndEnd),'r',timenow(1:IndEnd),Ep(1:IndEnd),'--g',timenow(1:IndEnd),Ek(1:IndEnd),'-.b');
    legend('Hamiltonian','Potential Energy','Kinetic Energy');
    xlabel('time[s]');ylabel('Energy[Nm]');
    %     xlim(xlimm);ylim(ylimm);%caxis(climm);
    plot_properties;
    subplot(2,1,2)
    plot(timenow(1:IndEnd),M(1:IndEnd),'r');
    xlabel('time[s]');ylabel('Momentum [Kg m/s]');
    %     xlim(xlimm);ylim(ylimm);%caxis(climm);
    plot_properties;
    
    
elseif ExecutePP==9
    temp=load(':\Workspace\1. Codes\developer\HaWaSSI\HAWASSI_AB2_160921\Testcases\Deltares T001\MeasDataT001.mat');
    if isstruct(temp)
        namevar = fieldnames(temp);
        DataMeas=temp.(namevar{1});
    else
        DataMeas=temp;
    end
    
%     Xmeas=DataMeas(1,[10 24 5 26 11 27 12]+1);
%     Ymeas=DataMeas(2,[10 24 5 26 11 27 12]+1);
%     Meas=DataMeas(3:end,[10 24 5 26 11 27 12]+1);Tmeas=DataMeas(3:end,1);
%     BuoyNr=[10 24 5 26 11 27 12];
    
%      Xmeas=DataMeas(1,[1 23 24 25 2]+1);
%     Ymeas=DataMeas(2,[1 23 24 25 2]+1);
%     Meas=DataMeas(3:end,[1 23 24 25 2]+1);Tmeas=DataMeas(3:end,1);
%     BuoyNr=[1 23 24 25 2];
    
%     Xmeas=DataMeas(1,[6 7 8 9]+1);
%     Ymeas=DataMeas(2,[6 7 8 9]+1);
%     Meas=DataMeas(3:end,[6 7 8 9]+1);Tmeas=DataMeas(3:end,1);
%     BuoyNr=[6 7 8 9];
    
    Xmeas=DataMeas(1,[13 14 15 16]+1);
    Ymeas=DataMeas(2,[13 14 15 16]+1);
    Meas=DataMeas(3:end,[13 14 15 16]+1);Tmeas=DataMeas(3:end,1);
    BuoyNr=[13 14 15 16];
    
    
    ti=50;tf=100;stept=1;
    MM=0.03;ylimm=[-MM MM];
    
    indtsi=funC_closest(T,ti);indtsf=funC_closest(T,tf);
    Nfig=length(Xmeas);
    hf3a=figure;FlagYlabel=0;
    Corr        = zeros(Nfig,1);
    OptCorr     = zeros(Nfig,2);
    for i=1:Nfig
        subaxis(Nfig,1,Nfig-i+1, 'Spacing', 0, 'Padding', 0, 'Margin', 0.13);
        indtmeas=funC_closest(Tmeas,T(indtsi:stept:indtsf));
        indx=funC_closest(X,Xmeas(i));indy=funC_closest(Y,Ymeas(i));
        Simul=squeeze(varb(:,indy,indx));
        Simulvar=Simul(indtsi:stept:indtsf);
        Measvar=Meas(indtmeas,i);
       % if i==1
          [Xcorel,lags]=xcorr(Simulvar,Measvar,'coeff');
% %        % end
         OptCorr(i,1)=max(Xcorel);
         indMaxCorel=funC_closest(Xcorel,OptCorr(i,1));
        OptCorr(i,2)=lags(indMaxCorel);
        
%         if GUIpp.PP3_timeshift==1
%             if GUIpp.PP3_timeshiftBest==1
                Measvar = circshift(Measvar,OptCorr(i,2));
%             else
             %    Measvar = circshift(Measvar,-floor(10.20/(T(2)-T(1)))-2);
%                 Corr(i)   = (simul'*meas)/sqrt(varmeas*varsimul)/(length(meas));%correl(simul,meas);
%             end
%         end
%         
        
        plot(Tmeas(indtmeas), Measvar,'b',T(indtsi:stept:indtsf), Simulvar,'--r');
        
        xlim([T(indtsi) T(indtsf)]);
        ylim(ylimm);
        
        
        if mod(i,2)==0
            set(gca,'ytick',[]);
        else
            try
            catch
                
            end
        end
        
        if i==1
            set(gca, 'box','off');
            xlabel('time [s]');
            %  legend('Meas','Simul')
        else
            set(gca,'xtick',[],'xcolor','w');
            set(gca,'box','off');
        end
        
        
        
        
        if i==round(Nfig/2)+1% &&  mod(i,2)~=0
%             if any(bathy>0)
%                 ylabel('$\eta-\bar{\eta}$ [m]','Interpreter','latex');
%             else
                ylabel('\eta [m]');
 %           end
            FlagYlabel=1;
        end
        if FlagYlabel==0 && i==round(Nfig/2)+1
%             if any(bathy>0)
%                 ylabel('$\eta-\bar{\eta}$ [m]','Interpreter','latex');
%             else
%                 ylabel('\eta [m]');
%             end
        end
        
        %th1=text(time(indti)+2, MM/2,['x=',num2str(roundn(MeasPos(i),-1))]);
        th1=text(T(indtsi)+1, MM/2,['W',num2str(BuoyNr(i))]);
        plot_properties;
        
    end
    
    
elseif ExecutePP==10
    
    ti=130;tf=T(end);stept=4;
    xi=X(1);xf=X(end);stepx=1;
    yi=Y(1);yf=Y(end);stepy=1;
    
    viewplot=2;%[40,10];%3
    maxvar=4;%abs(max(max(max(varb))));%max(max(max(var)));
    xlimm=[xi xf];ylimm=[yi yf];climm=[-maxvar maxvar];
    zlimm=[-maxvar maxvar];
    saveanim=1;
    colmap=2; %1=jet 2=winter
    
    GIF_delaytime=0.001;%SetAx.GIF.val(1);
    GIF_loopcount=inf;%SetAx.GIF.val(2);
    
    indti=funC_closest(T,ti);indtf=funC_closest(T,tf);
    dts=output.time(2)-output.time(1);
    Ntanim=floor((tf-ti)/dts/stept);
    indxi=funC_closest(X,xi);indxf=funC_closest(X,xf);
    indyi=funC_closest(Y,yi);indyf=funC_closest(Y,yf);
    
    
    SetAx.tcoarse.Id=stept;
    if SetAx.tcoarse.Id==1
        stepdt=stept;%SetAx.tcoarse.val;
    else
        
        stepdt=1;
    end
    
    
    filename = [Proj.Dir,'PP3_animation-',varname,'.gif'];
    hf3=figure;
    numframes=Ntanim;
    set(hf3, 'units','normalized','Position', [0.2, 0.2, 0.6,0.6]);
    
    if input.wall.option==1
        varb(:,dom.wall.char<1)=NaN;
    end
    
    XXs=XX(indyi:stepy:indyf,indxi:stepx:indxf);
    YYs=YY(indyi:stepy:indyf,indxi:stepx:indxf);
    dT=T(2)-T(1);
    for i=1:numframes
        set(0,'CurrentFigure',hf3);
        set(gcf,'Renderer','zbuffer'); %due to graphics driver
        hold off;
        indtti=indti+stept*(i-1);
        if indtti>=indtf, break; end;
        var_xy=squeeze(varb(indtti,indyi:stepy:indyf,indxi:stepx:indxf));
        stringtitle=['Evolution of ',varname,' @ time: ',num2str(roundn(T(indtti),-2)), ' [s]'];
        Maxxeta=max(max(var_xy));
        surf(XXs,YYs,var_xy,'edgecolor','none');
        grid off;
        if colmap==1
            colormap('jet');
        else
            colormap('winter');
        end
        shading interp;
        title(stringtitle)
        xlabel('x[m]');ylabel('y[m]');
        xlim(xlimm);ylim(ylimm);zlim(zlimm);caxis(climm);
        % axis equal
        view(viewplot);
        colorbar;
        plot_properties;
        
        if ~isempty(output.break_nodes)
            if strcmp(model.breaking.check,'Yes') &&T(indtti)>=output.break_nodes(1,1)...
                    &&T(indtti)<=output.break_nodes(end,1)
                indtBreak=funC_closest(output.break_nodes(:,1),T(indtti));
                
                
                
                nodesBreaktemp=output.break_nodes(indtBreak,2:end);
                nodesBreak=nodesBreaktemp(nodesBreaktemp>0);
                
                if all(abs(T(indtti)-output.break_nodes(:,1))>dT)
                else
                    var_xyB=squeeze(varb(indtti,:,:));
                    hold on;
                    plot3(dom.XX([nodesBreak]),dom.YY([nodesBreak]),var_xyB([nodesBreak]),'ow');
                end
                
            end
        end
        
        drawnow
        pause(GIF_delaytime)
        if ~ishandle(hf3)
            break;
        end
        
        try
            frame = getframe(hf3);
        catch
            
        end
        
        if saveanim==1
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 1;
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
            else
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
            end
        end
    end
    
elseif  ExecutePP==11  %%% Animation of a crosssection
    
    ti=0;tf=T(end);stept=1;
    
    xi=5;xf=15;
    stepx=1;
    yi=1;
    
    indyi=funC_closest(Y,yi);
    varb=squeeze(varb(:,indyi,:));
    
    maxvar=abs(min(min(varb)));
    xlimm=[xi xf];ylimm=[-maxvar maxvar];

    saveanim=1;
    
    GIF_delaytime=0.001;%SetAx.GIF.val(1);
    GIF_loopcount=inf;%SetAx.GIF.val(2);
    

    indti=funC_closest(T,ti);indtf=funC_closest(T,tf);
    dts=output.time(2)-output.time(1);
    Ntanim=floor((tf-ti)/dts/stept);
    indxi=funC_closest(X,xi);indxf=funC_closest(X,xf);
    
    SetAx.tcoarse.Id=stept;
    if SetAx.tcoarse.Id==1
        stepdt=stept;%SetAx.tcoarse.val;
    else

        stepdt=1;
    end
    
    
    filename = [Proj.Dir,'PP3_animation_line-',varname,'.gif'];
    hf3=figure;
    numframes=Ntanim;
    set(hf3, 'units','normalized','Position', [0.2, 0.2, 0.6,0.6]);
    
    if input.wall.option==1
       % var(:,dom.wall.char<1)=NaN;
    end
    
    XXs=X(indxi:stepx:indxf);
    dT=T(2)-T(1);
    for i=1:numframes
        set(0,'CurrentFigure',hf3);
        set(gcf,'Renderer','zbuffer'); %due to graphics driver
        hold off;
        indtti=indti+stept*(i-1);
        if indtti>=indtf, break; end;
        var_y=squeeze(varb(indtti,indxi:stepx:indxf));
        stringtitle=['Evolution of ',varname,'@y=',num2str(yi),'[m], @ time: ',num2str(roundn(T(indtti),-2)), ' [s]'];
        Maxeta=max(var_y);
        plot(XXs,var_y,'r');
        grid off;

        title(stringtitle)
        xlabel('x[m]');ylabel([varname,' [m]']);
        xlim(xlimm);ylim(ylimm);

        plot_properties;
        
%         if ~isempty(output.break_nodes)
%             if strcmp(model.breaking.check,'Yes') &&T(indtti)>=output.break_nodes(1,1)...
%                     &&T(indtti)<=output.break_nodes(end,1)
%                 indtBreak=funC_closest(output.break_nodes(:,1),T(indtti));
%                 
%                 
%                 
%                 nodesBreaktemp=output.break_nodes(indtBreak,2:end);
%                 nodesBreak=nodesBreaktemp(nodesBreaktemp>0);
%                 
%                 if all(abs(T(indtti)-output.break_nodes(:,1))>dT)
%                 else
%                     var_xyB=squeeze(var(indtti,:,:));
%                     hold on;
%                     plot3(dom.XX([nodesBreak]),dom.YY([nodesBreak]),var_xyB([nodesBreak]),'ow');
%                 end
%                 
%             end
%         end
        
        drawnow
        pause(GIF_delaytime)
        if ~ishandle(hf3)
            break;
        end
        
        try
            frame = getframe(hf3);
        catch
            
        end
        
        if saveanim==1
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 1;
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
            else
                imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
            end
        end
    end

    elseif  ExecutePP==12
    indti=funC_closest(output.time,10);
    indtf=funC_closest(output.time,120);
    stept=2;
    xobs=[10.7 21 16.4 28 24.5 17.2 20.8 24.5];
    yobs=[36 35.3 28 28 26 33.6 30 31.5];
    
    Hsbuoy=zeros(length(xobs),1);
    for i=1:length(xobs)
    indy=funC_closest(dom.Y,yobs(i));
    indx=funC_closest(dom.X,xobs(i));
    
    var_xy=squeeze(4*sqrt(var(varb(indti:stept:indtf,indy,indx))));
    Hsbuoy(i)=var_xy;
    end
   Hsbuoy 
end
