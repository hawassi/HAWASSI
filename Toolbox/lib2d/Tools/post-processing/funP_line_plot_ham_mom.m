function funP_line_plot_ham_mom(simuldata,setting,axesfig,h)
varId=setting.var;
output=simuldata.output;
input=simuldata.input;
par=simuldata.par;
dom=simuldata.dom;
model=simuldata.model;
bath=simuldata.bath;
body=input.body;
if strcmpi(input.wave.option,'No')
ivp=simuldata.ivp;
k_p=ivp.par.k_p;
meandepth=min(ivp.par.meandepth);
nupeak=ivp.par.nupeak;
else
influx=simuldata.influx;
k_p=influx.par.k_p;
meandepth=min(influx.par.meandepth);
nupeak=max(influx.par.nu_p);
end

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
if setting.tlim.check==1
indt1=funC_closest(T,setting.tlim.param(1));
indt2=funC_closest(T,setting.tlim.param(2));
else
indt1=1;indt2=length(T);   
end
xlimm=[T(indt1) T(indt2)];   
Tnow=T(indt1:stepT:indt2);
dts=T(2)-T(1);
Ntt=floor((T(indt2)-T(indt1))/dts/stepT);
    
cla(axesfig);
OperatorSetup2D;
Ep=zeros(Ntt,1);
Ek=zeros(Ntt,1);
Mx=zeros(Ntt,1);
My=zeros(Ntt,1);
xx=X(indx1:indx2);
yy=Y(indy1:indy2);   
stepy=1;stepx=1;
[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
set(jProgressBar,'Maximum',Ntt, 'Value',0);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Calculating...');
[jbStop]=funGui_Java_stopbutton(statusbarObj);
jProgressBar.setStringPainted( true );
tic;
iter=1;
for ii=1:Ntt
    IDstop=funGui_eventLoopStop(jbStop);
    if IDstop==1, break;end
    
indtti=indt1+stepT*(ii-1);
if indtti>=indt2, break; end
eta=squeeze(output.eta(indtti,:,:));
if model.phiForm==1
phi=squeeze(output.phi(indtti,:,:));
else
u=squeeze(output.u(indtti,:,:));  
v=squeeze(output.v(indtti,:,:));  
phi=intgrad2(u,v,dom.dx,dom.dy);
end
phi_hat=fft2(phi);
if varId==1
    if  strcmpi(bath.name, 'Flat')
        if strcmp(model.evol,'HS')
            Lphi_hat      =Oprt.L2d.*phi_hat;
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
            Ep(ii)=trapz(yy,trapz(xx,P(indy1:stepy:indy2,indx1:stepx:indx2),2));
            Ek(ii)=trapz(yy,trapz(xx,K(indy1:stepy:indy2,indx1:stepx:indx2),2));
        end
    else
        if strcmp(model.evol,'HS')
            Lphi_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,phi_hat,phi);
            Lphi=funC_ifft2(Lphi_hat);
            K2=0.5.*(phi.*Lphi);
            P=0.5.*par.g.*eta.^2;
             if model.nonlinear==1
                K=K2;
             else
                 etaLphi       = eta.*Lphi;
                 etaLphi_hat   = fft2(etaLphi);
                 LetaLphi_hat  = funOprt_L2d(g,Oprt.Om2dSqInterp,etaLphi_hat,etaLphi);
                 gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
                 gradphi2      = funOprt_innerproduct(gradphi,gradphi);
                 K3=eta.*gradphi2-phi.*funC_ifft2(LetaLphi_hat);
                 if model.nonlinear==2
                    K=K2+K3;
                 else
                     LetaLphi        =funC_ifft2(LetaLphi_hat);
                     etaLetaLphi     =eta.*LetaLphi;
                     etaLetaLphi_hat =fft2(etaLetaLphi);
                     LetaLetaLphi_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,etaLetaLphi_hat,etaLetaLphi);
                     eta2Lphi_hat=fft2(eta.^2.*Lphi);
                     gradeta2Lphi_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,eta2Lphi_hat);
                     DivGradeta2Lphi_hat=funOprt_div2d(dom.Kx,dom.Ky,gradeta2Lphi_hat);
                     Gradphi_hat    =funOprt_grad2d_hat(dom.Kx,dom.Ky,phi_hat);
                     DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,Gradphi_hat);
                     DivGradphi     =funC_ifft2(DivGradphi_hat);
                     eta2Divgradphi =eta.^2.*DivGradphi;
                     eta2Divgradphi_hat=fft2(eta2Divgradphi);
                     Leta2DivGradphi_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,eta2Divgradphi_hat,eta2Divgradphi);
                     delphiH3_hat  =LetaLetaLphi_hat+0.5.*(DivGradeta2Lphi_hat+Leta2DivGradphi_hat);
                     K4            =phi.*funC_ifft2(delphiH3_hat);
                    if model.nonlinear==3
                        K=K2+K3+K4;
                    end
                 end
             end
            Ep(ii)=trapz(yy,trapz(xx,P(indy1:stepy:indy2,indx1:stepx:indx2),2));
            Ek(ii)=trapz(yy,trapz(xx,K(indy1:stepy:indy2,indx1:stepx:indx2),2));
        end
    end
else
    Depth  =-dom.bathy.profile ;
    Momx=(Depth+eta).*funC_ifft2(1i.*dom.Kx.*phi_hat);
    Momy=(Depth+eta).*funC_ifft2(1i.*dom.Ky.*phi_hat);
    Mx(ii)=trapz(yy,trapz(xx,Momx(indy1:stepy:indy2,indx1:stepx:indx2),2));
    My(ii)=trapz(yy,trapz(xx,Momy(indy1:stepy:indy2,indx1:stepx:indx2),2));
end 
if mod(iter,floor(0.1*Ntt))==0 || iter==floor(0.01*Ntt)
        set(jProgressBar,'Maximum',Ntt, 'Value',iter);
        ETA=funGui_remain_time(iter,Ntt);
        statusbarObj.setText(['time remaining=', num2str(ETA)]);
end
    iter=iter+1;
end
set(jProgressBar,'Maximum',Ntt, 'Value',Ntt); 
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');


if varId==1
if Ep(end)==0
        IndEnd=find(Ep(2:end)==0,1,'first')-1;
    else
        IndEnd=Ntt;
end
Ham=Ep(1:IndEnd)+Ek(1:IndEnd);

plot(axesfig,Tnow(1:IndEnd),Ham,'r',Tnow(1:IndEnd),Ep(1:IndEnd),'--g',Tnow(1:IndEnd),Ek(1:IndEnd),'-.b');

legend(axesfig,'Hamiltonian','Potential Energy','Kinetic Energy');
xlabel(axesfig,'time[s]');ylabel(axesfig,'Energy[Nm]');
maxvar=max(Ham); 
minvar=min(Ham); 
varname='Hamiltonian';
else
if Mx(end)==0
        IndEnd=find(Mx(2:end)==0,1,'first')-1;
    else
        IndEnd=Ntt;
end
plot(axesfig,Tnow(1:IndEnd),Mx(1:IndEnd),'r',Tnow(1:IndEnd),My(1:IndEnd),'--b');
legend(axesfig,'Momentum (in x direction)','Momentum (in y direction)');
xlabel(axesfig,'time[s]');ylabel(axesfig,'Momentum [Kg m/s]');
maxvar=max(max(Mx),max(My)); 
minvar=min(min(Mx),min(My)); 
varname='Momentum';
end

set(jProgressBar,'Maximum',Ntt, 'Value',Ntt); 
jProgressBar.setVisible(0);
jbStop.setVisible(0);
statusbarObj.setText('done.');


if setting.zlim.check==0
zlimm=[minvar maxvar];
if minvar==maxvar
zlimm=[minvar-0.01 maxvar+0.01];    
end
else
zlimm=setting.zlim.param';    
end
axes_properties(axesfig,1);
xlim(axesfig,xlimm);ylim(axesfig,zlimm);
axes_properties(axesfig,1);

if CheckId.SaveFig==1
    if strcmp(CheckId.SaveFig_type,'.eps')
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',varname,CheckId.SaveFig_type]),'epsc');
        close
    else
        ax_old = axesfig;
        f_new = figure('unit','normalized','position',[1.1 1.1 0.7 0.7]);
        ax_new = copyobj(ax_old,f_new);
        axes_properties(ax_new,1);
        set(f_new,'unit','normalized','position',[1.1 1.1 0.7 0.7])
        set(ax_new,'outerposition',[0 0 1 1]);
        saveas(ax_new,cell2mat([Proj.workdir,'PP_line_',varname,CheckId.SaveFig_type]))
        close;
    end
end

end