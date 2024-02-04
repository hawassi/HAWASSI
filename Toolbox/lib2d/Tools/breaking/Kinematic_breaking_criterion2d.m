%%Kinematic Breaking Criteria
%%Wave will be break if particle speed exceeds crest speed
%%===============================================================

function [CB] = Kinematic_breaking_criterion2d(tnow,eta,gradphi,U,dom,param,g)
plotCheck=0;

global itCB iterdtcheck dataSpeedBreak ITERbspdt ITERbdt
CB=[];

XX=dom.XX;
YY=dom.YY;

if any(eta(param.char==1)>param.MaxEtaInit)
        EtaChecked=eta.*param.char;
        PeaksChar=imregionalmax(EtaChecked);
        PeaksChar(EtaChecked<param.MaxEtaInit)=0;
        IndPeaks=find(PeaksChar==1);
    if all(PeaksChar==0)
        flag=0;
    else
        flag=1;
    end
else
    flag=0;
end


if flag==1
   dx=dom.X(2)-dom.X(1);
   dy=dom.Y(2)-dom.Y(1);
   
    Npks=length(IndPeaks);
    %HilbEta=imag(hilbert(eta));
    HilbEta=funC_ifft2(-1i.*sign(dom.Kx+dom.Ky).*fft2(eta));
    
    
    [gradxHilbEta,gradyHilbEta]=gradient(HilbEta,dx,dy);
    [gradxEta,gradyEta]=gradient(eta,dx,dy);
    Kx=(eta.*gradxHilbEta-HilbEta.*gradxEta)./(eta.^2+HilbEta.^2);
    Ky=(eta.*gradyHilbEta-HilbEta.*gradyEta)./(eta.^2+HilbEta.^2);
    Kk=sqrt(Kx.^2+Ky.^2);
    

    itCB=1;

    Totdepth=-dom.bathy.profile+eta; 
    DirPropTheta=atan(gradphi.y./gradphi.x);
    for j=1:Npks
        indJ  =IndPeaks(j);
        Ccrest=funOprt_UpExact(Kk(indJ),Totdepth(indJ),g);  
        if abs(U(indJ))>=(param.KBC*abs(Ccrest));
           %  [Kk(indJ) Ccrest U(indJ)]            
            kwadran=funBr_check_kwadran_dir(gradphi.x(indJ),gradphi.y(indJ));
            CB.index(itCB)   =indJ;
            CB.Ucrest(itCB)  =abs(U(indJ));
            CB.DirProp.theta(itCB)=abs(DirPropTheta(indJ));
            CB.DirProp.kwadran(itCB)=kwadran;
            CB.position.x(itCB)=dom.XX(indJ);
            CB.position.y(itCB)=dom.YY(indJ); 
            CB.time(itCB)    =tnow;
            itCB=itCB+1;
        end
    end
    
        if tnow>ITERbdt*param.dt 
        Etanow=eta(IndPeaks);
        [~,ind]=max(Etanow);
        dataSpeedBreak(ITERbspdt,1)=tnow;
        dataSpeedBreak(ITERbspdt,2)=funOprt_UpExact(Kk(IndPeaks(ind)),Totdepth(IndPeaks(ind)),g);  
        dataSpeedBreak(ITERbspdt,3)=abs(U(IndPeaks(ind)));
        ITERbspdt=ITERbspdt+1;
        end   

    
    

    if plotCheck==1
            if  iterdtcheck==1
                figure;
            end
            CCcrest=funOprt_UpExact(Kk(IndPeaks),Totdepth(indJ),g);
            subplot(3,1,1)
            surf(XX,YY,EtaChecked,'Edgecolor','none')
            xlabel('x');ylabel('y');
            view([30,60])
            hold on
            plot3(XX(IndPeaks),YY(IndPeaks),EtaChecked(IndPeaks),'or')
            xlabel('x');ylabel('y');
            hold off
            subplot(3,1,2)
            surf(XX,YY,U,'Edgecolor','none')
            view([30,60])
            xlabel('x');ylabel('y');
            hold on
            plot3(XX(IndPeaks),YY(IndPeaks),U(IndPeaks),'or',XX(IndPeaks),YY(IndPeaks),CCcrest,'*g')
            xlabel('x');ylabel('y');
            hold off
            pause(0.001)
            iterdtcheck=iterdtcheck+1;
    end
end


