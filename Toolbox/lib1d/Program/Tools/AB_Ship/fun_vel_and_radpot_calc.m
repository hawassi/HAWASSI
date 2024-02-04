function [phiRad,Vel,dtetaS,FpsiS,dxi_K]=fun_vel_and_radpot_calc ...
    (shippar,chiship,chiwl,nutild,Pz,Px,Ptheta,beta,Smotion,Nship,dxi_K)


%%%%% normalized rad pot && added Mass
Vel.x        =zeros(Nship,1);
Vel.z        =zeros(Nship,1);
Vel.theta    =zeros(Nship,1);
dtetaS=0;phiRad=0;FpsiS=0;

for i=1:Nship
    if strcmp(Smotion(i),'Free')
        PiVect=[Px(i);Pz(i);Ptheta(i)];
        BetaVect=[beta.x(i);beta.z(i);beta.theta(i)];
        Svel= (shippar.form.MassS(i).Mat+shippar.rad.Ma.S(i).Mat)\(PiVect-BetaVect);
        Vel.x(i)=Svel(1);
        Vel.z(i)=Svel(2);
        Vel.theta(i)=Svel(3);
    elseif strcmp(Smotion(i),'Surge')
        Vel.x(i)=(Px(i)-beta.x(i))./(shippar.form.Mass(i)+shippar.rad.Ma.x(i));
    elseif strcmp(Smotion(i),'Heave')
        Vel.z(i)=(Pz(i)-beta.z(i))./(shippar.form.Mass(i)+shippar.rad.Ma.z(i));
    elseif strcmp(Smotion(i),'Pitch')
        Vel.theta(i)=(Ptheta(i)-beta.theta(i))./(shippar.form.MIner(i)+shippar.rad.Ma.theta(i));
    end
    
    if strcmp(Smotion(i),'Heave')||strcmp(Smotion(i),'Free')
        phiRad0     = Vel.z(i)*shippar.rad.psi.heave.*chiship(:,i);
        dtetaS0     = Vel.z(i).*chiship(:,i);
        phiRad      = phiRad+phiRad0;
        dtetaS      = dtetaS+dtetaS0;
        FpsiS       = FpsiS+shippar.rad.Fpsi.heave.'*Vel.z(i).*(chiship(:,i)+chiwl(:,i));
        dxi_K.z(i)  = dxi_K.z(i)*Vel.z(i).^2;
    end
    
    if strcmp(Smotion(i),'Surge')||strcmp(Smotion(i),'Free')
        phiRad0     = Vel.x(i)*shippar.rad.psi.surge.*chiship(:,i);
        dtetaS0     = Vel.x(i).*nutild.x(i,:).'.*chiship(:,i);
        phiRad      = phiRad+ phiRad0;
        dtetaS      = dtetaS+dtetaS0;
        FpsiS       = FpsiS+shippar.rad.Fpsi.surge.'*Vel.x(i).*(chiship(:,i)+chiwl(:,i));  
        dxi_K.x(i)  = dxi_K.x(i)*Vel.x(i).^2;
    end
    if strcmp(Smotion(i),'Pitch')||strcmp(Smotion(i),'Free')
        phiRad0     = Vel.theta(i)*shippar.rad.psi.pitch.*chiship(:,i);
        dtetaS0     = Vel.theta(i).*nutild.theta(i,:).'.*chiship(:,i);
        phiRad      = phiRad+ phiRad0;
        dtetaS      = dtetaS+ dtetaS0;
        FpsiS       = FpsiS+shippar.rad.Fpsi.pitch.'*Vel.theta(i).*(chiship(:,i)+chiwl(:,i));
        dxi_K.theta(i)  = dxi_K.theta(i)*Vel.theta(i).^2;
    end
    
end
end