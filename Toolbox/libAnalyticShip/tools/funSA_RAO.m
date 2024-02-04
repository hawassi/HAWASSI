function datRAO=funSA_RAO(proj,par,ship,wave,bottom,calc)
%tempF=load([proj.workdir,'\datCalForces.mat']);
F0x        =calc.difraddata.Forces(:,2);
F0z        =calc.difraddata.Forces(:,3);
F0theta    =calc.difraddata.Forces(:,4);
PhaseFx    =calc.difraddata.Forces(:,5);
PhaseFz    =calc.difraddata.Forces(:,6);
PhaseFtheta=calc.difraddata.Forces(:,7);

Aw_inc=par.Ainc;
we=par.w0;ke=par.k0;
M_mat=par.M_mat;
M=M_mat(1,1);
I=M_mat(3,3);
GM=ship.GM;
g=par.g;
B=ship.width;
draft=ship.draft;
Cd2=ship.visc_sway;
Cd4=ship.visc_roll;
L=1;

bnu2=Cd2*L*draft/2;
bnu4=Cd4*B^4/2;



Nk=length(ke);



%tempA=load([proj.workdir,'\datCalAddedMass.mat']);
addedMass=calc.difraddata.AddedMass;
%tempD=load([proj.workdir,'\datCalDampCoef.mat']);
dampCoef=calc.difraddata.DampCoef;
RestForces=[0;g*B;g*M*GM];

F0=[F0x F0z F0theta];
RAO=zeros(Nk,7);
Am=addedMass(:,[2 6 10]);
Dc=dampCoef(:,[2 6 10]);
RAO(:,1)=we;
for IdMot=1:3
    aa=Am(:,IdMot);
    bb=Dc(:,IdMot);
    cc=RestForces(IdMot);
    Mass=M_mat(IdMot,IdMot);
    if ship.mooring.check==1
        Tnnow= ship.mooring.Tn(IdMot);
        if Tnnow>0
            cc=(Mass+aa).*(2*pi/Tnnow)^2;
        end
    end
    
    
    F0exc=F0(:,IdMot);
    
    Xi=F0exc./sqrt((-(Mass+aa).*we.^2+cc).^2+(bb.*we).^2);
    RAO(:,IdMot+1)=Xi./Aw_inc;
    if IdMot==3
        RAO(:,IdMot+1)=RAO(:,IdMot+1)./ke;
    end
    RAO(:,IdMot+4)=atan(bb.*we)./(-(Mass+aa).*we.^2+cc);
end


options = optimoptions('fsolve','display','iter');%,'StepTolerance',1e-8);

RAOCM=zeros(Nk,7);
RAOCM(:,1)=we;

global Idstop

for ii=1:Nk
    X0=[0,0,0,0,0,0];%[0.5*Aw_inc(ii),0.5*Aw_inc(ii),0.5*Aw_inc(ii)*ke(ii),pi/2,pi/2,pi/2]*0;
    addMass_mat=[addedMass(ii,2:4);addedMass(ii,5:7);addedMass(ii,8:10)];
    dampCoef_mat=[dampCoef(ii,2:4);dampCoef(ii,5:7);dampCoef(ii,8:10)];
    Cvect=RestForces;
    if ship.mooring.check==1
        for IdMot=1:3
            Tnnow= ship.mooring.Tn(IdMot);
            if Tnnow~=0
                Mass=M_mat(IdMot,IdMot);
                Cvect(IdMot)=(Mass+addMass_mat(IdMot)).*(2*pi/Tnnow)^2;
            end
        end
    end
    
    F0ExcVect=F0(ii,:);
    PhaseExcVect=[PhaseFx(ii), PhaseFz(ii), PhaseFtheta(ii)];
    Xi=fsolve(@(XX)funSA_EOMsystem(XX,we(ii),M,I,addMass_mat, dampCoef_mat,Cvect,F0ExcVect,PhaseExcVect,bnu2,bnu4),X0,options);
    RAOCM(ii,2)=abs(Xi(1)/Aw_inc(ii));
    RAOCM(ii,3)=abs(Xi(2)/Aw_inc(ii));
    RAOCM(ii,4)=abs(Xi(3)/Aw_inc(ii)/ke(ii));
    RAOCM(ii,5)=Xi(4);
    RAOCM(ii,6)=Xi(5);
    RAOCM(ii,7)=Xi(6);
    
    if mod(ii,floor(0.1*Nk))==0
        set(par.jProgressBar,'Maximum',Nk, 'Value',ii);
        par.jProgressBar.setStringPainted( true );
        ETA=remain_time(ii,Nk);
        par.statusbarObj.setText(['time=', num2str(ETA)]);
    end
    
    fun_Stop_iter_button(par.jbStop);
    
    if Idstop==1
        set(par.jProgressBar,'Maximum',Nk, 'Value',Nk);
        par.jProgressBar.setStringPainted( true );
        break;
    end
end
if Idstop==1
    par.statusbarObj.setText(['Terminated.']);
else
    par.statusbarObj.setText(['']);
end

datRAO.notcoupled=RAO;
datRAO.coupled=RAOCM;



