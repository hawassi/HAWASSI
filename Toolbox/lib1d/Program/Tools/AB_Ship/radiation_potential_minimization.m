function radpar=radiation_potential_minimization(shippar,par,model,bath,Oprt,statusbarObj,jProgressBar)


set(jProgressBar,'Maximum',100, 'Value',10);
jProgressBar.setStringPainted( true );
statusbarObj.setText('Setting up the ship...');

chiship=shippar.form.chiXcZ0(:,end);
zeta=shippar.form.shapeXcZ0(end,:)';

Nship=shippar.Nship;
nu=0;
for ii=1:Nship
    if strcmp(shippar.data(ii,2),'Heave')
    nu=nu+shippar.form.nuXcZ0.z(ii,:)';
    elseif strcmp(shippar.data(ii,2),'Surge')
    nu=nu+shippar.form.nuXcZ0.x(ii,:)';  
    elseif strcmp(shippar.data(ii,2),'Fixed')
    chiship=chiship.*(1-shippar.form.chi(:,ii));
    end
end


% psi_hat=funOprt_Linv(par.g,par.k,Cp,zeta,nu,model.nonlinear);
% psi=Ifft(psi_hat);
% DtNpsi=Ifft(funOprt_L(par.g,par.k,Cp,zeta,psi,model.nonlinear));

Error_fun=@(x,nu,psi,DtNpsi) trapz(x,(DtNpsi-nu).*psi);
grad_Obj_fun=@(nu,DtNpsi,chiship) (DtNpsi-nu).*chiship;
tol=0.0001;
error=1;
psi=ones(size(par.k)).*chiship;
 Nship=shippar.Nship;
SL=cell2mat(shippar.data(1,3));
Upmin=UpExact(pi/min(SL),min(par.depth+zeta));
dtf=(2*min(SL)/Upmin)/100;
aal=1;%Oprt.aal;
plotnow=0;
if plotnow==1
figure
end
ModNonli=1;%model.nonlinear;

ii=1;

while error>tol    
    if strcmp(bath.type,'B')
    DtNpsi=Ifft(funOprtBathy_L(par.g,par.k,Oprt,zeta,psi,ModNonli,1).*aal);     
    else
    DtNpsi=Ifft(funOprt_L(par.k,Oprt.L0,zeta,psi,ModNonli,1).*aal); 
    end
    
    gradF=grad_Obj_fun(nu,DtNpsi,chiship);
    psi=psi-gradF*dtf;
    
    if strcmp(bath.type,'B')
    DtNpsi=Ifft(funOprtBathy_L(par.g,par.k,Oprt,zeta,psi,ModNonli,1));     
    else
    DtNpsi=Ifft(funOprt_L(par.k,Oprt.L0,zeta,psi,ModNonli,1)); 
    end
    
    error=abs(Error_fun(par.x,nu,psi,DtNpsi));
    
    if error<1-ii/8 &&ii<8
        set(jProgressBar,'Maximum',100, 'Value',10+10*ii);
        jProgressBar.setStringPainted( true );
        ii=ii+1;
    end
    
  if plotnow==1
    subplot(2,1,1)
    plot(par.x,psi,'r');
    title(['Error= ',num2str(error)])
   % xlim([-30 30])
    subplot(2,1,2)
    plot(par.x,DtNpsi,'r');
    %xlim([-30 30])
    pause(0.001)
  end
end

        set(jProgressBar,'Maximum',100, 'Value',90);
        jProgressBar.setStringPainted( true );

AddedMassZ=zeros(Nship,1);
AddedMassX=zeros(Nship,1);
for i=1:Nship
AddedMassZ(i)=trapz(par.x,psi.*shippar.form.nuXcZ0.z(i,:)');
AddedMassX(i)=trapz(par.x,psi.*shippar.form.nuXcZ0.x(i,:)');
end


radpar.psi=psi;
radpar.zMa=AddedMassZ;
radpar.xMa=AddedMassX;
% figure;
% subplot(2,1,1)
% plot(par.x,psi,'r')
% subplot(2,1,2)
% plot(par.x,DtNpsi,'r')
statusbarObj.setText('');
end