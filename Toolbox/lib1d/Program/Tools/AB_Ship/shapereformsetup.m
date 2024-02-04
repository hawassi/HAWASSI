% 26 Oct 2015 ship shapes
function [parS] = shapereformsetup(g,x,depth,nup,shippar,SXc,SZc,SThetac)
% draft=shippar.draft+shippar.Z;
% ell=shippar.length/2;
dx=x(2)-x(1);
Nship=shippar.Nship;
shapeXcZ0=shippar.form.shapeXcZ0;

shapeXcZc=zeros(Nship+1,length(x));
dxshapeXcZc=zeros(Nship+1,length(x));

nu.x=zeros(Nship,length(x)); % normal
nu.z=zeros(Nship,length(x)); % normal
nu.theta=zeros(Nship,length(x)); % normal
nutild.x=zeros(Nship,length(x)); % normal
nutild.z=zeros(Nship,length(x)); % normal
nutild.theta=zeros(Nship,length(x)); % normal


Chi=zeros(Nship+1,length(x));
ChiWl=zeros(Nship,length(x));
xShip  =zeros(Nship,6);
SChi=0;SshapeXcZc=0;SdxshapeXcZc=0;SshapeXcZ0=0;
SChiWl=0;SdxChi=0;
dxChi=zeros(Nship+1,length(x));
nuz=0;nux=0;nutheta=0;
nuztild=0;nuxtild=0;nuthetatild=0;
cfSA_ship=0;
XZrot0=shippar.form.XZrot0;
Zrot0=XZrot0(:,2);
Xrot0=XZrot0(:,1);
kappaS=zeros((shippar.Evmode+1),shippar.Nship*2);
if length(depth)>1
    gradDepth=gradient(depth,dx);
end
Dwl=zeros(shippar.Nship,2);
dxDwl=zeros(shippar.Nship,2);
Sdampcoef=zeros(length(x),1);

iterx=1;
for i=1:Nship
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%% update shapedat
    try
        shapeXcZcdati(:,1)=shippar.form.shapedat.data(:,1)+SXc(i);
        shapeXcZcdati(:,2)=shippar.form.shapedat.data(:,2)+SZc(i);
        
        [shapeXcZcdati(:,1), shapeXcZcdati(:,2)]=fun_rotation(shapeXcZcdati(:,1),shapeXcZcdati(:,2),SThetac(i),SXc(i)+shippar.form.xShip(i,2));
        
        shapeXcZc_shifted=zeros(size(x));
        indx1=closest(x,min(shapeXcZcdati(:,1)));
        indx2=closest(x,max(shapeXcZcdati(:,1)));
        
        if x(indx1)<shapeXcZcdati(1,1)
            indx1=indx1+1;
        end
        
        if x(indx2)>shapeXcZcdati(end,1)
            indx2=indx2-1;
        end
        
        shapeXcZc_shifted(indx1:indx2)=interp1(shapeXcZcdati(:,1),shapeXcZcdati(:,2),x(indx1:indx2));
        
        Chi(i,:)           = heaviside(-sign(shapeXcZc_shifted(i,:))-10^(-12));
        shapeXcZc(i,:)     = shapeXcZc_shifted.*Chi(i,:);
        dxshapeXcZc(i,:)   = funOprt_FDGradient1dShip(shapeXcZc(i,:),x,SXc(i)+shippar.form.xShip(i,2));%gradient(shapeXcZc(i,:),dx);
        indxL       = find(Chi(i,:),1,'first')-1;
        xShip(i,1)  = x(indxL);
        indxR       = find(Chi(i,:),1,'last')+1;
        xShip(i,3)  = x(indxR);
        indxC       = indxL+floor((indxR-indxL)/2);
        xShip(i,2)  = x(indxC);
        Zrot0(i)    = Zrot0(i)+SZc(i);
        Xrot0(i)    = x(indxC);
    catch
       %%% no rotation
    %   disp('no rotation applied')
        shapeXcZcdati(:,1)=shippar.form.shapedat.data(:,1)+SXc(i);
        shapeXcZcdati(:,2)=shippar.form.shapedat.data(:,2)+SZc(i);
       
        try
            shapeXcZ0I         = circshift(shapeXcZ0(i,:),floor(SXc(i)/dx));
        catch
             [SXc(i) floor(SXc(i)/dx)]
            shapeXcZ0I         = shapeXcZ0(i,:);
        end
        shapeXcZc_shifted  = (shapeXcZ0I+SZc(i));
        
        Chi(i,:)           = heaviside(-sign(shapeXcZ0I(i,:))-10^(-12));
        shapeXcZc(i,:)     = shapeXcZc_shifted.*Chi(i,:);
        dxshapeXcZc(i,:)   = funOprt_FDGradient1dShip(shapeXcZc(i,:),x,SXc(i)+shippar.form.xShip(i,2));%gradient(shapeXcZc(i,:),dx);
        indxL       = find(Chi(i,:),1,'first')-1;
        xShip(i,1)  = x(indxL);
        indxR       = find(Chi(i,:),1,'last')+1;
        xShip(i,3)  = x(indxR);
        indxC       = indxL+floor((indxR-indxL)/2);
        xShip(i,2)  = x(indxC);
        Zrot0(i)    = Zrot0(i)+SZc(i);
        Xrot0(i)    = x(indxC);
    end
  
    
    
    shapeN=shapeXcZc(i,:); shapeN([indxL indxR])=shapeXcZc(i,[indxL+1 indxR-1]);
    
    ChiWl(i,indxL)=1;
    ChiWl(i,indxR)=1;
    SChiWl        = SChiWl+ChiWl(i,:);
    dxChi(i,indxL)= 1/dx;
    dxChi(i,indxR)= -1/dx;
    
    SChi            = SChi+Chi(i,:);
    SdxChi          = SdxChi+dxChi(i,:);
    SshapeXcZc      = SshapeXcZc+shapeXcZc(i,:);
    SdxshapeXcZc    = SdxshapeXcZc+dxshapeXcZc(i,:);
    bf0            = shippar.form.Slength(i)/3;
    cfSA_ship      = cfSA_ship+cf(x,x(indxL+1),bf0)-cf(x,x(indxR-1)-bf0,bf0); % cfSimulationArea
    
    NN=sqrt((Chi(i,:)+ChiWl(i,:)).^2+dxshapeXcZc(i,:).^2);
    NN(NN==0)=1;
    nu.z(i,:)=(Chi(i,:)+ChiWl(i,:))./NN;
    nu.x(i,:)=-dxshapeXcZc(i,:)./NN;
    nu.theta(i,:)=(dxshapeXcZc(i,:).*(shapeN-Zrot0(i))+(x-x(indxC))).*(Chi(i,:)+ChiWl(i,:))./NN;
    nuz      =nuz+nu.z(i,:);
    nux      =nux+nu.x(i,:);
    nutheta  =nutheta+nu.theta(i,:);
    nutild.z(i,:)=Chi(i,:)+ChiWl(i,:);
    nutild.x(i,:)=-dxshapeXcZc(i,:);
    nutild.theta(i,:)=(dxshapeXcZc(i,:).*(shapeN-Zrot0(i))+(x-x(indxC))).*(Chi(i,:)+ChiWl(i,:));
    nuztild      =nuztild+nutild.z(i,:);
    nuxtild      =nuxtild+nutild.x(i,:);
    nuthetatild  =nuthetatild+nutild.theta(i,:);
    
    if length(depth)==1
        DL=depth;DR=depth;
    else
        DL=depth(indxL);DR=depth(indxR);
        dxDL=gradDepth(indxL);dxDR=gradDepth(indxR);
        dxDwl(i,:)=[dxDL dxDR];
    end
    Dwl(i,:)=[DL DR];
    
    kpL=invOmExact(nup,DL);
    kpR=invOmExact(nup,DR);
    sig0=nup.^2./g;
    if shippar.Evmode>=1
        kappaDL=invOmEvanescent(sig0*DL,shippar.Evmode)';
        kappaDR=invOmEvanescent(sig0*DR,shippar.Evmode)';
        kappaL=[-1i*kpL kappaDL./DL].';
        kappaR=[-1i*kpR kappaDR./DR].';
        kappaS(:,iterx:iterx+1)=[kappaL kappaR];
    else
        kappaS(1,iterx:iterx+1)=[-1i*kpL -1i*kpR];
    end
    iterx=iterx+2;
    
    DS=(DL+DR)/2;
    kpS=invOmExact(nup,DS);
    CpS=UpExact(kpS,DS);
    LS=shippar.form.Slength(i);
    Sdampcoef(indxL:indxR)=7*CpS/(LS/3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ndom=shippar.form.shapedat.Ndom; 
    shapedat(i).Ndom=Ndom;
    shapedat(i).data=shapeXcZcdati;
    IdDraftMax=closest(shapeXcZcdati(:,2),min(shapeXcZcdati(:,2)));
    shapedat(i).XshipL=sort(shapeXcZcdati(1:IdDraftMax,1),'descend');
    shapedat(i).ZshipL=sort(shapeXcZcdati(1:IdDraftMax,2));
    shapedat(i).XshipR=shapeXcZcdati(IdDraftMax:end,1);
    shapedat(i).ZshipR=shapeXcZcdati(IdDraftMax:end,2);
    shapedat(i).dxZshipL=gradient(shapedat(i).ZshipL)./gradient(shapedat(i).XshipL);
    shapedat(i).dxZshipR=gradient(shapedat(i).ZshipR)./gradient(shapedat(i).XshipR);
    
    indxBdy=zeros(1,Ndom-1);
    indxBdy(1)=indxL;indxBdy(end)=indxR;
    indxBdy(Ndom/2)=indxC;
    dIndL=floor((indxC-indxL)/((Ndom-2)/2));
    dIndR=floor((indxR-indxC)/((Ndom-2)/2));
    
    for ii=2:Ndom/2-1
        indxBdy(Ndom/2-(ii-1))=indxC-dIndL*(ii-1);
        indxBdy(Ndom/2+(ii-1))=indxC+dIndR*(ii-1);
    end
    shapedat(i).indxBdy=indxBdy;
    
    shipformOnGrid_aprox=shapeXcZc(i,:).';
    
    for ii=1:Ndom/2-1
        shipformOnGrid_aprox(indxBdy(ii)+1:indxBdy(ii+1))= shapeXcZc(i,indxBdy(ii+1));
        shipformOnGrid_aprox(indxBdy(end-ii):indxBdy(end-ii+1)-1)= shapeXcZc(i,indxBdy(end-ii));
    end
    shapedat(i).shapeapproxOnGrid=shipformOnGrid_aprox;
    
    Xbdy=x(indxBdy);
    shapedat(i).Xbdy=Xbdy;
    
    dr=zeros(size(indxBdy));
    
    for ii=1:(Ndom-2)/2
        dr(ii)=-shipformOnGrid_aprox(indxBdy(ii+1));% shipform(indxC) shipform(indxC) shipform(indxC) shipform(indxcr+1)];
        dr(end-ii+1)=-shipformOnGrid_aprox(indxBdy(end-ii));
    end
    dr(Ndom/2)=-shipformOnGrid_aprox(indxBdy(Ndom/2));
    
    shapedat(i).dr_approx=dr;
    
    
end
Chi(i+1,:)=SChi;
ChiWl(i+1,:)=SChiWl;
dxChi(i+1,:)=SdxChi;
shapeXcZc(i+1,:)=SshapeXcZc;
dxshapeXcZc(i+1,:)=SdxshapeXcZc;
nu.z(i+1,:)=nuz;
nu.x(i+1,:)=nux;


parS=shippar.form; %take the initial setup
parS.shapedat=shapedat;
parS.shapeXcZc=shapeXcZc;
parS.dxshapeXcZc=dxshapeXcZc;
parS.chi=Chi';
parS.dxchi=dxChi';%
parS.chiWl=ChiWl';
parS.xShip=xShip;
parS.nu=nu;
parS.nutild=nutild;
parS.cfSA_ship=cfSA_ship;
parS.kappa=kappaS;
parS.XZrot0=[Xrot0 Zrot0];
parS.kappa=kappaS;
parS.Dwl=Dwl;
parS.dxDwl=dxDwl;
parS.Sdampcoef=Sdampcoef;

% eta=eta.*(1-chiship(:,end))+shippar.reform.shapeXcZc(end,:)';
% figure
% plot(x,SshapeXcZc,'b',x,parS.chi(:,end),'k--',x,SdxshapeXcZc,'r-.');
% title(['Mass (total)=',num2str(sum(Mass)),'; LWL (total)= ',num2str(sum(LWL))])
% plot_properties;

