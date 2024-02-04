% 26 Oct 2015 ship shapes
function [parS] = shapeformsetup(x,shippar,Xc,Zc,Thetac,bathName)
% draft=shippar.draft+shippar.Z;
% ell=shippar.length/2;
Nship=shippar.Nship;
%shipform=shippar.shipform;
shapeXcZc=zeros(Nship+1,length(x));
dxshapeXcZc=zeros(Nship+1,length(x));
shapeXcZ0=zeros(Nship+1,length(x));
dxshapeXcZ0=zeros(Nship+1,length(x));
nuXcZ0.x=zeros(Nship,length(x)); % normal
nuXcZ0.z=zeros(Nship,length(x)); % normal
nuXcZ0.The=zeros(Nship,length(x)); % normal

nu.x=zeros(Nship,length(x)); % normal
nu.z=zeros(Nship,length(x)); % normal
nu.theta=zeros(Nship,length(x)); % normal
nutild.x=zeros(Nship,length(x)); % normal
nutild.z=zeros(Nship,length(x)); % normal
nutild.theta=zeros(Nship,length(x)); % normal


Mass=zeros(Nship,1);MIner=zeros(Nship,1);
LWL=zeros(Nship,1);
dx=x(2)-x(1);
Chi=zeros(Nship+1,length(x));
dxChi=zeros(Nship+1,length(x));
ChiWl=zeros(Nship,length(x));
ChiSl=zeros(Nship,length(x));
xShip  =zeros(Nship,3);
SChi=0;SChiWl=0;SChiSl=0;SdxChi=0;
Slength=zeros(Nship,1);
Sdraft=zeros(Nship,1);
Zrot0=zeros(Nship,1);Xrot0=zeros(Nship,1);
Gyradius=zeros(Nship,1);
GM     =zeros(Nship,1);
Ndom=zeros(Nship,1);
ZCoG=zeros(Nship,1);
CdX=zeros(Nship,1);
CdZ=zeros(Nship,1);
CdTheta=zeros(Nship,1);
cfSA_ship=zeros(size(x));
SshapeXcZc=0;SdxshapeXcZc=0;SshapeXcZ0=0;
nuz=0;nux=0;nutheta=0;
nuztild=0;nuxtild=0;nuthetatild=0;
for i=1:Nship
    shipform=shippar.data{i,1};
    try
        Slength(i)=str2num(cell2mat(shippar.data(i,3)));
        Sdraft(i)=str2num(cell2mat(shippar.data(i,4)));
    catch
        Slength(i)=cell2mat(shippar.data(i,3));
        Sdraft(i)=cell2mat(shippar.data(i,4));
    end
     if strcmpi(shippar.data(i,2),'Pitch')||strcmpi(shippar.data(i,2),'Free') ||strcmpi(shippar.data(i,2),'Fixed')
         try
            Gyradius(i)=str2num(cell2mat(shippar.data(i,8)));
         catch
            Gyradius(i)=cell2mat(shippar.data(i,8));
         end       
         try
             ZCoG(i)=str2num(cell2mat(shippar.data(i,9)));
         catch
             ZCoG(i)=cell2mat(shippar.data(i,9));
         end
         try
             GM(i)=str2num(cell2mat(shippar.data(i,10)));
         catch
             GM(i)=cell2mat(shippar.data(i,10));
         end
    else
       Gyradius(i)=0;
       Zrot0(i)=0;
       GM(i)=0;
     end

%      KB =Sdraft(i)/2;
%      BM =Slength(i)^3/12/Slength(i)/Sdraft(i);
%      KM =KB+BM;
%      KG =KM-GM(i);
%      CoG=KG-Sdraft(i)
%      BC =Sdraft(i)/2;
     Xrot0(i)= 0;
     Zrot0(i)= ZCoG(i);
    
     try
        CdX(i)=str2num(cell2mat(shippar.data(i,15)));
    catch
        CdX(i)=cell2mat(shippar.data(i,15));
    end
    try
        CdZ(i)=str2num(cell2mat(shippar.data(i,16)));
    catch
        CdZ(i)=cell2mat(shippar.data(i,16));
    end
    try
        CdTheta(i)=str2num(cell2mat(shippar.data(i,17)));
    catch
        CdTheta(i)=cell2mat(shippar.data(i,17));
    end
    
    SXc=Xc(i);
    SZc=Zc(i);
    Sthetac=Thetac(i);
    [shapeXcZcTheta0I,~,~]=shapeformLib(x,SXc,SZc,0,Sdraft(i),Slength(i)/2,[],shipform);
     
    shapeXcZ0I=zeros(size(shapeXcZcTheta0I));
    shapeXcZ0I(shapeXcZcTheta0I<0)=shapeXcZcTheta0I(shapeXcZcTheta0I<0)-SZc;
   
    [shapeXcZcI,dxshapeXcZcO,shapeXcZcdati]=shapeformLib(x,SXc,SZc,Sthetac,Sdraft(i),Slength(i)/2,[],shipform);
   
    shapeXcZc(i,:)=shapeXcZcI.';
    dxshapeXcZc(i,:)= dxshapeXcZcO.';%% for barge: delta function outside for calculating fluid momentum
    shapeXcZ0(i,:)=shapeXcZ0I.';
    Chi(i,:)    = heaviside(-sign(shapeXcZcI)-10^(-12));
  
    SChi        = SChi+Chi(i,:);
    indxL       = find(Chi(i,:),1,'first')-1;
    xShip(i,1)  = x(indxL);
    indxR       = find(Chi(i,:),1,'last')+1;
    xShip(i,3)  = x(indxR);
    indxC       = indxL+floor((indxR-indxL)/2);
    xShip(i,2)  = x(indxC);
    Slength(i)  = x(indxR)-x(indxL)-dx; %for barge
  
    dxChi(i,indxL)= 1/dx;
    dxChi(i,indxR)= -1/dx;
    bf0           = Slength(i)/3;
    cfSA_ship     = cfSA_ship+cf(x,x(indxL+1),bf0)-cf(x,x(indxR-1)-bf0,bf0); % cfSimulationArea
    
     
    SdxChi        = SdxChi+dxChi(i,:);
    
    ChiWl(i,indxL)=1;
    ChiWl(i,indxR)=1;
    SChiWl        = SChiWl+ChiWl(i,:);
    ChiSl(i,indxL+1)=1;
    ChiSl(i,indxR-1)=1;
    SChiSl        = SChiSl+ChiSl(i,:);
    
    Mass(i)     = -trapz(x,min(0,shapeXcZ0I.'));
    LWL(i)      = trapz(x,heaviside(-sign(shapeXcZ0I.')-10^-12));
    if strcmpi(shippar.data(i,2),'Pitch')||strcmpi(shippar.data(i,2),'Free')         
    MIner(i)=Mass(i) .*Gyradius(i).^2;
    end
    SshapeXcZc      = SshapeXcZc+shapeXcZc(i,:);
    SdxshapeXcZc    = SdxshapeXcZc+dxshapeXcZc(i,:);
    SshapeXcZ0      = SshapeXcZ0+shapeXcZ0I(i,:);
    shapeN=shapeXcZc(i,:);%shapeN([indxL indxR])=shapeXcZc(i,[indxL+1 indxR-1]);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmpi(shipform,'barge') %&& strcmpi(bathName,'Flat')
        Ndom=6;
    else
        NdxS=(indxR-indxL)+1;
        Ndom=(NdxS-1)/2+2; % 10; % number of domain 16
        if mod(Ndom,2)~=0
            Ndom=Ndom+1;
        end
    end
    
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
% figure;
% plot(x,nuthetatild)
Chi(i+1,:)=SChi;
dxChi(i+1,:)=SdxChi;
ChiWl(i+1,:)=SChiWl;
ChiSl(i+1,:)=SChiSl;
shapeXcZc(i+1,:)=SshapeXcZc;
dxshapeXcZc(i+1,:)=SdxshapeXcZc;
shapeXcZ0(i+1,:)=SshapeXcZ0;
nu.z(i+1,:)=nuz;
nu.x(i+1,:)=nux;
nu.theta(i+1,:)  =nutheta;
nutild.z(i+1,:)=nuztild;
nutild.x(i+1,:)=nuxtild;
nutild.theta(i+1,:)=nuthetatild;

parS.Cd.X=CdX;
parS.Cd.Z=CdZ;
parS.Cd.Theta=CdTheta;

parS.GM=GM;
parS.shapedat=shapedat;
parS.shapeXcZc=shapeXcZc;parS.dxshapeXcZc=dxshapeXcZc;
parS.shapeXcZ0=shapeXcZ0;
parS.Mass=Mass;parS.LWL=LWL;
parS.MIner=MIner;
for ii=1:Nship
parS.MassS(ii).Mat=diag([Mass(ii) Mass(ii) MIner(ii)]);
end
parS.cfSA_ship=cfSA_ship;
parS.chi=Chi';parS.xShip=xShip;parS.dxchi=dxChi';
parS.chiWl=ChiWl';parS.chiSl=ChiSl';
parS.nu=nu;parS.nutild=nutild;
parS.Slength=Slength;
parS.Sdraft=Sdraft;
parS.XZrot0=[Xrot0 Zrot0];
parS.Ndom=Ndom;

% figure;
% plot(x,parS.dxchi)

% figure
% plot(x,SshapeXcZc,'b',x,(1-parS.chi(:,end))+parS.chiSl(:,end),'k--',x,SdxshapeXcZc,'r-.');
% title(['Mass (total)=',num2str(sum(Mass)),'; LWL (total)= ',num2str(sum(LWL))])
% plot_properties;

