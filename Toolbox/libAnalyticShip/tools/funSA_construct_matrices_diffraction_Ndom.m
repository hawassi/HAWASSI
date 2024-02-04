%% function for constructing matrix system for solving radiation
%  17 December 2019
%  applying matching condition at ship boundary;
%  domain composition:
%  layoclcut: % ....(-(b+B))------(-b)....0....(b)------(b+B)....
%  domain: %   1           2           3          4         5

%Amat.*Xvect=Bvect
%Xvect=[A1 A3L A3R A5 tau0 tau1 tau2 tau3 eps1m eps21m eps22m eps31m eps32m eps41m eps42m eps5m]
function [Amat,Bvect] =funSA_construct_matrices_diffraction_Ndom(Ndom,psibdy,dxpsibdy,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,Evmodes,kappaS,etaWl,dxetaWl,Xbdy,A)

Ns=1;
Nm=2*Ndom-2; % number of unknown for propagation term only;
N=Nm*(Evmodes+1);

Nbdy=length(Dbdy);
eta1L=etaWl(1);eta1R=etaWl(2);
dxeta1L=dxetaWl(1);dxeta1R=dxetaWl(2);
X1L=Xbdy(1);X1R=Xbdy(end);

Amat=zeros(N,N);
Bvect=zeros(N,1);

for n=0:Evmodes
    muS1L_n=n*pi/(Dbdy(1)-dr(1)); %% mu for 1st ship
    muS1R_n=n*pi/(Dbdy(end)-dr(end));
    
    IVp0Vpn1L=fun_IntegVarphimVarphin(0,n,Dbdy(1),-dr(1));
    IVp0Vpn1R=fun_IntegVarphimVarphin(0,n,Dbdy(end),-dr(end));
    
    
    NvarI=2*(Nbdy-3)+1;
    muSIn_n =zeros(NvarI,1);
    IdC=(NvarI-1)/2+1;IdbdyC=(Nbdy-1)/2+1;
    muSIn_n(IdC)=n*pi/(Dbdy(IdbdyC)-dr(IdbdyC));
    Idbdy=2;
    
    IVp0VpnIn=zeros(NvarI,1);
    IVp0VpnIn(IdC)=fun_IntegVarphimVarphin(0,n,Dbdy(IdbdyC),-dr(IdbdyC));
    IQnVp0In=zeros(NvarI,1);
    IQnVp0In(IdC)=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,0,Dbdy(IdbdyC),-dr(IdbdyC));
    IdbdyC =IdbdyC+1;
    
    
    for ii=1:(NvarI-1)/2
        if mod(ii,2)~=0
            muSIn_n(ii)=n*pi/(Dbdy(Idbdy)-dr(Idbdy-1));
            muSIn_n(IdC+ii)=n*pi/(Dbdy(IdbdyC)-dr(IdbdyC));
            IVp0VpnIn(ii)=fun_IntegVarphimVarphin(0,n,Dbdy(Idbdy),-dr(Idbdy-1));
            IVp0VpnIn(IdC+ii)=fun_IntegVarphimVarphin(0,n,Dbdy(IdbdyC),-dr(IdbdyC));
            IQnVp0In(ii)=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,Dbdy(Idbdy),-dr(Idbdy-1));
            IQnVp0In(IdC+ii)=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,Dbdy(IdbdyC),-dr(IdbdyC));
        else
            muSIn_n(ii)=n*pi/(Dbdy(Idbdy)-dr(Idbdy));
            muSIn_n(IdC+ii)=n*pi/(Dbdy(IdbdyC)-dr(IdbdyC+1));
            IVp0VpnIn(ii)=fun_IntegVarphimVarphin(0,n,Dbdy(Idbdy),-dr(Idbdy));
            IVp0VpnIn(IdC+ii)=fun_IntegVarphimVarphin(0,n,Dbdy(IdbdyC),-dr(IdbdyC+1));
            IQnVp0In(ii)=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,Dbdy(Idbdy),-dr(Idbdy));
            IQnVp0In(IdC+ii)=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,Dbdy(IdbdyC),-dr(IdbdyC+1));
    
            Idbdy=Idbdy+1;IdbdyC=IdbdyC+1;
        end
        
    end
    
    IQ0Vpn1L=fun_IntegQmVarphin(0,n,kappaS(1,1),muS1L_n,eta1L,Dbdy(1),-dr(1));
    IQ0Vpn1R=fun_IntegQmVarphin(0,n,kappaS(1,2),muS1R_n,eta1R,Dbdy(end),-dr(end));
    
    IQnVp01L=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,Dbdy(1),-dr(1));
    IQnVp01R=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,Dbdy(end),-dr(end));
    
    IQ0Qn1L=fun_IntegQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dbdy(1),eta1L);
    IQ0Qn1R=fun_IntegQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dbdy(end),eta1R);
    
    IdxQ0Qn1L=fun_IntegdxQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dbdy(1),dxeta1L,dxDbdy(1),eta1L);
    IdxQ0Qn1R=fun_IntegdxQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dbdy(end),dxeta1R,dxDbdy(end),eta1R);
    
    Amat(Nm*Ns*n+1,1:3)        =[ IQ0Vpn1L  , -IVp0Vpn1L , -IVp0Vpn1L*tauGbdy(1)];
    Amat(Nm*Ns*n+Nm/2,Nm-2:Nm) =[ IVp0Vpn1R , IVp0Vpn1R*tauGbdy(end)      , -IQ0Vpn1R  ];
    Amat(Nm*Ns*n+Nm/2+1,1:3)   =[ -kappaS(1,1)*IQ0Qn1L+IdxQ0Qn1L , 0          , -IQnVp01L*dxtauGbdy(1)];
    Amat(Nm*Ns*n+Nm,Nm-2:Nm)   =[ 0           ,  IQnVp01R*dxtauGbdy(end)     , -(kappaS(1,2)*IQ0Qn1R+IdxQ0Qn1R) ];
   
    idcol=2;IdIn=1;
    
    for ii=2:Nm/2-1
        if ii<(Nm/2-1)/2+1
             Amat(Nm*Ns*n+ii,idcol:idcol+3) =[IVp0VpnIn(IdIn+1), IVp0VpnIn(IdIn+1)*tauGbdy(ii), -IVp0VpnIn(IdIn+1)    , -IVp0VpnIn(IdIn+1)*tauGbdy(ii)];
             Amat(Nm*Ns*n+Nm/2+ii,idcol:idcol+3) =[ 0   , IQnVp0In(IdIn)*dxtauGbdy(ii),  0 , -IQnVp0In(IdIn+1)*dxtauGbdy(ii)];
      
        elseif ii>(Nm/2-1)/2+1
             Amat(Nm*Ns*n+ii,idcol:idcol+3) =[IVp0VpnIn(IdIn-1), IVp0VpnIn(IdIn-1)*tauGbdy(ii), -IVp0VpnIn(IdIn-1)    , -IVp0VpnIn(IdIn-1)*tauGbdy(ii)];
             Amat(Nm*Ns*n+Nm/2+ii,idcol:idcol+3) =[ 0   , IQnVp0In(IdIn-1)*dxtauGbdy(ii),  0 , -IQnVp0In(IdIn)*dxtauGbdy(ii)];
    
        else
             Amat(Nm*Ns*n+ii,idcol:idcol+3)=[IVp0VpnIn(IdIn), IVp0VpnIn(IdIn)*tauGbdy(ii), -IVp0VpnIn(IdIn)    , -IVp0VpnIn(IdIn)*tauGbdy(ii)];
             Amat(Nm*Ns*n+Nm/2+ii,idcol:idcol+3) =[ 0   , IQnVp0In(IdIn)*dxtauGbdy(ii),  0 , -IQnVp0In(IdIn)*dxtauGbdy(ii)];
      
        end
        idcol=idcol+2;IdIn=IdIn+2;
    end
%     assignin('base','AmatI',AmatI);
%     pause;
    %Xvect=[Al                            sig20         sig21                     sig30         sig31                         sig40          sig41                       sig50         sig51                        Ar ]
    
      
    
    if Evmodes>0
        for m=1:Evmodes
            
            NvarI=2*(Nbdy-3)+1;
            muSIn_m =zeros(NvarI,1);
            dxmuS1In_m=zeros(NvarI,1);
            IdC=(NvarI-1)/2+1;IdbdyC=(Nbdy-1)/2+1;
            muSIn_m(IdC)=m*pi/(Dbdy(IdbdyC)-dr(IdbdyC));
           % dxmuS1In_m(IdC)=-m*pi*dxDZetabdy(IdbdyC)./(DZetabdy(IdbdyC)^2);
            Idbdy=2;
            
            IVpmVpnIn=zeros(NvarI,1);
            IVpmVpnIn(IdC)=fun_IntegVarphimVarphin(m,n,Dbdy(IdbdyC),-dr(IdbdyC));
            
            IQnVpmIn=zeros(NvarI,1);
            IQnVpmIn(IdC)=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muSIn_m(IdC),0,Dbdy(IdbdyC),-dr(IdbdyC));
            IdxVpmQnIn=zeros(NvarI,1);
            IdxVpmQnIn(IdC)=fun_IntegdxVarphimQn(muSIn_m(IdC),kappaS(n+1,1),0,Dbdy(IdbdyC),dxDbdy(IdbdyC),-dr(IdbdyC));
            IdbdyC=IdbdyC+1;
            for ii=1:(NvarI-1)/2
                if mod(ii,2)~=0
                    muSIn_m(ii)=m*pi/(Dbdy(Idbdy)-dr(Idbdy-1));
                    muSIn_m(IdC+ii)=m*pi/(Dbdy(IdbdyC)-dr(IdbdyC));
                    IVpmVpnIn(ii)=fun_IntegVarphimVarphin(m,n,Dbdy(Idbdy),-dr(Idbdy-1));
                    IVpmVpnIn(IdC+ii)=fun_IntegVarphimVarphin(m,n,Dbdy(IdbdyC),-dr(IdbdyC));
                    IQnVpmIn(ii)=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muSIn_m(ii),eta1L,Dbdy(Idbdy),-dr(Idbdy-1));
                    IQnVpmIn(IdC+ii)=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muSIn_m(IdC+ii),eta1R,Dbdy(IdbdyC),-dr(IdbdyC));
                    IdxVpmQnIn(ii)=fun_IntegdxVarphimQn(muSIn_m(ii),kappaS(n+1,1),eta1L,Dbdy(Idbdy),dxDbdy(Idbdy),-dr(Idbdy-1));
                    IdxVpmQnIn(IdC+ii)=fun_IntegdxVarphimQn(muSIn_m(IdC+ii),kappaS(n+1,2),eta1R,Dbdy(IdbdyC),dxDbdy(IdbdyC),-dr(IdbdyC));
                else
                    muSIn_m(ii)=m*pi/(Dbdy(Idbdy)-dr(Idbdy));
                    muSIn_m(IdC+ii)=m*pi/(Dbdy(IdbdyC)-dr(IdbdyC+1));
                    IVpmVpnIn(ii)=fun_IntegVarphimVarphin(m,n,Dbdy(Idbdy),-dr(Idbdy));
                    IVpmVpnIn(IdC+ii)=fun_IntegVarphimVarphin(m,n,Dbdy(IdbdyC),-dr(IdbdyC+1));
                    IQnVpmIn(ii)=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muSIn_m(ii),eta1L,Dbdy(Idbdy),-dr(Idbdy));
                    IQnVpmIn(IdC+ii)=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muSIn_m(IdC+ii),eta1R,Dbdy(IdbdyC),-dr(IdbdyC+1));
                     IdxVpmQnIn(ii)=fun_IntegdxVarphimQn(muSIn_m(ii),kappaS(n+1,1),eta1L,Dbdy(Idbdy),dxDbdy(Idbdy),-dr(Idbdy));
                     IdxVpmQnIn(IdC+ii)=fun_IntegdxVarphimQn(muSIn_m(IdC+ii),kappaS(n+1,1),eta1R,Dbdy(IdbdyC),dxDbdy(IdbdyC),-dr(IdbdyC+1));
                    
                    Idbdy=Idbdy+1;IdbdyC=IdbdyC+1;
                end
                
            end
            
            muS1L_m=m*pi/(Dbdy(1)-dr(1));
            muS1R_m=m*pi/(Dbdy(end)-dr(end));
            
            dxmuS1L_m=-m*pi*dxDZetabdy(1)./(DZetabdy(1)^2);
            dxmuS1R_m=-m*pi*dxDZetabdy(end)./(DZetabdy(end)^2);
            
            IVpmVpn1L=fun_IntegVarphimVarphin(m,n,Dbdy(1),-dr(1));
            IVpmVpn1R=fun_IntegVarphimVarphin(m,n,Dbdy(end),-dr(end));
            
            IQmVpn1L=fun_IntegQmVarphin(m,n,kappaS(m+1,1),muS1L_n,eta1L,Dbdy(1),-dr(1));
            IQmVpn1R=fun_IntegQmVarphin(m,n,kappaS(m+1,2),muS1R_n,eta1R,Dbdy(end),-dr(end));
            
            IQnVpm1L=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muS1L_m,eta1L,Dbdy(1),-dr(1));
            IQnVpm1R=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muS1R_m,eta1R,Dbdy(end),-dr(end));
            
            IQmQn1L=fun_IntegQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dbdy(1),eta1L);
            IQmQn1R=fun_IntegQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dbdy(end),eta1R);
            
            IdxQmQn1L=fun_IntegdxQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dbdy(1),dxeta1L,dxDbdy(1),eta1L);
            IdxQmQn1R=fun_IntegdxQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dbdy(end),dxeta1R,dxDbdy(end),eta1R);
            
            IdxVpmQn1L  =fun_IntegdxVarphimQn(muS1L_m,kappaS(n+1,1),eta1L,Dbdy(1),dxDbdy(1),-dr(1));
            IdxVpmQn1R=fun_IntegdxVarphimQn(muS1R_m,kappaS(n+1,2),eta1R,Dbdy(end),dxDbdy(end),-dr(end));
            
            Amat(Nm*Ns*n+1,Nm*Ns*m+(1:3)) =[IQmVpn1L , -IVpmVpn1L, -exp(-muS1L_m*(Xbdy(1)-Xbdy(2)))*IVpmVpn1L ];
            Amat(Nm*Ns*n+Nm/2,Nm*Ns*(m+1)-2:Nm*Ns*(m+1))  =[exp(muS1R_m*(Xbdy(end)-Xbdy(end-1)))*IVpmVpn1R, IVpmVpn1R,-IQmVpn1R ];
            Amat(Nm*Ns*n+Nm/2+1,Nm*Ns*m+(1:3))=[kappaS(m+1,1)*IQmQn1L+IdxQmQn1L , -(muS1L_m*IQnVpm1L+IdxVpmQn1L)   , -exp(-muS1L_m*(Xbdy(1)-Xbdy(2)))*((-muS1L_m-dxmuS1L_m*(Xbdy(1)-Xbdy(2)))*IQnVpm1L+IdxVpmQn1L)];
            Amat(Nm*Ns*n+Nm,Nm*Ns*(m+1)-2:Nm*Ns*(m+1))   =[ exp(muS1R_m*(Xbdy(end)-Xbdy(end-1)))*((muS1R_m+dxmuS1R_m*(Xbdy(end)-Xbdy(end-1)))*IQnVpm1R+IdxVpmQn1R),(-muS1R_m*IQnVpm1R+IdxVpmQn1R),kappaS(m+1,2)*IQmQn1R-IdxQmQn1R  ];
                                                           
            idcol=2;IdIn=1;
            for ii=2:Nm/2-1
               % ii
                if ii<(Nm/2-1)/2+1
                    Amat(Nm*Ns*n+ii,Nm*Ns*m+(idcol:idcol+3)) =[exp(muSIn_m(IdIn+1)*(Xbdy(ii)-Xbdy(ii-1)))*IVpmVpnIn(IdIn+1), IVpmVpnIn(IdIn+1), -IVpmVpnIn(IdIn+1)    , -exp(-muSIn_m(IdIn+1)*(Xbdy(ii)-Xbdy(ii+1)))*IVpmVpnIn(IdIn+1)];
                    Amat(Nm*Ns*n+Nm/2+ii,Nm*Ns*m+(idcol:idcol+3)) =[ exp(muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii-1)))*((muSIn_m(IdIn)+dxmuS1In_m(IdIn)*(Xbdy(ii)-Xbdy(ii-1)))*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn)),  -muSIn_m(IdIn)*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn) , -muSIn_m(IdIn+1)*IQnVpmIn(IdIn+1)-IdxVpmQnIn(IdIn+1), -exp(-muSIn_m(IdIn+1)*(Xbdy(ii)-Xbdy(ii+1)))*((-muSIn_m(IdIn+1)-dxmuS1In_m(IdIn+1)*(Xbdy(ii)-Xbdy(ii+1)))*IQnVpmIn(IdIn+1)+IdxVpmQnIn(IdIn+1))];
                elseif ii>(Nm/2-1)/2+1
                    Amat(Nm*Ns*n+ii,Nm*Ns*m+(idcol:idcol+3)) =[exp(muSIn_m(IdIn-1)*(Xbdy(ii)-Xbdy(ii-1)))*IVpmVpnIn(IdIn-1), IVpmVpnIn(IdIn-1), -IVpmVpnIn(IdIn-1)    , -exp(-muSIn_m(IdIn-1)*(Xbdy(ii)-Xbdy(ii+1)))*IVpmVpnIn(IdIn-1)];
                    Amat(Nm*Ns*n+Nm/2+ii,Nm*Ns*m+(idcol:idcol+3)) =[ exp(muSIn_m(IdIn-1)*(Xbdy(ii)-Xbdy(ii-1)))*((muSIn_m(IdIn-1)+dxmuS1In_m(IdIn-1)*(Xbdy(ii)-Xbdy(ii-1)))*IQnVpmIn(IdIn-1)+IdxVpmQnIn(IdIn-1)),  -muSIn_m(IdIn-1)*IQnVpmIn(IdIn-1)+IdxVpmQnIn(IdIn-1) , -muSIn_m(IdIn)*IQnVpmIn(IdIn)-IdxVpmQnIn(IdIn), -exp(-muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii+1)))*((-muSIn_m(IdIn)-dxmuS1In_m(IdIn)*(Xbdy(ii)-Xbdy(ii+1)))*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn))];
                else
                    Amat(Nm*Ns*n+ii,Nm*Ns*m+(idcol:idcol+3)) =[exp(muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii-1)))*IVpmVpnIn(IdIn), IVpmVpnIn(IdIn), -IVpmVpnIn(IdIn)    , -exp(-muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii+1)))*IVpmVpnIn(IdIn)];
                    Amat(Nm*Ns*n+Nm/2+ii,Nm*Ns*m+(idcol:idcol+3)) =[ exp(muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii-1)))*((muSIn_m(IdIn)+dxmuS1In_m(IdIn)*(Xbdy(ii)-Xbdy(ii-1)))*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn)),  -muSIn_m(IdIn)*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn) , -muSIn_m(IdIn)*IQnVpmIn(IdIn+1)-IdxVpmQnIn(IdIn), -exp(-muSIn_m(IdIn)*(Xbdy(ii)-Xbdy(ii+1)))*((-muSIn_m(IdIn)-dxmuS1In_m(IdIn)*(Xbdy(ii)-Xbdy(ii+1)))*IQnVpmIn(IdIn)+IdxVpmQnIn(IdIn))];
                end
                idcol=idcol+2;IdIn=IdIn+2;
            end
         %   pause;
        end
    end
    
    DT1L  =Dbdy(1)-dr(1);
    DT1R  =Dbdy(end)-dr(end);
    
    
    IZD2varphinDT1L = fun_IntegZD2Varphin(muS1L_n,Dbdy(1),-Dbdy(1),-dr(1));
    IZD2varphinDT1R = fun_IntegZD2Varphin(muS1R_n,Dbdy(end),-Dbdy(end),-dr(end));
    
    
    IntQnDT1In=zeros(NvarI,1);
    IZDQnDT1In=zeros(NvarI,1);
    IZD2QnDT1In=zeros(NvarI,1);
    
    DT1In=zeros(NvarI,1);
    Idbdy=2;
    IdC=(NvarI-1)/2+1;IdbdyC=(Nbdy-1)/2+1+1;
    for ii=1:(NvarI-1)/2
        if mod(ii,2)~=0
            DT1In(ii)=Dbdy(Idbdy)-dr(Idbdy-1);
            DT1In(IdC+ii)=Dbdy(IdbdyC)-dr(IdbdyC);
            
            IntQnDT1In(ii)=fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr( Idbdy-1));
            IntQnDT1In(IdC+ii)=fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC));
            IZDQnDT1In(ii)=fun_IntegZDQn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr(Idbdy-1));
            IZDQnDT1In(IdC+ii)=fun_IntegZDQn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC));
            IZD2QnDT1In(ii)=fun_IntegZD2Qn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr(Idbdy-1));
            IZD2QnDT1In(IdC+ii)=fun_IntegZD2Qn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC));
            
        else
            DT1In(ii)=Dbdy(Idbdy)-dr(Idbdy);
            DT1In(IdC+ii)=Dbdy(IdbdyC)-dr(IdbdyC+1);
            
            IntQnDT1In(ii)=fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr( Idbdy));
            IntQnDT1In(IdC+ii)=fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC+1));
            IZDQnDT1In(ii)=fun_IntegZDQn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr(Idbdy));
            IZDQnDT1In(IdC+ii)=fun_IntegZDQn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC+1));
            IZD2QnDT1In(ii)=fun_IntegZD2Qn(kappaS(n+1,1),eta1L,Dbdy(Idbdy),-Dbdy(Idbdy),-dr(Idbdy));
            IZD2QnDT1In(IdC+ii)=fun_IntegZD2Qn(kappaS(n+1,2),eta1R,Dbdy(IdbdyC),-Dbdy(IdbdyC),-dr(IdbdyC+1));
            
            
            Idbdy=Idbdy+1;IdbdyC=IdbdyC+1;
        end
       
    end
       Bvect(Nm*Ns*n+1)=-A*IQ0Vpn1L;
       Bvect(Nm*Ns*n+Nm/2+1)=-A*(kappaS(1,1)*IQ0Qn1L+IdxQ0Qn1L);
     
%      assignin('base','AmatI',AmatI);
%      assignin('base','BvectHI',BvectHI);
%      pause;
    
    
end