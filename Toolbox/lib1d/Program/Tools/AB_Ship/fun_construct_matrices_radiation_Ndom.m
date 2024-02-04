%% function for constructing matrix system for solving radiation
%  17 December 2019
%  applying matching condition at ship boundary;
%  domain composition:
%  layoclcut: % ....(-(b+B))------(-b)....0....(b)------(b+B)....
%  domain: %   1           2           3          4         5

%Amat.*Xvect=Bvect
%Xvect=[A1 A3L A3R A5 tau0 tau1 tau2 tau3 eps1m eps21m eps22m eps31m eps32m eps41m eps42m eps5m]
function [Amat,BvectS,BvectH,BvectP] =fun_construct_matrices_radiation_Ndom(Ndom,psibdy,dxpsibdy,tauGbdy,dxtauGbdy,Dbdy,dxDbdy,DZetabdy,dxDZetabdy,dr,Evmodes,kappaS,Idmotion,X0,Z0,etaWl,dxetaWl,Xbdy)
Ns=1;
Nm=2*Ndom-2; % number of unknown for propagation term only;
N=Nm*(Evmodes+1);

Nbdy=length(Dbdy);
eta1L=etaWl(1);eta1R=etaWl(2);
dxeta1L=dxetaWl(1);dxeta1R=dxetaWl(2);
X1L=Xbdy(1);X1R=Xbdy(end);

Amat=zeros(N,N);
BvectH=zeros(N,1);
BvectP=zeros(N,1);
BvectS=zeros(N,1);

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
    
    
    
    IvarphinDT1L = fun_IntegVarphin(muS1L_n,Dbdy(1),-Dbdy(1),-dr(1));
    IvarphinDT1R = fun_IntegVarphin(muS1R_n,Dbdy(end),-Dbdy(end),-dr(end));
    
    IntQnDT1L=fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(1),-Dbdy(1),-dr(1));
    IntQnDT1R=fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(end),-Dbdy(end),-dr(end));
    
    IZDQnDT1L = fun_IntegZDQn(kappaS(n+1,1),eta1L,Dbdy(1),-Dbdy(1),-dr(1));
    IZDQnDT1R = fun_IntegZDQn(kappaS(n+1,2),eta1R,Dbdy(end),-Dbdy(end),-dr(end));
    
    
    IZD2QnDT1L = fun_IntegZD2Qn(kappaS(n+1,1),eta1L,Dbdy(1),-Dbdy(1),-dr(1));
    IZD2QnDT1R = fun_IntegZD2Qn(kappaS(n+1,2),eta1R,Dbdy(end),-Dbdy(end),-dr(end));
    
    if strcmpi(Idmotion,'Surge') || strcmpi(Idmotion,'Free') % heave
       
       InQnTEtaL=fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
       InQnTEtaR=fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
    
        BvectS(Nm*Ns*n+Nm/2+1)=InQnTEtaL;
        BvectS(Nm*Ns*n+Nm)=-InQnTEtaR;
            
       IdIn=1;
        for ii=2:Nm/2-1
            if ii<(Nm/2-1)/2+1
                BvectS(Nm*Ns*n+Nm/2+ii)=fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(ii),-dr(ii),-dr(ii-1));
            elseif ii>(Nm/2-1)/2+1
                BvectS(Nm*Ns*n+Nm/2+ii)=-fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(ii),-dr(ii),-dr(ii+1));
            end
            IdIn=IdIn+2;
        end
    end
 
    
    if strcmpi(Idmotion,'Heave') || strcmpi(Idmotion,'Free') % heave
        InGamma31LVarphin=(IZD2varphinDT1L-(X1L-X0(1))^2*IvarphinDT1L)/2/DT1L;
        InGamma31RVarphin=(IZD2varphinDT1R-(X1R-X0(1))^2*IvarphinDT1R)/2/DT1R;
        
        IndxGamma31LQn  =-((X1L-X0(1))./DT1L- (X1L-X0(1)).^2*dxDZetabdy(1)./2./DT1L/DT1L)*IntQnDT1L ....
            +dxDbdy(1)*IZDQnDT1L./DT1L-dxDZetabdy(1)*IZD2QnDT1L/2/DT1L/DT1L;
        IndxGamma31RQn  =-((X1R-X0(1))./DT1R -(X1R-X0(1)).^2*dxDZetabdy(end)./2./DT1R/DT1R)*IntQnDT1R ....
            +dxDbdy(end)*IZDQnDT1R./DT1R-dxDZetabdy(end)*IZD2QnDT1R/2/DT1R/DT1R;
        %
       
        Idbdy=2;
        IdC=(NvarI-1)/2+1;IdbdyC=(Nbdy-1)/2+1;
        IndxGamma31InQn=zeros(NvarI,1);
 
        for ii=1:(NvarI-1)/2
                       
            
            IndxGamma31InQn(ii)= -((Xbdy(Idbdy)-X0(1))./DT1In(ii)-(Xbdy(Idbdy)-X0(1)).^2*dxDZetabdy(Idbdy)./2./DT1In(ii)/DT1In(ii))*IntQnDT1In(ii) ...
                +dxDbdy(Idbdy)*IZDQnDT1In(ii)./DT1In(ii)-dxDZetabdy(Idbdy)*IZD2QnDT1In(ii)/2/DT1In(ii)/DT1In(ii);
            
            IndxGamma31InQn(IdC+ii)= -((Xbdy(IdbdyC+1)-X0(1))./DT1In(IdC+ii)-(Xbdy(IdbdyC+1)-X0(1)).^2*dxDZetabdy(IdbdyC+1)./2./DT1In(IdC+ii)/DT1In(ii))*IntQnDT1In(IdC+ii) ...
                +dxDbdy(IdbdyC+1)*IZDQnDT1In(IdC+ii)./DT1In(IdC+ii)-dxDZetabdy(IdbdyC+1)*IZD2QnDT1In(IdC+ii)/2/DT1In(IdC+ii)/DT1In(IdC+ii);
            
            if mod(ii,2)==0
                Idbdy=Idbdy+1;IdbdyC=IdbdyC+1;
            end
            
        end
      
        BvectH(Nm*Ns*n+1)=InGamma31LVarphin;
        BvectH(Nm*Ns*n+Nm/2)=-InGamma31RVarphin;
        BvectH(Nm*Ns*n+Nm/2+1)=IndxGamma31LQn;
        BvectH(Nm*Ns*n+Nm)=-IndxGamma31RQn;
        
       IdIn=1;
        for ii=2:Nm/2-1
            if ii<(Nm/2-1)/2+1
               BvectH(Nm*Ns*n+Nm/2+ii)=IndxGamma31InQn(IdIn+1)-IndxGamma31InQn(IdIn);
            elseif ii>(Nm/2-1)/2+1
               BvectH(Nm*Ns*n+Nm/2+ii)=IndxGamma31InQn(IdIn)-IndxGamma31InQn(IdIn-1);
            end
            IdIn=IdIn+2;
        end
        
    end
    
   if strcmpi(Idmotion,'Pitch') || strcmpi(Idmotion,'Free') 
        InGamma51LVarphin=(IZD2varphinDT1L*(X1L-X0(1))-(X1L-X0(1))^3/3*IvarphinDT1L)/2/DT1L;
        InGamma51RVarphin=(IZD2varphinDT1R*(X1R-X0(1))-(X1R-X0(1))^3/3*IvarphinDT1R)/2/DT1R;
        
        IndxGamma51LQn  =(IZD2QnDT1L+2*IZDQnDT1L*dxDbdy(1)*(X1L-X0(1))-(X1L-X0(1))^2*IntQnDT1L)/2/DT1L-...
                         (IZD2QnDT1L*(X1L-X0(1))-(X1L-X0(1))^3/3*IntQnDT1L)*dxDZetabdy(1)/2/DT1L/DT1L; 
        IndxGamma51RQn  =(IZD2QnDT1R+2*IZDQnDT1R*dxDbdy(end)*(X1R-X0(1))-(X1R-X0(1))^2*IntQnDT1R)/2/DT1R-...
                         (IZD2QnDT1R*(X1R-X0(1))-(X1R-X0(1))^3/3*IntQnDT1R)*dxDZetabdy(end)/2/DT1R/DT1R;
        Idbdy=2;
        IdC=(NvarI-1)/2+1;IdbdyC=(Nbdy-1)/2+1;
        IndxGamma51InQn=zeros(NvarI,1);
        
        for ii=1:(NvarI-1)/2                
                 IndxGamma51InQn(ii)=(IZD2QnDT1In(ii)+2*IZDQnDT1In(ii)*dxDbdy(Idbdy)*(Xbdy(Idbdy)-X0(1))-(Xbdy(Idbdy)-X0(1))^2*IntQnDT1In(ii))/2/DT1In(ii)-...
                         (IZD2QnDT1In(ii)*(Xbdy(Idbdy)-X0(1))-(Xbdy(Idbdy)-X0(1))^3/3*IntQnDT1In(ii))*dxDZetabdy(Idbdy)/2/DT1In(ii)/DT1In(ii);  
                 IndxGamma51InQn(IdC+ii)=(IZD2QnDT1In(IdC+ii)+2*IZDQnDT1In(IdC+ii)*dxDbdy(IdbdyC+1)*(Xbdy(IdbdyC+1)-X0(1))-(Xbdy(IdbdyC+1)-X0(1))^2*IntQnDT1In(IdC+ii))/2/DT1In(IdC+ii)-...
                         (IZD2QnDT1In(IdC+ii)*(Xbdy(IdbdyC+1)-X0(1))-(Xbdy(IdbdyC+1)-X0(1))^3/3*IntQnDT1In(IdC+ii))*dxDZetabdy(IdbdyC+1)/2/DT1In(IdC+ii)/DT1In(IdC+ii);  
                       
            if mod(ii,2)==0
               Idbdy=Idbdy+1;IdbdyC=IdbdyC+1;
            end
        end
        
        IZQnTeta1L = fun_IntegZQn(kappaS(n+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
        IZQnTeta1R = fun_IntegZQn(kappaS(n+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
        IQnTeta1L  = fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(1),-dr(1),eta1L);
        IQnTeta1R  = fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(end),-dr(end),eta1R);
       
        Intdxpsi5QnTeta1L=-IZQnTeta1L+Z0(1).*IQnTeta1L;
        Intdxpsi5QnTeta1R=-IZQnTeta1R+Z0(1).*IQnTeta1R;
        
      
        BvectP(Nm*Ns*n+1)=InGamma51LVarphin;
        BvectP(Nm*Ns*n+Nm/2)=-InGamma51RVarphin;
        BvectP(Nm*Ns*n+Nm/2+1)=Intdxpsi5QnTeta1L+IndxGamma51LQn;
        BvectP(Nm*Ns*n+Nm)=-Intdxpsi5QnTeta1R-IndxGamma51RQn;
        
       IdIn=1;
        for ii=2:Nm/2-1
            if ii<(Nm/2-1)/2+1
               IZQnTetaIn = fun_IntegZQn(kappaS(n+1,1),eta1L,Dbdy(ii),-dr(ii),-dr(ii-1));
               IQnTetaIn = fun_IntegQn(kappaS(n+1,1),eta1L,Dbdy(ii),-dr(ii),-dr(ii-1));
               Intdxpsi5QnTetaIn=-IZQnTetaIn+Z0(1).*IQnTetaIn;
               BvectP(Nm*Ns*n+Nm/2+ii)=Intdxpsi5QnTetaIn*0+IndxGamma51InQn(IdIn+1)-IndxGamma51InQn(IdIn);
               
            elseif ii>(Nm/2-1)/2+1
                IZQnTetaIn = fun_IntegZQn(kappaS(n+1,2),eta1R,Dbdy(ii),-dr(ii),-dr(ii+1));
                IQnTetaIn = fun_IntegQn(kappaS(n+1,2),eta1R,Dbdy(ii),-dr(ii),-dr(ii+1));
                Intdxpsi5QnTetaIn=-IZQnTetaIn+Z0(1).*IQnTetaIn;
                BvectP(Nm*Ns*n+Nm/2+ii)=-Intdxpsi5QnTetaIn*0+IndxGamma51InQn(IdIn)-IndxGamma51InQn(IdIn-1);
            end
            IdIn=IdIn+2;
        end
        
    end 
   
    
%      assignin('base','AmatI',AmatI);
%      assignin('base','BvectHI',BvectHI);
%      pause;
    
    
end