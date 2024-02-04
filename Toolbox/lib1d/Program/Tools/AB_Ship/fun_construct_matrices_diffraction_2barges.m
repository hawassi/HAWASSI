%% function for constructing matrix system for solving diffraction
%  17 December 2019
%  applying matching condition at ship boundary;
%  domain composition: 
%  layout: % ....(-(b+B))------(-b)....0....(b)------(b+B)....
%  domain: %   1           2           3          4         5        

%Amat.*Xvect=Bvect
%Xvect=[R1 R3 Tr3 Tr5 tau0 tau1 tau2 tau3 eps1m eps21m eps22m eps31m eps32m eps41m eps42m eps5m]
function [Amat,Bvect]=fun_construct_matrices_diffraction_2barges(tauAzwl,tauGwl,dxtauAzwl,dxtauGwl,Dwl,dxDwl,dr,Evmodes,kappaS,X0,phiWl,dxphiWl,etaWl,dxetaWl,XWl)
Ns=2;
N=4+4*Ns*Evmodes;
eta1L=etaWl(1);eta1R=etaWl(2);
eta2L=etaWl(3);eta2R=etaWl(4);
dxeta1L=dxetaWl(1);dxeta1R=dxetaWl(2);
dxeta2L=dxetaWl(3);dxeta2R=dxetaWl(4);
X1L=XWl(1);X1R=XWl(2);
X2L=XWl(3);X2R=XWl(4);

Amat=zeros(N,N);
Bvect=zeros(N,1);

IVp0Vp01L=fun_IntegVarphimVarphin(0,0,Dwl(1),-dr(1));
IVp0Vp01R=fun_IntegVarphimVarphin(0,0,Dwl(2),-dr(1));
IVp0Vp02L=fun_IntegVarphimVarphin(0,0,Dwl(3),-dr(2));
IVp0Vp02R=fun_IntegVarphimVarphin(0,0,Dwl(4),-dr(2));

%Xvect=              [tau0 tau1 tau2 tau3 ]
Amat(1,1:4)        =[-IVp0Vp01L , -IVp0Vp01L*tauGwl(1), 0         , 0];
Amat(2,1:4)        =[-IVp0Vp01R , -IVp0Vp01R*tauGwl(2), 0         , 0];
Amat(3,1:4)        =[0          , 0                     ,-IVp0Vp02L ,-IVp0Vp02L*tauGwl(3)];
Amat(4,1:4)        =[0          , 0                     ,-IVp0Vp02R ,-IVp0Vp02R*tauGwl(4)];


IQ0Vp01L=fun_IntegQmVarphin(0,0,kappaS(1,1),0,etaWl(1),Dwl(1),-dr(1));
IQ0Vp01R=fun_IntegQmVarphin(0,0,kappaS(1,2),0,etaWl(2),Dwl(2),-dr(1));
IQ0Vp02L=fun_IntegQmVarphin(0,0,kappaS(1,3),0,etaWl(3),Dwl(3),-dr(2));
IQ0Vp02R=fun_IntegQmVarphin(0,0,kappaS(1,4),0,etaWl(4),Dwl(4),-dr(2));

Bvect(1)=-phiWl(1).*IQ0Vp01L+tauAzwl(1)*IVp0Vp01L;
Bvect(2)=-phiWl(2).*IQ0Vp01R+tauAzwl(2)*IVp0Vp01R;
Bvect(3)=-phiWl(3).*IQ0Vp02L+tauAzwl(3)*IVp0Vp02L;
Bvect(4)=-phiWl(4).*IQ0Vp02R+tauAzwl(4)*IVp0Vp02R;

if Evmodes>=1
    for n=1:Evmodes
        muS1L_n=n*pi/(Dwl(1)-dr(1)); %% mu for 1st ship
        muS1R_n=n*pi/(Dwl(2)-dr(1));
        muS2L_n=n*pi/(Dwl(3)-dr(2));
        muS2R_n=n*pi/(Dwl(4)-dr(2));
        
        IVp0Vpn1L=fun_IntegVarphimVarphin(0,n,Dwl(1),-dr(1));
        IVp0Vpn1R=fun_IntegVarphimVarphin(0,n,Dwl(2),-dr(1));
        IVp0Vpn2L=fun_IntegVarphimVarphin(0,n,Dwl(3),-dr(2));
        IVp0Vpn2R=fun_IntegVarphimVarphin(0,n,Dwl(4),-dr(2));
        
        IQnVp01L=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,Dwl(1),-dr(1));
        IQnVp01R=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,Dwl(2),-dr(1));
        IQnVp02L=fun_IntegQmVarphin(n,0,kappaS(n+1,3),0,eta2L,Dwl(3),-dr(2));
        IQnVp02R=fun_IntegQmVarphin(n,0,kappaS(n+1,4),0,eta2R,Dwl(4),-dr(2));
        
        
                       %Xvect=[ tau0         tau1                   tau2         tau3]
        Amat(4*Ns*(n-1)+5,1:4)=[-IVp0Vpn1L , -IVp0Vpn1L*tauGwl(1), 0         , 0];
        Amat(4*Ns*(n-1)+6,1:4)=[-IVp0Vpn1R , -IVp0Vpn1R*tauGwl(2), 0         , 0];
        Amat(4*Ns*(n-1)+7,1:4)=[0          , 0                     ,-IVp0Vpn2L ,-IVp0Vpn2L*tauGwl(3)];
        Amat(4*Ns*(n-1)+8,1:4)=[0          , 0                     ,-IVp0Vpn2R ,-IVp0Vpn2R*tauGwl(4)];
        %Xvect=[tau0   tau1   tau2 tau3]
        Amat(4*Ns*(n-1)+9,1:4) =[0, -IQnVp01L*dxtauGwl(1), 0, 0];
        Amat(4*Ns*(n-1)+10,1:4)=[0, -IQnVp01R*dxtauGwl(2), 0, 0];
        Amat(4*Ns*(n-1)+11,1:4)=[0, 0        , 0,-IQnVp02L*dxtauGwl(3)];
        Amat(4*Ns*(n-1)+12,1:4)=[0, 0        , 0,-IQnVp02R*dxtauGwl(4)];
        
            for m=1:Evmodes
                muS1L_m=m*pi/(Dwl(1)-dr(1));
                muS1R_m=m*pi/(Dwl(2)-dr(1));
                muS2L_m=m*pi/(Dwl(3)-dr(2));
                muS2R_m=m*pi/(Dwl(4)-dr(2));
                
                IVpmVpn1L=fun_IntegVarphimVarphin(m,n,Dwl(1),-dr(1));
                IVpmVpn1R=fun_IntegVarphimVarphin(m,n,Dwl(2),-dr(1));
                IVpmVpn2L=fun_IntegVarphimVarphin(m,n,Dwl(3),-dr(2));
                IVpmVpn2R=fun_IntegVarphimVarphin(m,n,Dwl(4),-dr(2));
                
                IQmVpn1L=fun_IntegQmVarphin(m,n,kappaS(m+1,1),muS1L_n,eta1L,Dwl(1),-dr(1));
                IQmVpn1R=fun_IntegQmVarphin(m,n,kappaS(m+1,2),muS1R_n,eta1R,Dwl(2),-dr(1));
                IQmVpn2L=fun_IntegQmVarphin(m,n,kappaS(m+1,3),muS2L_n,eta2L,Dwl(3),-dr(2));
                IQmVpn2R=fun_IntegQmVarphin(m,n,kappaS(m+1,4),muS2R_n,eta2R,Dwl(4),-dr(2));
                
                IQnVpm1L=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muS1L_m,eta1L,Dwl(1),-dr(1));
                IQnVpm1R=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muS1R_m,eta1R,Dwl(2),-dr(1));
                IQnVpm2L=fun_IntegQmVarphin(n,m,kappaS(n+1,3),muS2L_m,eta2L,Dwl(3),-dr(2));
                IQnVpm2R=fun_IntegQmVarphin(n,m,kappaS(n+1,4),muS2R_m,eta2R,Dwl(4),-dr(2));
                
                IQmQn1L=fun_IntegQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dwl(1),eta1L);
                IQmQn1R=fun_IntegQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dwl(2),eta1R);
                IQmQn2L=fun_IntegQmQn(m,n,kappaS(m+1,3),kappaS(n+1,3),eta2L,Dwl(3),eta2L);
                IQmQn2R=fun_IntegQmQn(m,n,kappaS(m+1,4),kappaS(n+1,4),eta2R,Dwl(4),eta2R);
                
                IdxQmQn1L=fun_IntegdxQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dwl(1),dxeta1L,dxDwl(1),eta1L);
                IdxQmQn1R=fun_IntegdxQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dwl(2),dxeta1R,dxDwl(2),eta1R);
                IdxQmQn2L=fun_IntegdxQmQn(m,n,kappaS(m+1,3),kappaS(n+1,3),eta2L,Dwl(3),dxeta2L,dxDwl(3),eta2L);
                IdxQmQn2R=fun_IntegdxQmQn(m,n,kappaS(m+1,4),kappaS(n+1,4),eta2R,Dwl(4),dxeta2R,dxDwl(4),eta2R);
                
                IdxVpmQn1L=fun_IntegdxVarphimQn(muS1L_m,kappaS(n+1,1),eta1L,Dwl(1),dxDwl(1),-dr(1));
                IdxVpmQn1R=fun_IntegdxVarphimQn(muS1R_m,kappaS(n+1,2),eta1R,Dwl(2),dxDwl(2),-dr(1));
                IdxVpmQn2L=fun_IntegdxVarphimQn(muS2L_m,kappaS(n+1,3),eta2L,Dwl(3),dxDwl(3),-dr(2));
                IdxVpmQn2R=fun_IntegdxVarphimQn(muS2R_m,kappaS(n+1,4),eta2R,Dwl(4),dxDwl(4),-dr(2));
                
                 
                if n==1
                    IQmVp01L=fun_IntegQmVarphin(m,0,kappaS(m+1,1),0,eta1L,Dwl(1),-dr(1));
                    IQmVp01R=fun_IntegQmVarphin(m,0,kappaS(m+1,2),0,eta1R,Dwl(2),-dr(1));
                    IQmVp02L=fun_IntegQmVarphin(m,0,kappaS(m+1,3),0,eta2L,Dwl(3),-dr(2));
                    IQmVp02R=fun_IntegQmVarphin(m,0,kappaS(m+1,4),0,eta2R,Dwl(4),-dr(2));
                    IVpmVp01L=fun_IntegVarphimVarphin(m,0,Dwl(1),-dr(1));
                    IVpmVp01R=fun_IntegVarphimVarphin(m,0,Dwl(2),-dr(1));
                    IVpmVp02L=fun_IntegVarphimVarphin(m,0,Dwl(3),-dr(2));
                    IVpmVp02R=fun_IntegVarphimVarphin(m,0,Dwl(4),-dr(2));
                                                          %Xvect=[  eps1m    eps21m                              eps22m                                eps31m                                 eps32m                               eps41m                              eps42m                               eps5m]
                    Amat(1,4*Ns*(m-1)+(5:12))=[IQmVp01L , -IVpmVp01L                         , -exp(-muS1L_m*(X1L-X1R))*IVpmVp01L  , 0                                     , 0                                    , 0                                , 0                                 , 0];
                    Amat(2,4*Ns*(m-1)+(5:12))=[0        , -exp(muS1R_m*(X1R-X1L))*IVpmVp01R  , -IVpmVp01R                          , IQmVp01R                              , exp(kappaS(m+1,2)*(X1R-X2L))*IQmVp01R, 0                                , 0                                 , 0];
                    Amat(3,4*Ns*(m-1)+(5:12))=[0        , 0                                  , 0                                   , exp(-kappaS(m+1,3)*(X2L-X1R))*IQmVp02L, IQmVp02L                             , -IVpmVp02L                       , -exp(-muS2L_m*(X2L-X2R))*IVpmVp02L, 0];
                    Amat(4,4*Ns*(m-1)+(5:12))=[0        , 0                                  , 0                                   , 0                                     , 0                                    , -exp(muS2R_m*(X2R-X2L))*IVpmVp02R, -IVpmVp02R                       , IQmVp02R];

                
                end
                                                %Xvect=[  eps1m    eps21m                              eps22m                                eps31m                                 eps32m                               eps41m                              eps42m                               eps5m]
                Amat(4*Ns*(n-1)+5,4*Ns*(m-1)+(5:12))=[IQmVpn1L , -IVpmVpn1L                         , -exp(-muS1L_m*(X1L-X1R))*IVpmVpn1L  , 0                                     , 0                                    , 0                                , 0                                 , 0];
                Amat(4*Ns*(n-1)+6,4*Ns*(m-1)+(5:12))=[0        , -exp(muS1R_m*(X1R-X1L))*IVpmVpn1R  , -IVpmVpn1R                          , IQmVpn1R                              , exp(kappaS(m+1,2)*(X1R-X2L))*IQmVpn1R, 0                                , 0                                 , 0];
                Amat(4*Ns*(n-1)+7,4*Ns*(m-1)+(5:12))=[0        , 0                                  , 0                                   , exp(-kappaS(m+1,3)*(X2L-X1R))*IQmVpn2L, IQmVpn2L                             , -IVpmVpn2L                       , -exp(-muS2L_m*(X2L-X2R))*IVpmVpn2L, 0];
                Amat(4*Ns*(n-1)+8,4*Ns*(m-1)+(5:12))=[0        , 0                                  , 0                                   , 0                                     , 0                                    , -exp(muS2R_m*(X2R-X2L))*IVpmVpn2R, -IVpmVpn2R                       , IQmVpn2R];
                %                                %Xvect=[  eps1m                          eps21m                                              eps22m                                                    eps31m                                                            eps32m                                                           eps41m                                               eps42m                                                      eps5m]
                Amat(4*Ns*(n-1)+9 ,4*Ns*(m-1)+(5:12))=[kappaS(m+1,1)*IQmQn1L+IdxQmQn1L, -(muS1L_m*IQnVpm1L+IdxVpmQn1L)                        , -exp(-muS1L_m*(X1L-X1R))*(-muS1L_m*IQnVpm1L+IdxVpmQn1L), 0                                                                , 0                                                             , 0                                                    , 0                                                        , 0    ];
                Amat(4*Ns*(n-1)+10,4*Ns*(m-1)+(5:12))=[0                              , -exp(muS1R_m*(X1R-X1L))*(muS1R_m*IQnVpm1R+IdxVpmQn1R) , -(-muS1R_m*IQnVpm1R+IdxVpmQn1R)                        , -kappaS(m+1,2)*IQmQn1R+IdxQmQn1R                                 , exp(kappaS(m+1,2)*(X1R-X2L))*(kappaS(m+1,2)*IQmQn1R+IdxQmQn1R), 0                                                    , 0                                                        , 0    ];
                Amat(4*Ns*(n-1)+11,4*Ns*(m-1)+(5:12))=[0                              , 0                                                     , 0                                                      , exp(-kappaS(m+1,3)*(X2L-X1R))*(-kappaS(m+1,3)*IQmQn2L+IdxQmQn2L) , kappaS(m+1,3)*IQmQn2L+IdxQmQn2L                               , -(muS2L_m*IQnVpm2L+IdxVpmQn2L)                       , -exp(-muS2L_m*(X2L-X2R))*(-muS2L_m*IQnVpm2L+IdxVpmQn2L)  , 0    ];
                Amat(4*Ns*(n-1)+12,4*Ns*(m-1)+(5:12))=[0                              , 0                                                     , 0                                                      , 0                                                                , 0                                                             , -exp(muS2R_m*(X2R-X2L))*(muS2R_m*IQnVpm2R+IdxVpmQn2R), -(-muS2R_m*IQnVpm2R+IdxVpmQn2R)                          , -kappaS(m+1,4)*IQmQn2R+IdxQmQn2R];
            end
     
            
            IQ0Vpn1L=fun_IntegQmVarphin(0,n,kappaS(1,1),muS1L_n,eta1L,Dwl(1),-dr(1));
            IQ0Vpn1R=fun_IntegQmVarphin(0,n,kappaS(1,2),muS1R_n,eta1R,Dwl(2),-dr(1));
            IQ0Vpn2L=fun_IntegQmVarphin(0,n,kappaS(1,3),muS2L_n,eta2L,Dwl(3),-dr(2));
            IQ0Vpn2R=fun_IntegQmVarphin(0,n,kappaS(1,4),muS2R_n,eta2R,Dwl(4),-dr(2));
            
            IQ0Qn1L=fun_IntegQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dwl(1),eta1L);
            IQ0Qn1R=fun_IntegQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dwl(2),eta1R);
            IQ0Qn2L=fun_IntegQmQn(0,n,kappaS(1,3),kappaS(n+1,3),eta2L,Dwl(3),eta2L);
            IQ0Qn2R=fun_IntegQmQn(0,n,kappaS(1,4),kappaS(n+1,4),eta2R,Dwl(4),eta2R);
            
            IdxQ0Qn1L=fun_IntegdxQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dwl(1),dxeta1L,dxDwl(1),eta1L);
            IdxQ0Qn1R=fun_IntegdxQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dwl(2),dxeta1R,dxDwl(2),eta1R);
            IdxQ0Qn2L=fun_IntegdxQmQn(0,n,kappaS(1,3),kappaS(n+1,3),eta2L,Dwl(1),dxeta2L,dxDwl(3),eta2L);
            IdxQ0Qn2R=fun_IntegdxQmQn(0,n,kappaS(1,4),kappaS(n+1,4),eta2R,Dwl(2),dxeta2R,dxDwl(4),eta2R);


            
        Bvect(4*Ns*(n-1)+5)=-phiWl(1).*IQ0Vpn1L+tauAzwl(1)*IVp0Vpn1L;
        Bvect(4*Ns*(n-1)+6)=-phiWl(2).*IQ0Vpn1R+tauAzwl(2)*IVp0Vpn1R;
        Bvect(4*Ns*(n-1)+7)=-phiWl(3).*IQ0Vpn2L+tauAzwl(3)*IVp0Vpn2L;
        Bvect(4*Ns*(n-1)+8)=-phiWl(4).*IQ0Vpn2R+tauAzwl(4)*IVp0Vpn2R;
        
        Bvect(4*Ns*(n-1)+9) =-(dxphiWl(1).*IQ0Qn1L+phiWl(1)*IdxQ0Qn1L)+dxtauAzwl(1)*IQnVp01L;
        Bvect(4*Ns*(n-1)+10)=-(dxphiWl(2).*IQ0Qn1R+phiWl(2)*IdxQ0Qn1R)+dxtauAzwl(2)*IQnVp01R;
        Bvect(4*Ns*(n-1)+11)=-(dxphiWl(3).*IQ0Qn2L+phiWl(3)*IdxQ0Qn2L)+dxtauAzwl(3)*IQnVp02L;
        Bvect(4*Ns*(n-1)+12)=-(dxphiWl(4).*IQ0Qn2R+phiWl(4)*IdxQ0Qn2R)+dxtauAzwl(4)*IQnVp02R;
        
   
    end
end