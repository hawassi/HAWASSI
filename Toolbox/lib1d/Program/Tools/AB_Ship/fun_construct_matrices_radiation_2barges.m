%% function for constructing matrix system for solving radiation
%  17 December 2019
%  applying matching condition at ship boundary;
%  domain composition: 
%  layout: % ....(-(b+B))------(-b)....0....(b)------(b+B)....
%  domain: %   1           2           3          4         5        

%Amat.*Xvect=Bvect
%Xvect=[A1 A3L A3R A5 tau0 tau1 tau2 tau3 eps1m eps21m eps22m eps31m eps32m eps41m eps42m eps5m]
function [Amat,BvectS,BvectH,BvectP] =fun_construct_matrices_radiation_2barges(D,dxD,dr,Evmodes,kappaS,Idmotion,X0,Z0,etaWl,dxetaWl,XWl)
Ns=2;
N=4*Ns+4*Ns*Evmodes;
eta1L=etaWl(1);eta1R=etaWl(2);
eta2L=etaWl(3);eta2R=etaWl(4);
dxeta1L=dxetaWl(1);dxeta1R=dxetaWl(2);
dxeta2L=dxetaWl(3);dxeta2R=dxetaWl(4);
X1L=XWl(1);X1R=XWl(2);
X2L=XWl(3);X2R=XWl(4);

Amat=zeros(N,N);
BvectH=zeros(N,1);
BvectP=zeros(N,1);
BvectS=zeros(N,1);

for n=0:Evmodes
    muS1L_n=n*pi/(D(1)-dr(1)); %% mu for 1st ship 
    muS1R_n=n*pi/(D(2)-dr(1));
    muS2L_n=n*pi/(D(3)-dr(2));
    muS2R_n=n*pi/(D(4)-dr(2));
    
    IVp0Vpn1L=fun_IntegVarphimVarphin(0,n,D(1),-dr(1));
    IVp0Vpn1R=fun_IntegVarphimVarphin(0,n,D(2),-dr(1));
    IVp0Vpn2L=fun_IntegVarphimVarphin(0,n,D(3),-dr(2));
    IVp0Vpn2R=fun_IntegVarphimVarphin(0,n,D(4),-dr(2));
    
    IQ0Vpn1L=fun_IntegQmVarphin(0,n,kappaS(1,1),muS1L_n,eta1L,D(1),-dr(1));
    IQ0Vpn1R=fun_IntegQmVarphin(0,n,kappaS(1,2),muS1R_n,eta1R,D(2),-dr(1));
    IQ0Vpn2L=fun_IntegQmVarphin(0,n,kappaS(1,3),muS2L_n,eta2L,D(3),-dr(2));
    IQ0Vpn2R=fun_IntegQmVarphin(0,n,kappaS(1,4),muS2R_n,eta2R,D(4),-dr(2));
    
    
    IQnVp01L=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,D(1),-dr(1));
    IQnVp01R=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,D(2),-dr(1));
    IQnVp02L=fun_IntegQmVarphin(n,0,kappaS(n+1,3),0,eta2L,D(3),-dr(2));
    IQnVp02R=fun_IntegQmVarphin(n,0,kappaS(n+1,4),0,eta2R,D(4),-dr(2));
    
    IQ0Qn1L=fun_IntegQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,D(1),eta1L);
    IQ0Qn1R=fun_IntegQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,D(2),eta1R);
    IQ0Qn2L=fun_IntegQmQn(0,n,kappaS(1,3),kappaS(n+1,3),eta2L,D(3),eta2L);
    IQ0Qn2R=fun_IntegQmQn(0,n,kappaS(1,4),kappaS(n+1,4),eta2R,D(4),eta2R);
    
    IdxQ0Qn1L=fun_IntegdxQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,D(1),dxeta1L,dxD(1),eta1L);
    IdxQ0Qn1R=fun_IntegdxQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,D(2),dxeta1R,dxD(2),eta1R);
    IdxQ0Qn2L=fun_IntegdxQmQn(0,n,kappaS(1,3),kappaS(n+1,3),eta2L,D(1),dxeta2L,dxD(3),eta2L);
    IdxQ0Qn2R=fun_IntegdxQmQn(0,n,kappaS(1,4),kappaS(n+1,4),eta2R,D(2),dxeta2R,dxD(4),eta2R);
                     
                   %Xvect=[A1         A3L                                 A3R                                  A5         tau0         tau1                  tau2        tau3 ]
    Amat(4*Ns*n+1,1:4*Ns)=[IQ0Vpn1L , 0                                  , 0                                   , 0       , -IVp0Vpn1L , -IVp0Vpn1L*(X1L-X0(1)), 0         , 0];
    Amat(4*Ns*n+2,1:4*Ns)=[0        , IQ0Vpn1R                           , exp(-kappaS(1,2)*(X1R-X2L))*IQ0Vpn1R, 0       , -IVp0Vpn1R , -IVp0Vpn1R*(X1R-X0(1)), 0         , 0];
    Amat(4*Ns*n+3,1:4*Ns)=[0        , exp(kappaS(1,3)*(X2L-X1R))*IQ0Vpn2L, IQ0Vpn2L                            , 0       , 0          , 0                     ,-IVp0Vpn2L ,-IVp0Vpn2L*(X2L-X0(2))];
    Amat(4*Ns*n+4,1:4*Ns)=[0        , 0                                  ,  0                                  , IQ0Vpn2R, 0          , 0                     ,-IVp0Vpn2R ,-IVp0Vpn2R*(X2R-X0(2))];
   
    Amat(4*Ns*n+5,1:4*Ns)=[-kappaS(1,1)*IQ0Qn1L+IdxQ0Qn1L , 0                                                          , 0                                                           , 0                              , 0, -IQnVp01L, 0, 0];
    Amat(4*Ns*n+6,1:4*Ns)=[0                              , kappaS(1,2)*IQ0Qn1R+IdxQ0Qn1R                              , exp(-kappaS(1,2)*(X1R-X2L))*(-kappaS(1,2)*IQ0Qn1R+IdxQ0Qn1R), 0                              , 0, -IQnVp01R, 0, 0];
    Amat(4*Ns*n+7,1:4*Ns)=[0                              , exp(kappaS(1,3)*(X2L-X1R))*(kappaS(1,3)*IQ0Qn2L+IdxQ0Qn2L) , -kappaS(1,3)*IQ0Qn2L+IdxQ0Qn2L                              , 0                              , 0, 0        , 0,-IQnVp02L];
    Amat(4*Ns*n+8,1:4*Ns)=[0                              , 0                                                          , 0                                                           , kappaS(1,4)*IQ0Qn2R+IdxQ0Qn2R  , 0, 0        , 0,-IQnVp02R];
    
   if Evmodes>0
       for m=1:Evmodes
           muS1L_m=m*pi/(D(1)-dr(1)); 
           muS1R_m=m*pi/(D(2)-dr(1));
           muS2L_m=m*pi/(D(3)-dr(2));
           muS2R_m=m*pi/(D(4)-dr(2));
           
           IVpmVpn1L=fun_IntegVarphimVarphin(m,n,D(1),-dr(1));
           IVpmVpn1R=fun_IntegVarphimVarphin(m,n,D(2),-dr(1));
           IVpmVpn2L=fun_IntegVarphimVarphin(m,n,D(3),-dr(2));
           IVpmVpn2R=fun_IntegVarphimVarphin(m,n,D(4),-dr(2));
           
           IQmVpn1L=fun_IntegQmVarphin(m,n,kappaS(m+1,1),muS1L_n,eta1L,D(1),-dr(1));
           IQmVpn1R=fun_IntegQmVarphin(m,n,kappaS(m+1,2),muS1R_n,eta1R,D(2),-dr(1));
           IQmVpn2L=fun_IntegQmVarphin(m,n,kappaS(m+1,3),muS2L_n,eta2L,D(3),-dr(2));
           IQmVpn2R=fun_IntegQmVarphin(m,n,kappaS(m+1,4),muS2R_n,eta2R,D(4),-dr(2));
           
           IQnVpm1L=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muS1L_m,eta1L,D(1),-dr(1));
           IQnVpm1R=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muS1R_m,eta1R,D(2),-dr(1));
           IQnVpm2L=fun_IntegQmVarphin(n,m,kappaS(n+1,3),muS2L_m,eta2L,D(3),-dr(2));
           IQnVpm2R=fun_IntegQmVarphin(n,m,kappaS(n+1,4),muS2R_m,eta2R,D(4),-dr(2));
           
           IQmQn1L=fun_IntegQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,D(1),eta1L);
           IQmQn1R=fun_IntegQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,D(2),eta1R);
           IQmQn2L=fun_IntegQmQn(m,n,kappaS(m+1,3),kappaS(n+1,3),eta2L,D(3),eta2L);
           IQmQn2R=fun_IntegQmQn(m,n,kappaS(m+1,4),kappaS(n+1,4),eta2R,D(4),eta2R);
           
           IdxQmQn1L=fun_IntegdxQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,D(1),dxeta1L,dxD(1),eta1L);
           IdxQmQn1R=fun_IntegdxQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,D(2),dxeta1R,dxD(2),eta1R);
           IdxQmQn2L=fun_IntegdxQmQn(m,n,kappaS(m+1,3),kappaS(n+1,3),eta2L,D(3),dxeta2L,dxD(3),eta2L);
           IdxQmQn2R=fun_IntegdxQmQn(m,n,kappaS(m+1,4),kappaS(n+1,4),eta2R,D(4),dxeta2R,dxD(4),eta2R);
           
           IdxVpmQn1L=fun_IntegdxVarphimQn(muS1L_m,kappaS(n+1,1),eta1L,D(1),dxD(1),-dr(1));
           IdxVpmQn1R=fun_IntegdxVarphimQn(muS1R_m,kappaS(n+1,2),eta1R,D(2),dxD(2),-dr(1));
           IdxVpmQn2L=fun_IntegdxVarphimQn(muS2L_m,kappaS(n+1,3),eta2L,D(3),dxD(3),-dr(2));
           IdxVpmQn2R=fun_IntegdxVarphimQn(muS2R_m,kappaS(n+1,4),eta2R,D(4),dxD(4),-dr(2));
          
                                   %Xvect=[  eps1m    eps21m                              eps22m                                eps31m                                 eps32m                               eps41m                              eps42m                               eps5m]
           Amat(4*Ns*n+1,4*Ns*m+(1:4*Ns))=[IQmVpn1L , -IVpmVpn1L                         , -exp(-muS1L_m*(X1L-X1R))*IVpmVpn1L  , 0                                    , 0                                    , 0                                , 0                                 , 0];
           Amat(4*Ns*n+2,4*Ns*m+(1:4*Ns))=[0        , -exp(muS1R_m*(X1R-X1L))*IVpmVpn1R  , -IVpmVpn1R                          , IQmVpn1R                             , exp(kappaS(m+1,2)*(X1R-X2L))*IQmVpn1R, 0                                , 0                                 , 0];
           Amat(4*Ns*n+3,4*Ns*m+(1:4*Ns))=[0        , 0                                  , 0                                   , exp(-kappaS(m+1,3)*(X2L-X1R))*IQmVpn2L, IQmVpn2L                             , -IVpmVpn2L                       , -exp(-muS2L_m*(X2L-X2R))*IVpmVpn2L, 0];
           Amat(4*Ns*n+4,4*Ns*m+(1:4*Ns))=[0        , 0                                  , 0                                   , 0                                    , 0                                    , -exp(muS2R_m*(X2R-X2L))*IVpmVpn2R, -IVpmVpn2R                       , IQmVpn2R];
%                                    %Xvect=[  eps1m                          eps21m                                              eps22m                                                    eps31m                                                         eps32m                                                           eps41m                                               eps42m                                                      eps5m]
           Amat(4*Ns*n+5,4*Ns*m+(1:4*Ns))=[kappaS(m+1,1)*IQmQn1L+IdxQmQn1L, -(muS1L_m*IQnVpm1L+IdxVpmQn1L)                        , -exp(-muS1L_m*(X1L-X1R))*(-muS1L_m*IQnVpm1L+IdxVpmQn1L), 0                                                            , 0                                                             , 0                                                    , 0                                                        , 0    ];
           Amat(4*Ns*n+6,4*Ns*m+(1:4*Ns))=[0                              , -exp(muS1R_m*(X1R-X1L))*(muS1R_m*IQnVpm1R+IdxVpmQn1R) , -(-muS1R_m*IQnVpm1R+IdxVpmQn1R)                        , -kappaS(m+1,2)*IQmQn1R+IdxQmQn1R                             , exp(kappaS(m+1,2)*(X1R-X2L))*(kappaS(m+1,2)*IQmQn1R+IdxQmQn1R), 0                                                    , 0                                                        , 0    ];
           Amat(4*Ns*n+7,4*Ns*m+(1:4*Ns))=[0                              , 0                                                     , 0                                                      , exp(-kappaS(m+1,3)*(X2L-X1R))*(-kappaS(m+1,3)*IQmQn2L+IdxQmQn2L) , kappaS(m+1,3)*IQmQn2L+IdxQmQn2L                                 , -(muS2L_m*IQnVpm2L+IdxVpmQn2L)                 , -exp(-muS2L_m*(X2L-X2R))*(-muS2L_m*IQnVpm2L+IdxVpmQn2L)  , 0    ];
           Amat(4*Ns*n+8,4*Ns*m+(1:4*Ns))=[0                              , 0                                                     , 0                                                      , 0                                                            , 0                                                             , -exp(muS2R_m*(X2R-X2L))*(muS2R_m*IQnVpm2R+IdxVpmQn2R), -(-muS2R_m*IQnVpm2R+IdxVpmQn2R)                          , -kappaS(m+1,4)*IQmQn2R+IdxQmQn2R];
     
         
       
       end
   end
   
   DT1L=D(1)-dr(1);DT1R=D(2)-dr(1);
   DT2L=D(3)-dr(2);DT2R=D(4)-dr(2);
   
   IntQnTeta1L=fun_IntegQn(kappaS(n+1,1),eta1L,D(1),-dr(1),eta1L);
   IntQnTeta1R=fun_IntegQn(kappaS(n+1,2),eta1R,D(2),-dr(1),eta1R);
   IntQnTeta2L=fun_IntegQn(kappaS(n+1,3),eta2L,D(3),-dr(2),eta2L);
   IntQnTeta2R=fun_IntegQn(kappaS(n+1,4),eta2R,D(4),-dr(2),eta2R);
   
   IZD2varphin1L = fun_IntegZD2Varphin(muS1L_n,D(1),-D(1),-dr(1));
   IZD2varphin1R = fun_IntegZD2Varphin(muS1R_n,D(2),-D(2),-dr(1));
   IZD2varphin2L = fun_IntegZD2Varphin(muS2L_n,D(3),-D(3),-dr(2));
   IZD2varphin2R = fun_IntegZD2Varphin(muS2R_n,D(4),-D(4),-dr(2));
   
   Ivarphin1L = fun_IntegVarphin(muS1L_n,D(1),-D(1),-dr(1));
   Ivarphin1R = fun_IntegVarphin(muS1R_n,D(2),-D(2),-dr(1));
   Ivarphin2L = fun_IntegVarphin(muS2L_n,D(3),-D(3),-dr(2));
   Ivarphin2R = fun_IntegVarphin(muS2R_n,D(4),-D(4),-dr(2));
  
   IntQnDT1L=fun_IntegQn(kappaS(n+1,1),eta1L,D(1),-D(1),-dr(1));
   IntQnDT1R=fun_IntegQn(kappaS(n+1,2),eta1R,D(2),-D(2),-dr(1));
   IntQnDT2L=fun_IntegQn(kappaS(n+1,3),eta2L,D(3),-D(3),-dr(2));
   IntQnDT2R=fun_IntegQn(kappaS(n+1,4),eta2R,D(4),-D(4),-dr(2));
   
   IZDQnDT1L = fun_IntegZDQn(kappaS(n+1,1),eta1L,D(1),-D(1),-dr(1));
   IZDQnDT1R = fun_IntegZDQn(kappaS(n+1,2),eta1R,D(2),-D(2),-dr(1));
   IZDQnDT2L = fun_IntegZDQn(kappaS(n+1,3),eta2L,D(3),-D(3),-dr(2));
   IZDQnDT2R = fun_IntegZDQn(kappaS(n+1,4),eta2R,D(4),-D(4),-dr(2));
   
   IZD2QnDT1L = fun_IntegZD2Qn(kappaS(n+1,1),eta1L,D(1),-D(1),-dr(1));
   IZD2QnDT1R = fun_IntegZD2Qn(kappaS(n+1,2),eta1R,D(2),-D(2),-dr(1));
   IZD2QnDT2L = fun_IntegZD2Qn(kappaS(n+1,3),eta2L,D(3),-D(3),-dr(2));
   IZD2QnDT2R = fun_IntegZD2Qn(kappaS(n+1,4),eta2R,D(4),-D(4),-dr(2));
   
   if  strcmpi(Idmotion,'Surge') || strcmpi(Idmotion,'Free')% Surge
        BvectS(4*Ns*n+1)=0;
        BvectS(4*Ns*n+2)=0;
        BvectS(4*Ns*n+3)=0;
        BvectS(4*Ns*n+4)=0;
        BvectS(4*Ns*n+5)=IntQnTeta1L;
        BvectS(4*Ns*n+6)=IntQnTeta1R;
        BvectS(4*Ns*n+7)=IntQnTeta2L;
        BvectS(4*Ns*n+8)=IntQnTeta2R;    
    end
    
    if strcmpi(Idmotion,'Heave') || strcmpi(Idmotion,'Free') % heave  
        InGamma31LVarphin=(IZD2varphin1L-(X1L-X0(1))^2*Ivarphin1L)/2/DT1L;
        InGamma31RVarphin=(IZD2varphin1R-(X1R-X0(1))^2*Ivarphin1R)/2/DT1R;
        InGamma32LVarphin=(IZD2varphin2L-(X2L-X0(2))^2*Ivarphin2L)/2/DT2L;
        InGamma32RVarphin=(IZD2varphin2R-(X2R-X0(2))^2*Ivarphin2R)/2/DT2R;
        
        IndxGamma31LQn=(-IntQnDT1L*(X1L-X0(1))/DT1L)*(1-(X1L-X0(1))*dxD(1)/2/DT1L)+(dxD(1)/DT1L)*(IZDQnDT1L-IZD2QnDT1L/2/DT1L) ;
        IndxGamma31RQn=(-IntQnDT1R*(X1R-X0(1))/DT1R)*(1-(X1R-X0(1))*dxD(2)/2/DT1R)+(dxD(2)/DT1R)*(IZDQnDT1R-IZD2QnDT1R/2/DT1R) ;
        IndxGamma32LQn=(-IntQnDT2L*(X2L-X0(2))/DT2L)*(1-(X2L-X0(2))*dxD(3)/2/DT2L)+(dxD(3)/DT2L)*(IZDQnDT2L-IZD2QnDT2L/2/DT2L) ;
        IndxGamma32RQn=(-IntQnDT2R*(X2R-X0(2))/DT2R)*(1-(X2R-X0(2))*dxD(4)/2/DT2R)+(dxD(4)/DT2R)*(IZDQnDT2R-IZD2QnDT2R/2/DT2R) ;
     
        BvectH(4*Ns*n+1)=InGamma31LVarphin;
        BvectH(4*Ns*n+2)=InGamma31RVarphin;
        BvectH(4*Ns*n+3)=InGamma32LVarphin;
        BvectH(4*Ns*n+4)=InGamma32RVarphin;
        BvectH(4*Ns*n+5)=IndxGamma31LQn;
        BvectH(4*Ns*n+6)=IndxGamma31RQn;
        BvectH(4*Ns*n+7)=IndxGamma32LQn;
        BvectH(4*Ns*n+8)=IndxGamma32RQn; 
    end  
    
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')% pitch
        IZQnTeta1L = fun_IntegZQn(kappaS(n+1,1),eta1L,D(1),-dr(1),eta1L);
        IZQnTeta1R = fun_IntegZQn(kappaS(n+1,2),eta1R,D(2),-dr(1),eta1R);
        IZQnTeta2L = fun_IntegZQn(kappaS(n+1,3),eta2L,D(3),-dr(2),eta2L);
        IZQnTeta2R = fun_IntegZQn(kappaS(n+1,4),eta2R,D(4),-dr(2),eta2R);
        
        Intdxpsi5QnTeta1L=-IZQnTeta1L+Z0(1).*IntQnTeta1L;
        Intdxpsi5QnTeta1R=-IZQnTeta1R+Z0(1).*IntQnTeta1R;
        Intdxpsi5QnTeta2L=-IZQnTeta2L+Z0(2).*IntQnTeta2L;
        Intdxpsi5QnTeta2R=-IZQnTeta2R+Z0(2).*IntQnTeta2R;
            
        IntdxGamma5QnDT1L=(IZD2QnDT1L-(X1L-X0(1))^2*IntQnDT1L)./DT1L/2+dxD(1)*((X1L-X0(1)).*IZDQnDT1L-((X1L-X0(1)).*IZD2QnDT1L-(X1L-X0(1))^3/3*IntQnDT1L)./2./DT1L)./DT1L;
        IntdxGamma5QnDT1R=(IZD2QnDT1R-(X1R-X0(1))^2*IntQnDT1R)./DT1R/2+dxD(2)*((X1R-X0(1)).*IZDQnDT1R-((X1R-X0(1)).*IZD2QnDT1R-(X1R-X0(1))^3/3*IntQnDT1R)./2./DT1R)./DT1R;
        IntdxGamma5QnDT2L=(IZD2QnDT2L-(X2L-X0(2))^2*IntQnDT2L)./DT2L/2+dxD(3)*((X2L-X0(2)).*IZDQnDT2L-((X2L-X0(2)).*IZD2QnDT2L-(X2L-X0(2))^3/3*IntQnDT2L)./2./DT2L)./DT2L;
        IntdxGamma5QnDT2R=(IZD2QnDT2R-(X2R-X0(2))^2*IntQnDT2R)./DT2R/2+dxD(4)*((X2R-X0(2)).*IZDQnDT2R-((X2R-X0(2)).*IZD2QnDT2R-(X2R-X0(2))^3/3*IntQnDT2R)./2./DT2R)./DT2R;
        
         
        IntGamma5varphinDT1L=((X1L-X0(1)).*IZD2varphin1L-(X1L-X0(1))^3*Ivarphin1L/3)./DT1L/2;
        IntGamma5varphinDT1R=((X1R-X0(1)).*IZD2varphin1R-(X1R-X0(1))^3*Ivarphin1R/3)./DT1R/2;
        IntGamma5varphinDT2L=((X2L-X0(2)).*IZD2varphin2L-(X2L-X0(2))^3*Ivarphin2L/3)./DT2L/2;
        IntGamma5varphinDT2R=((X2R-X0(2)).*IZD2varphin2R-(X2R-X0(2))^3*Ivarphin2R/3)./DT2R/2;
        
        BvectP(4*Ns*n+1)=IntGamma5varphinDT1L;
        BvectP(4*Ns*n+2)=IntGamma5varphinDT1R;
        BvectP(4*Ns*n+3)=IntGamma5varphinDT2L;
        BvectP(4*Ns*n+4)=IntGamma5varphinDT2R;
        
        BvectP(4*Ns*n+5)=Intdxpsi5QnTeta1L+IntdxGamma5QnDT1L;
        BvectP(4*Ns*n+6)=Intdxpsi5QnTeta1R+IntdxGamma5QnDT1R;
        BvectP(4*Ns*n+7)=Intdxpsi5QnTeta2L+IntdxGamma5QnDT2L;
        BvectP(4*Ns*n+8)=Intdxpsi5QnTeta2R+IntdxGamma5QnDT2R;   
    end   
    
   
end