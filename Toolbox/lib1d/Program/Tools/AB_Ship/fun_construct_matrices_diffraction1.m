function [Amat, Bvect]=fun_construct_matrices_diffraction1(tauAzwl,tauGwl,dxtauAzwl,dxtauGwl,Dwl,dxDwl,dr,Evmodes,kappaS,X0,phiWl,dxphiWl,etaWl,dxetaWl,XWl)
N=2+4*Evmodes;
eta1L=etaWl(1);eta1R=etaWl(2);
dxeta1L=dxetaWl(1);dxeta1R=dxetaWl(2);
X1L=XWl(1);X1R=XWl(2);

Amat=zeros(N,N);
Bvect=zeros(N,1);

IVp0Vp01L=fun_IntegVarphimVarphin(0,0,Dwl(1),-dr);
IVp0Vp01R=fun_IntegVarphimVarphin(0,0,Dwl(2),-dr);

                     % sig0           sig1
Amat(1,1:2)        =[-IVp0Vp01L , -IVp0Vp01L*tauGwl(1)];
Amat(2,1:2)        =[-IVp0Vp01R , -IVp0Vp01R*tauGwl(2)];

IQ0Vp01L=fun_IntegQmVarphin(0,0,kappaS(1,1),0,etaWl(1),Dwl(1),-dr);
IQ0Vp01R=fun_IntegQmVarphin(0,0,kappaS(1,2),0,etaWl(2),Dwl(2),-dr);

Bvect(1)=-phiWl(1).*IQ0Vp01L+tauAzwl(1)*IVp0Vp01L;
Bvect(2)=-phiWl(2).*IQ0Vp01R+tauAzwl(2)*IVp0Vp01R;


if Evmodes>=1
    for n=1:Evmodes
        muS1L_n=n*pi/(Dwl(1)-dr); 
        muS1R_n=n*pi/(Dwl(2)-dr);
        
        IVp0Vpn1L=fun_IntegVarphimVarphin(0,n,Dwl(1),-dr);
        IVp0Vpn1R=fun_IntegVarphimVarphin(0,n,Dwl(2),-dr);
        
        IQnVp01L=fun_IntegQmVarphin(n,0,kappaS(n+1,1),0,eta1L,Dwl(1),-dr);
        IQnVp01R=fun_IntegQmVarphin(n,0,kappaS(n+1,2),0,eta1R,Dwl(2),-dr);
                            % sig0           sig1
        Amat(4*(n-1)+3,1:2)=[-IVp0Vpn1L , -IVp0Vpn1L*tauGwl(1)];
        Amat(4*(n-1)+4,1:2)=[-IVp0Vpn1R , -IVp0Vpn1R*tauGwl(2)];
                            % sig0           sig1
        Amat(4*(n-1)+5,1:2) =[0 , -IQnVp01L*dxtauGwl(1)];
        Amat(4*(n-1)+6,1:2) =[0 , -IQnVp01R*dxtauGwl(2)];
        
        
        for m=1:Evmodes
            muS1L_m=m*pi/(Dwl(1)-dr);
            muS1R_m=m*pi/(Dwl(2)-dr);
            
            IVpmVpn1L=fun_IntegVarphimVarphin(m,n,Dwl(1),-dr(1));
            IVpmVpn1R=fun_IntegVarphimVarphin(m,n,Dwl(2),-dr(1));

            IQmVpn1L=fun_IntegQmVarphin(m,n,kappaS(m+1,1),muS1L_n,eta1L,Dwl(1),-dr(1));
            IQmVpn1R=fun_IntegQmVarphin(m,n,kappaS(m+1,2),muS1R_n,eta1R,Dwl(2),-dr(1));
   
            IQnVpm1L=fun_IntegQmVarphin(n,m,kappaS(n+1,1),muS1L_m,eta1L,Dwl(1),-dr(1));
            IQnVpm1R=fun_IntegQmVarphin(n,m,kappaS(n+1,2),muS1R_m,eta1R,Dwl(2),-dr(1));
            
            IQmQn1L=fun_IntegQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dwl(1),eta1L);
            IQmQn1R=fun_IntegQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dwl(2),eta1R);

            IdxQmQn1L=fun_IntegdxQmQn(m,n,kappaS(m+1,1),kappaS(n+1,1),eta1L,Dwl(1),dxeta1L,dxDwl(1),eta1L);
            IdxQmQn1R=fun_IntegdxQmQn(m,n,kappaS(m+1,2),kappaS(n+1,2),eta1R,Dwl(2),dxeta1R,dxDwl(2),eta1R);

            IdxVpmQn1L=fun_IntegdxVarphimQn(muS1L_m,kappaS(n+1,1),eta1L,Dwl(1),dxDwl(1),-dr(1));
            IdxVpmQn1R=fun_IntegdxVarphimQn(muS1R_m,kappaS(n+1,2),eta1R,Dwl(2),dxDwl(2),-dr(1));

                
            if n==1
                IQmVp01L=fun_IntegQmVarphin(m,0,kappaS(m+1,1),0,eta1L,Dwl(1),-dr(1));
                IQmVp01R=fun_IntegQmVarphin(m,0,kappaS(m+1,2),0,eta1R,Dwl(2),-dr(1));
                IVpmVp01L=fun_IntegVarphimVarphin(m,0,Dwl(1),-dr(1));
                IVpmVp01R=fun_IntegVarphimVarphin(m,0,Dwl(2),-dr(1));
                                 %Xvect=[  eps1m       eps21m                              eps22m                             eps3m]
                Amat(1,4*(m-1)+(3:6))=[IQmVp01L , -IVpmVp01L                         , -exp(-muS1L_m*(X1L-X1R))*IVpmVp01L  ,  0];
                Amat(2,4*(m-1)+(3:6))=[0        , -exp(muS1R_m*(X1R-X1L))*IVpmVp01R  , -IVpmVp01R                          ,  IQmVp01R];
      
            end
                                         %Xvect=[  eps1m       eps21m                              eps22m                             eps3m]
            Amat(4*(n-1)+3,4*(m-1)+(3:6))=[IQmVpn1L , -IVpmVpn1L                         , -exp(-muS1L_m*(X1L-X1R))*IVpmVpn1L  , 0    ];
            Amat(4*(n-1)+4,4*(m-1)+(3:6))=[0        , -exp(muS1R_m*(X1R-X1L))*IVpmVpn1R  , -IVpmVpn1R                          , IQmVpn1R  ];
            Amat(4*(n-1)+5,4*(m-1)+(3:6))=[kappaS(m+1,1)*IQmQn1L+IdxQmQn1L, -(muS1L_m*IQnVpm1L+IdxVpmQn1L)                        , -exp(-muS1L_m*(X1L-X1R))*(-muS1L_m*IQnVpm1L+IdxVpmQn1L), 0  ];
            Amat(4*(n-1)+6,4*(m-1)+(3:6))=[0                              , -exp(muS1R_m*(X1R-X1L))*(muS1R_m*IQnVpm1R+IdxVpmQn1R) , -(-muS1R_m*IQnVpm1R+IdxVpmQn1R)                        , -kappaS(m+1,2)*IQmQn1R+IdxQmQn1R];

           
        end
        
        IQ0Vpn1L=fun_IntegQmVarphin(0,n,kappaS(1,1),muS1L_n,eta1L,Dwl(1),-dr(1));
        IQ0Vpn1R=fun_IntegQmVarphin(0,n,kappaS(1,2),muS1R_n,eta1R,Dwl(2),-dr(1));
        IQ0Qn1L=fun_IntegQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dwl(1),eta1L);
        IQ0Qn1R=fun_IntegQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dwl(2),eta1R);
        IdxQ0Qn1L=fun_IntegdxQmQn(0,n,kappaS(1,1),kappaS(n+1,1),eta1L,Dwl(1),dxeta1L,dxDwl(1),eta1L);
        IdxQ0Qn1R=fun_IntegdxQmQn(0,n,kappaS(1,2),kappaS(n+1,2),eta1R,Dwl(2),dxeta1R,dxDwl(2),eta1R);

        Bvect(4*(n-1)+3)=-phiWl(1).*IQ0Vpn1L+tauAzwl(1)*IVp0Vpn1L;
        Bvect(4*(n-1)+4)=-phiWl(2).*IQ0Vpn1R+tauAzwl(2)*IVp0Vpn1R;
        Bvect(4*(n-1)+5)=-(dxphiWl(1).*IQ0Qn1L+phiWl(1)*IdxQ0Qn1L)+dxtauAzwl(1)*IQnVp01L;
        Bvect(4*(n-1)+6)=-(dxphiWl(2).*IQ0Qn1R+phiWl(2)*IdxQ0Qn1R)+dxtauAzwl(2)*IQnVp01R;
        
    end
end
