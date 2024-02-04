function [Amat, Bvect]=fun_construct_matrices_diffraction(DL,DR,dxDL,dxDR,phiL,phiR,dxphiL,dxphiR,etaL,etaR,dxetaL,dxetaR,draft,b,modes,KappaL,KappaR)

N=2+4*modes;
Amat=zeros(N,N);
Bvect=zeros(N,1);

kappa0L=KappaL(1);
kappaL=KappaL(2:end);
kappa0R=KappaR(1);
kappaR=KappaR(2:end);

I200L=fun_IntegVarphimVarphin(0,0,DL,-draft);%fun_diff_Integ2mn(0,0,D,draft);
I200R=fun_IntegVarphimVarphin(0,0,DR,-draft);

Amat(1,1:2)        =[I200L         , -I200L ];
Amat(2,1:2)        =[I200R         ,  I200R ];

I300L=fun_IntegQmVarphin(0,0,kappa0L,0,etaL,DL,-draft);%fun_diff_Integ3mn(0,0,kappa0,0,D,draft);
I300R=fun_IntegQmVarphin(0,0,kappa0R,0,etaR,DR,-draft);
Bvect(1)=phiL.*I300L;
Bvect(2)=phiR.*I300R;

if modes>=1
    for n=1:modes
        muL_n=n*pi/(DL-draft);
        muR_n=n*pi/(DR-draft);
        
        
        %     Amat(4*(n-1)+3,1:2)=[0            , 0       ];
        %     Amat(4*(n-1)+4,1:2)=[0            , 0       ];
        
        I3n0L=fun_IntegQmVarphin(n,0,kappaL(n),0,etaL,DL,-draft);%fun_diff_Integ3mn(n,0,kappa(n),0,D,draft);
        I3n0R=fun_IntegQmVarphin(n,0,kappaR(n),0,etaR,DR,-draft);
        
        Amat(4*(n-1)+5,1:2)=[0            , I3n0L/b      ];
        Amat(4*(n-1)+6,1:2)=[0            , I3n0R/b      ];
        
        for m=1:modes
            muL_m=m*pi/(DL-draft);
            muR_m=m*pi/(DR-draft);
            
            
            if n==1
                I3m0L=fun_IntegQmVarphin(m,0,kappaL(m),0,etaL,DL,-draft);%fun_diff_Integ3mn(m,0,kappa(m),0,D,draft);
                I3m0R=fun_IntegQmVarphin(m,0,kappaR(m),0,etaR,DR,-draft);
                
                Amat(1,4*(m-1)+(3:6))  =[-I3m0L        , 0  , 0    ,  0              ];
                Amat(2,4*(m-1)+(3:6))  =[ 0            , 0  , 0    , -I3m0R            ];
            end
            
            I1mnL=fun_IntegQmQn(m,n,kappaL(m),kappaL(n),etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
            I1mnR=fun_IntegQmQn(m,n,kappaR(m),kappaR(n),etaR,DR,etaR);
            I1dxmnL=fun_IntegdxQmQn(m,n,kappaL(m),kappaL(n),etaL,DL,dxetaL,dxDL,etaL);
            I1dxmnR=fun_IntegdxQmQn(m,n,kappaR(m),kappaR(n),etaR,DR,dxetaR,dxDR,etaR);
            
            I2mnL=fun_IntegVarphimVarphin(m,n,DL,-draft);%fun_diff_Integ2mn(m,n,D,draft);
            I2mnR=fun_IntegVarphimVarphin(m,n,DR,-draft);
            
            I3mnL=fun_IntegQmVarphin(m,n,kappaL(m),muL_n,etaL,DL,-draft);%fun_diff_Integ3mn(m,n,kappa(m),mu_n,D,draft);
            I3mnR=fun_IntegQmVarphin(m,n,kappaR(m),muR_n,etaR,DR,-draft);
            
            I3nmL=fun_IntegQmVarphin(n,m,kappaL(n),muL_m,etaL,DL,-draft);%fun_diff_Integ3mn(n,m,kappa(n),mu_m,D,draft);
            I3nmR=fun_IntegQmVarphin(n,m,kappaR(n),muR_m,etaL,DR,-draft);
            
            Amat(4*(n-1)+3,4*(m-1)+(3:6))=[-I3mnL                       , exp(-muL_m*b)*I2mnL      , exp(muL_m*b)*I2mnL      ,  0                          ];
            Amat(4*(n-1)+4,4*(m-1)+(3:6))=[ 0                           , exp(muR_m*b)*I2mnR       , exp(-muR_m*b)*I2mnR     , -I3mnR                      ];
            Amat(4*(n-1)+5,4*(m-1)+(3:6))=[-kappaL(m)*(I1mnL)-I1dxmnL    , muL_m*exp(-muL_m*b)*I3nmL , -muL_m*exp(muL_m*b)*I3nmL, 0                         ];
            Amat(4*(n-1)+6,4*(m-1)+(3:6))=[ 0                           , muR_m*exp(muR_m*b)*I3nmR  , -muR_m*exp(-muR_m*b)*I3nmR, kappaR(m)*(I1mnR)-I1dxmnR ];
        end
        
        
        I30nL=fun_IntegQmVarphin(0,n,kappa0L,muL_n,etaL,DL,-draft);%fun_diff_Integ3mn(0,n,kappa0,mu_n,D,draft);
        I30nR=fun_IntegQmVarphin(0,n,kappa0R,muR_n,etaR,DR,-draft);
        I10nL=fun_IntegQmQn(0,n,kappa0L,kappaL(n),etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
        I10nR=fun_IntegQmQn(0,n,kappa0R,kappaR(n),etaR,DR,etaR);
        I1dx0nL=fun_IntegdxQmQn(0,n,kappa0L,kappaL(n),etaL,DL,dxetaL,dxDL,etaL);
        I1dx0nR=fun_IntegdxQmQn(0,n,kappa0R,kappaR(n),etaR,DR,dxetaR,dxDR,etaR);
        Bvect(4*(n-1)+3)=phiL.*I30nL;
        Bvect(4*(n-1)+4)=phiR.*I30nR;
        Bvect(4*(n-1)+5)=dxphiL*I10nL+phiL*I1dx0nL;
        Bvect(4*(n-1)+6)=dxphiR*I10nR+phiR*I1dx0nR;
        
    end
end

end