function [Amat, BvectS, BvectH, BvectP]=fun_construct_matrices_radiation(DL,DR,dxDL,dxDR,etaL,etaR,dxetaL,dxetaR,draft,b,modes,KappaL,KappaR,Idmotion,Z0)

N=4+4*modes;
Amat=zeros(N,N);
BvectH=zeros(N,1);
BvectP=zeros(N,1);
BvectS=zeros(N,1);

kappa0L=KappaL(1);
kappaL=KappaL(2:end);
kappa0R=KappaR(1);
kappaR=KappaR(2:end);


for n=0:modes
    muL_n=n*pi/(DL-draft);
    muR_n=n*pi/(DR-draft);
    I20nL=fun_IntegVarphimVarphin(0,n,DL,-draft);%fun_diff_Integ2mn(0,0,D,draft);
    I20nR=fun_IntegVarphimVarphin(0,n,DR,-draft);
    I30nL=fun_IntegQmVarphin(0,n,kappa0L,muL_n,etaL,DL,-draft);%fun_diff_Integ3mn(0,n,kappa0,mu_n,D,draft);
    I30nR=fun_IntegQmVarphin(0,n,kappa0R,muR_n,etaR,DR,-draft);


        
    if n>0
        I3n0L=fun_IntegQmVarphin(n,0,kappaL(n),0,etaL,DL,-draft);%fun_diff_Integ3mn(n,0,kappa(n),0,D,draft);
        I3n0R=fun_IntegQmVarphin(n,0,kappaR(n),0,etaR,DR,-draft);
        I10nL=fun_IntegQmQn(0,n,kappa0L,kappaL(n),etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
        I10nR=fun_IntegQmQn(0,n,kappa0R,kappaR(n),etaR,DR,etaR);
        I1dx0nL=fun_IntegdxQmQn(0,n,kappa0L,kappaL(n),etaL,DL,dxetaL,dxDL,etaL);
        I1dx0nR=fun_IntegdxQmQn(0,n,kappa0R,kappaR(n),etaR,DR,dxetaR,dxDR,etaR);

    else
        I3n0L=fun_IntegQmVarphin(0,0,kappa0L,0,etaL,DL,-draft);%fun_diff_Integ3mn(n,0,kappa(n),0,D,draft);
        I3n0R=fun_IntegQmVarphin(0,0,kappa0R,0,etaR,DR,-draft);
        I10nL=fun_IntegQmQn(0,0,kappa0L,kappa0L,etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
        I10nR=fun_IntegQmQn(0,0,kappa0R,kappa0R,etaR,DR,etaR);
        I1dx0nL=fun_IntegdxQmQn(0,0,kappa0L,kappa0L,etaL,DL,dxetaL,dxDL,etaL);
        I1dx0nR=fun_IntegdxQmQn(0,0,kappa0R,kappa0R,etaR,DR,dxetaR,dxDR,etaR);
    end
    
    
    Amat(4*n+1,1:4)=[I30nL                  , 0                        ,-I20nL ,I20nL*b];
    Amat(4*n+2,1:4)=[0                      , I30nR                    ,-I20nR ,-I20nR*b];
    Amat(4*n+3,1:4)=[kappa0L*I10nL+I1dx0nL  , 0                        ,0     ,-I3n0L];
    Amat(4*n+4,1:4)=[0                      , -kappa0R*I10nR+I1dx0nR   ,0     ,-I3n0R];
    
    if modes>=1
        for m=1:modes
            muL_m=m*pi/(DL-draft);
            muR_m=m*pi/(DR-draft);

            
            I2mnL=fun_IntegVarphimVarphin(m,n,DL,-draft);%fun_diff_Integ2mn(m,n,D,draft);
            I2mnR=fun_IntegVarphimVarphin(m,n,DR,-draft);
            
            I3mnL=fun_IntegQmVarphin(m,n,kappaL(m),muL_n,etaL,DL,-draft);%fun_diff_Integ3mn(m,n,kappa(m),mu_n,D,draft);
            I3mnR=fun_IntegQmVarphin(m,n,kappaR(m),muR_n,etaR,DR,-draft);
            
            if n>0
                I3nmL=fun_IntegQmVarphin(n,m,kappaL(n),muL_m,etaL,DL,-draft);%fun_diff_Integ3mn(n,m,kappa(n),mu_m,D,draft);
                I3nmR=fun_IntegQmVarphin(n,m,kappaR(n),muR_m,etaR,DR,-draft);
                
                I1mnL=fun_IntegQmQn(m,n,kappaL(m),kappaL(n),etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
                I1mnR=fun_IntegQmQn(m,n,kappaR(m),kappaR(n),etaR,DR,etaR);
                I1dxmnL=fun_IntegdxQmQn(m,n,kappaL(m),kappaL(n),etaL,DL,dxetaL,dxDL,etaL);
                I1dxmnR=fun_IntegdxQmQn(m,n,kappaR(m),kappaR(n),etaR,DR,dxetaR,dxDR,etaR);

            else
                I3nmL=fun_IntegQmVarphin(0,m,kappa0L,muL_m,etaL,DL,-draft);%fun_diff_Integ3mn(n,m,kappa(n),mu_m,D,draft);
                I3nmR=fun_IntegQmVarphin(0,m,kappa0R,muR_m,etaR,DR,-draft);
                
                I1mnL=fun_IntegQmQn(m,0,kappaL(m),kappa0L,etaL,DL,etaL);%fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
                I1mnR=fun_IntegQmQn(m,0,kappaR(m),kappa0R,etaR,DR,etaR);
                I1dxmnL=fun_IntegdxQmQn(m,0,kappaL(m),kappa0L,etaL,DL,dxetaL,dxDL,etaL);
                I1dxmnR=fun_IntegdxQmQn(m,0,kappaR(m),kappa0R,etaR,DR,dxetaR,dxDR,etaR);    
            end
            
            Amat(4*n+1,4*m+(1:4))=[I3mnL                    , -I2mnL                     , -exp(2*muL_m*b)*I2mnL        , 0                          ];
            Amat(4*n+2,4*m+(1:4))=[0                        , -exp(2*muR_m*b)*I2mnR       , -I2mnR                      ,  I3mnR                     ];
            Amat(4*n+3,4*m+(1:4))=[kappaL(m)*I1mnL+I1dxmnL   , -muL_m*I3nmL                , muL_m*exp(2*muL_m*b)*I3nmL  ,  0                         ];
            Amat(4*n+4,4*m+(1:4))=[0                        , -muR_m*exp(2*muR_m*b)*I3nmR  , muR_m*I3nmR                   , -kappaR(m)*I1mnR+I1dxmnR  ];
            
        end
    end
    
    kLn=KappaL(n+1);
    kRn=KappaR(n+1);
      
    DTL=DL-draft;
    DTR=DR-draft;
   
    IntQnTetaL=fun_IntegQn(kLn,etaL,DL,-draft,etaL);
    IntQnTetaR=fun_IntegQn(kRn,etaR,DR,-draft,etaR);
    IZD2varphinL = fun_IntegZD2Varphin(muL_n,DL,-DL,-draft);
    IZD2varphinR = fun_IntegZD2Varphin(muR_n,DR,-DR,-draft);
    IvarphinL = fun_IntegVarphin(muL_n,DL,-DL,-draft);
    IvarphinR = fun_IntegVarphin(muR_n,DR,-DR,-draft);
    
    IntQnDTL=fun_IntegQn(kLn,etaL,DL,-DL,-draft);
    IntQnDTR=fun_IntegQn(kRn,etaR,DR,-DR,-draft);
    
    IZDQnDTL = fun_IntegZDQn(kLn,etaL,DL,-DL,-draft);
    IZD2QnDTL = fun_IntegZD2Qn(kLn,etaL,DL,-DL,-draft);
    IZDQnDTR = fun_IntegZDQn(kRn,etaR,DR,-DR,-draft);
    IZD2QnDTR = fun_IntegZD2Qn(kRn,etaR,DR,-DR,-draft);
    
    
    if  strcmpi(Idmotion,'Surge') || strcmpi(Idmotion,'Free')% Surge
        BvectS(4*n+1)=0;
        BvectS(4*n+2)=0;
        BvectS(4*n+3)=IntQnTetaL;
        BvectS(4*n+4)=IntQnTetaR;
    end
        
    if strcmpi(Idmotion,'Heave') || strcmpi(Idmotion,'Free') % heave
         
        InGamma3LVarphin=(IZD2varphinL-b^2*IvarphinL)/2/DTL;
        InGamma3RVarphin=(IZD2varphinR-b^2*IvarphinR)/2/DTR;
        
        IndxGamma3LQn=(IntQnDTL*b/DTL)*(1+b*dxDL/2/DTL)+(dxDL/DTL)*(IZDQnDTL-IZD2QnDTL/2/DTL) ;
        IndxGamma3RQn=(IntQnDTR*b/DTR)*(-1+b*dxDR/2/DTR)+(dxDR/DTR)*(IZDQnDTR-IZD2QnDTR/2/DTR);
       
        BvectH(4*n+1)=InGamma3LVarphin;
        BvectH(4*n+2)=InGamma3RVarphin;
        BvectH(4*n+3)=IndxGamma3LQn;
        BvectH(4*n+4)=IndxGamma3RQn;
    end    
        
    if strcmpi(Idmotion,'Pitch')|| strcmpi(Idmotion,'Free')% pitch
        IZQnTetaL = fun_IntegZQn(kLn,etaL,DL,-draft,etaL);
        IZQnTetaR = fun_IntegZQn(kRn,etaR,DR,-draft,etaR);
        Intdxpsi5QnTetaL=-IZQnTetaL+Z0.*IntQnTetaL;
        Intdxpsi5QnTetaR=-IZQnTetaR+Z0.*IntQnTetaR;
       
        
            
        IntdxGamma5QnDTL=(IZD2QnDTL-b^2*IntQnDTL)./DTL/2+dxDL*(-b.*IZDQnDTL-(-b.*IZD2QnDTL+b^3/3*IntQnDTL)./2./DTL)./DTL;
        IntdxGamma5QnDTR=(IZD2QnDTR-b^2*IntQnDTR)./DTR/2+dxDR*(b.*IZDQnDTL-(b.*IZD2QnDTL-b^3/3*IntQnDTL)./2./DTL)./DTL;
        
        IntGamma5varphinDTL=(-b.*IZD2varphinL+b^3*IvarphinL/3)./DTL/2;
        IntGamma5varphinDTR=(b.*IZD2varphinR-b^3*IvarphinR/3)./DTR/2;
        
        BvectP(4*n+1)=IntGamma5varphinDTL;
        BvectP(4*n+2)=IntGamma5varphinDTR;
        BvectP(4*n+3)=Intdxpsi5QnTetaL+IntdxGamma5QnDTL;
        BvectP(4*n+4)=Intdxpsi5QnTetaR+IntdxGamma5QnDTR;
    end   
end




