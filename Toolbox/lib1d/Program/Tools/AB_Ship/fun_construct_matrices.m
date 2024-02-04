function [Amat, Bvect]=fun_construct_matrices(D,draft,b,modes,Kappa)

N=2+4*modes;
Amat=zeros(N,N);
Bvect=zeros(N,1);

kappa0=Kappa(1);
kappa=Kappa(2:end);

for n=1:modes
    mu_n=n*pi/(D-draft);   
    if n==1
        I200=fun_diff_Integ2mn(0,0,D,draft);
        
        Amat(1,1:2)        =[I200         , -I200    ];
        Amat(2,1:2)        =[I200         ,  I200    ];
    end
%     Amat(4*(n-1)+3,1:2)=[0            , 0       ];
%     Amat(4*(n-1)+4,1:2)=[0            , 0       ];
    
    I3n0=fun_diff_Integ3mn(n,0,kappa(n),0,D,draft);
  
    Amat(4*(n-1)+5,1:2)=[0            , I3n0/b      ];
    Amat(4*(n-1)+6,1:2)=[0            , I3n0/b      ];
    
    for m=1:modes
        mu_m=m*pi/(D-draft);
        
        if n==1
            I3m0=fun_diff_Integ3mn(m,0,kappa(m),0,D,draft);
            
            Amat(1,4*(m-1)+(3:6))  =[-I3m0        , 0  , 0    ,  0              ];
            Amat(2,4*(m-1)+(3:6))  =[ 0           , 0  , 0    , -I3m0            ];
        end  
        I0mn=fun_diff_Integ0mn(m,n,kappa(m),kappa(n),D,draft);   
        I1mn=fun_diff_Integ1mna(m,n,kappa(m),kappa(n),D,draft);%Sulisz
        %I1mn=fun_diff_Integ1mn(n-1,kappa(m),D,draft);
        I2mn=fun_diff_Integ2mn(m,n,D,draft);
        I3mn=fun_diff_Integ3mn(m,n,kappa(m),mu_n,D,draft);
        I3nm=fun_diff_Integ3mn(n,m,kappa(n),mu_m,D,draft);
 
        Amat(4*(n-1)+3,4*(m-1)+(3:6))=[-I3mn                , exp(-mu_m*b)*I2mn      , exp(mu_m*b)*I2mn      ,  0                     ];
        Amat(4*(n-1)+4,4*(m-1)+(3:6))=[ 0                   , exp(mu_m*b)*I2mn       , exp(-mu_m*b)*I2mn     , -I3mn                  ];
        Amat(4*(n-1)+5,4*(m-1)+(3:6))=[-kappa(m)*(I0mn+I1mn), mu_m*exp(-mu_m*b)*I3nm , -mu_m*exp(mu_m*b)*I3nm, 0                      ];
        Amat(4*(n-1)+6,4*(m-1)+(3:6))=[ 0                   , mu_m*exp(mu_m*b)*I3nm  , -mu_m*exp(-mu_m*b)*I3nm, kappa(m)*(I0mn+I1mn)  ];
    end
    I300=fun_diff_Integ3mn(0,0,kappa0,0,D,draft);
    I30n=fun_diff_Integ3mn(0,n,kappa0,mu_n,D,draft);
    
    if n==1
        Bvect(1)=I300;
        Bvect(2)=I300;
    end
   
    Bvect(4*(n-1)+3)=I30n;
    Bvect(4*(n-1)+4)=I30n;
%     Bvect(4*(n-1)+5)=0*(I0mn+I1mn);
%     Bvect(4*(n-1)+6)=0*(I0mn+I1mn);
end


end