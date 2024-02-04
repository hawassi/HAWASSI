function [AmatW,BvectW] =fun_construct_matrices_reflecting_wall(D,Evmodes,Kappa,alphaR,phiL,uL)
N=1+Evmodes;
AmatW=zeros(N,N);BvectW=zeros(N,1);

kappa0=Kappa(1);
kappa=Kappa(2:end);

for n=0:Evmodes
    if n>0
        IntQ0Qn=fun_IntegQmQn(0,n,kappa0,kappa(n),0,D,0);
    else
        IntQ0Qn=fun_IntegQmQn(0,n,kappa0,kappa0,0,D,0);
    end
   AmatW(n+1,1)=[-kappa0*(1+alphaR)*IntQ0Qn];
   
   if Evmodes>=1
       for m=1:Evmodes
           if n>0
           IntQmQn=fun_IntegQmQn(m,n,kappa(m),kappa(n),0,D,0);
           else
            IntQmQn=fun_IntegQmQn(m,n,kappa(m),kappa0,0,D,0);     
           end
           AmatW(n+1,m+1)=[(kappa(m)-alphaR*kappa0)*IntQmQn];
       end
   end
   BvectW(n+1)=(alphaR*phiL*kappa0-uL)*IntQ0Qn;%%(alphaR*phiL*kappa0-uL)*IntQ0Qn;%%
end
% 
% N=Evmodes;
% AmatW=zeros(N,N);BvectW=zeros(N,1);
% 
% kappa0=Kappa(1);
% kappa=Kappa(2:end);
% 
% for n=1:Evmodes
%     
%    if Evmodes>=1
%        for m=1:Evmodes
%            IntQmQn=fun_IntegQmQn(m,n,kappa(m),kappa(n),0,D,0);
%            AmatW(n,m)=[(kappa(m)-alphaR*kappa0)*IntQmQn];
%        end
%    end
%        IntQ0Qn=fun_IntegQmQn(0,n,kappa0,kappa(n),0,D,0);
% 
%    
%    BvectW(n)=(alphaR*phiL*kappa0-uL)*IntQ0Qn;%%uL*(alphaR-1)*IntQ0Qn;
% end