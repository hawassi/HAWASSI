%% findex indexes for mapping profile in mattrix x,y to matrix X,Y

function index=funC_mappedmatrixindexes(X,Y,x,y)
Np=length(x);
index=ones(Np,1);
[XX,YY]=meshgrid(X,Y);
for i=1:Np
 dist=sqrt((XX-x(i)).^2+(YY-y(i)).^2);
[~,index(i)]=min(reshape(dist,1,[]));
end




