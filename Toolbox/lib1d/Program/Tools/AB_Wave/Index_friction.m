function [index,Cf]=Index_friction(x,cfdata)
Ni=length(cfdata.interval);
dx=x(2)-x(1);
Cf=zeros(Ni,1);
Iter=1;
for i=1:Ni
    Interv=str2num(cfdata.interval{i});
    xinterv=[Interv(1):dx:Interv(2)];
    Nxnow=length(xinterv);
    index(Iter:Iter+Nxnow-1)=closest(x,xinterv);
    Cf(Iter:Iter+Nxnow-1,1)=str2num(cfdata.cf{i});
    Iter=Iter+Nxnow;
end