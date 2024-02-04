function cfun=funC_tapered2d(X,Y,xs,xf,ys,yf,del) 
xl   =X(funC_closest(X,xs+del.l)); xr   =X(funC_closest(X,xf-del.r));
yb   =Y(funC_closest(Y,ys+del.b)); yt   =Y(funC_closest(Y,yf-del.t));


cfXl = funC_Heaviside(X-xs).*max(sign(X-xl),(1-cos((X-xs)*pi/del.l))/2);
cfXr = funC_Heaviside(X-xr).*max(sign(X-xf),(1-cos((X-xr)*pi/del.r))/2);
if del.l==0, cfXl=1;
elseif del.r==0, cfXr=0;end
cfX=cfXl-cfXr;
cfX2d=repmat(cfX,length(Y),1);

cfYb = funC_Heaviside(Y-ys).*max(sign(Y-yb),(1-cos((Y-ys)*pi/del.b))/2);
cfYt = funC_Heaviside(Y-yt).*max(sign(Y-yf),(1-cos((Y-yt)*pi/del.t))/2);
if del.b==0, cfYb=1;
elseif del.t==0, cfYt=0;end
cfY=cfYb-cfYt;

cfY2d=repmat(cfY',1,length(X));


cfun=cfX2d.*cfY2d;

