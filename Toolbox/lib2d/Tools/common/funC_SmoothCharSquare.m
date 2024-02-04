function cfun=funC_SmoothCharSquare(X,Y,X1,X2,Y1,Y2,del)

xl   =X(funC_closest(X,X1+del.l)); xr   =X(funC_closest(X,X2-del.r));
yb   =Y(funC_closest(Y,Y1+del.b)); yt   =Y(funC_closest(Y,Y2-del.t));


cfXl = funC_Heaviside(X-X1).*max(sign(X-xl),(1-cos((X-X1)*pi/del.l))/2);
cfXr = funC_Heaviside(X-xr).*max(sign(X-X2),(1-cos((X-xr)*pi/del.r))/2);

if del.l==0, cfXl=1;
elseif del.r==0, cfXr=0;end;
cfX=cfXl-cfXr;
cfX2d=repmat(cfX,length(Y),1);

cfYb = funC_Heaviside(Y-Y1).*max(sign(Y-yb),(1-cos((Y-Y1)*pi/del.b))/2);
cfYt = funC_Heaviside(Y-yt).*max(sign(Y-Y2),(1-cos((Y-yt)*pi/del.t))/2);
if del.b==0, cfYb=1;
elseif del.t==0, cfYt=0;end;
cfY=cfYb-cfYt;

cfY2d=repmat(cfY',1,length(X));

cfun=cfX2d.*cfY2d;