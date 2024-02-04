function funP_contour(axesfig,XX,YY,ZZ,N,Zlev,col,lw)
tcon = hgtransform('Parent',axesfig);
if ~isempty(N)
[~,h]=contour(XX,YY,ZZ,N,'parent',tcon);
else
[~,h]=contour(XX,YY,ZZ,'parent',tcon);    
end
if ~isempty(col)
set(h,'linecolor',col,'linewidth',lw)    
end
set(tcon, 'Matrix',makehgtform('translate',[0 0 Zlev]));