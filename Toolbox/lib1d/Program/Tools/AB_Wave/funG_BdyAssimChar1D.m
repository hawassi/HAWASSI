function [charupdate,ChiAdj,propdir]=funG_BdyAssimChar1D(par,assimdat)
Xleft=assimdat.data.xinterv(1);
Xright=assimdat.data.xinterv(2);
bf0=assimdat.data.smoothfact*par.dx;
propdir=assimdat.data.direction;
charupdate    = cf(par.x,Xleft,bf0)-cf(par.x,Xright-bf0,bf0); % cfSimulationArea
ChiAdj=1;
if assimdat.data.checkbox_nonlin_adj==1
    LAdj=assimdat.data.nonlinadj_length;  
    if propdir==1
        ChiAdj=cfAdj(par.x,Xright,LAdj,1,'Uni+');
    else
        ChiAdj=cfAdj(par.x,Xleft,LAdj,1,'Uni-');
    end
end

end