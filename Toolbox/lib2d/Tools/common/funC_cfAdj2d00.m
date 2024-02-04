%%2d characteristic function
function ChiAdj=funC_cfAdj2d00(dom,lambda_p,spatinf,Infdirection,adjcoef)

  if adjcoef~=0 
    inflPos.x=0;%spatinf.x(1);
    if strcmp(Infdirection,'UniR')
        ChiAdj   = funC_cf2d(dom.X,inflPos.x,adjcoef*lambda_p,length(dom.Y));
     elseif strcmp(Infdirection,'UniL')
        ChiAdj   = 1-funC_cf2d(dom.X,inflPos.x-adjcoef*lambda_p,adjcoef*lambda_p,length(dom.Y));
     elseif strcmp(Infdirection,'Bi')
        ChiAdjL    =funC_cf2d(dom.X,inflPos.x,adjcoef*lambda_p.y,length(dom.Y));
        ChiAdjR    =1-funC_cf2d(dom.X,inflPos.x-adjcoef*lambda_p,adjcoef*lambda_p.y,length(dom.Y));
        ChiAdj    = ChiAdjR+ChiAdjL;
    end

    inflPos.y=0;%spatinf.y(1);
    if strcmp(Infdirection,'UniT')
        ChiAdj   = funC_cf2d(dom.Y,inflPos.y,adjcoef*lambda_p,length(dom.X))';
    elseif strcmp(Infdirection,'UniB')
        ChiAdj   = 1-funC_cf2d(dom.Y,inflPos.y-adjcoef*lambda_p,adjcoef*lambda_p,length(dom.X))';
    elseif strcmp(Infdirection,'Bi')
        ChiAdjB    =1-funC_cf2d(dom.Y,inflPos.y,adjcoef*lambda_p,length(dom.X))';
        ChiAdjT    =funC_cf2d(dom.Y,inflPos.y-adjcoef*lambda_p,adjcoef*lambda_p.x,length(dom.X))';
        ChiAdj   =ChiAdjB+ChiAdjT;
    end
  else
      ChiAdj=ones(dom.Ny,dom.Nx);
  end
    
%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(2)
    
%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(3)
