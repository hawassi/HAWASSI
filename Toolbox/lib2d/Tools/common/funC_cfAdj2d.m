%%2d characteristic function
function ChiAdj=funC_cfAdj2d(dom,lambda_p,spatinf,Infdirection,adjcoef)

if adjcoef~=0
    inflPosX=spatinf.x; inflPosY=spatinf.y;
    ChiAdj=ones(dom.Ny,dom.Nx);
    for jj=1:length(inflPosY)
        indyinfl=funC_closest(dom.Y,inflPosY(jj));
        %         if strcmp(Infdirection,'UniR')
        %             ChiAdj(indyinfl,:)   = ChiAdj(indyinfl,:).*(funC_cf(dom.X,inflPosX(jj),adjcoef*lambda_p));
        %         elseif strcmp(Infdirection,'UniL')
        %             ChiAdj(indyinfl,:)   = ChiAdj(indyinfl,:).*(1-funC_cf(dom.X,inflPosX(jj)-adjcoef*lambda_p,adjcoef*lambda_p));
        %         elseif strcmp(Infdirection,'BiRL')
        %             ChiAdjL    =funC_cf(dom.X,inflPosX(jj),adjcoef*lambda_p);
        %             ChiAdjR    =1-funC_cf(dom.X,inflPosX(jj)-adjcoef*lambda_p,adjcoef*lambda_p);
        %             ChiAdj(indyinfl,:)    = ChiAdj(indyinfl,:).*(ChiAdjR+ChiAdjL);
        %         end
        if strcmp(Infdirection,'UniR') || strcmp(Infdirection,'UniL')
            ChiAdjL    =funC_cf(dom.X,inflPosX(jj),adjcoef*lambda_p);
            ChiAdjR    =1-funC_cf(dom.X,inflPosX(jj)-adjcoef*lambda_p,adjcoef*lambda_p);
            ChiAdj(indyinfl,:)    = ChiAdj(indyinfl,:).*(ChiAdjR+ChiAdjL);
        end
    end
    
    for jj=1:length(inflPosX)
        indxinfl=funC_closest(dom.X,inflPosX(jj));
%         if strcmp(Infdirection,'UniT')
%             ChiAdj(:,indxinfl)   = ChiAdj(:,indxinfl).*funC_cf(dom.Y,inflPosY(jj),adjcoef*lambda_p)';
%         elseif strcmp(Infdirection,'UniB')
%             ChiAdj(:,indxinfl)   = ChiAdj(:,indxinfl).*(1-funC_cf(dom.Y,inflPosY(jj)-adjcoef*lambda_p,adjcoef*lambda_p)');
%         elseif strcmp(Infdirection,'BiTB')
%             ChiAdjB    =funC_cf(dom.Y,inflPosY(jj),adjcoef*lambda_p)';
%             ChiAdjT    =1-funC_cf(dom.Y,inflPosY(jj)-adjcoef*lambda_p,adjcoef*lambda_p)';
%             ChiAdj(:,indxinfl)=ChiAdj(:,indxinfl).*(ChiAdjB+ChiAdjT);
%         end
        if strcmp(Infdirection,'UniT') || strcmp(Infdirection,'UniB')
            ChiAdjB    =funC_cf(dom.Y,inflPosY(jj),adjcoef*lambda_p)';
            ChiAdjT    =1-funC_cf(dom.Y,inflPosY(jj)-adjcoef*lambda_p,adjcoef*lambda_p)';
            ChiAdj(:,indxinfl)=ChiAdj(:,indxinfl).*(ChiAdjB+ChiAdjT);
        end
    end
else
    ChiAdj=ones(dom.Ny,dom.Nx);
end

%
%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(2)

%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(3)
