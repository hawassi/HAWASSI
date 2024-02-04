%%2d characteristic function
function ChiAdj=funC_cfNonLinAdj2d(dom,lambda_p,spatinf,Infdirection,adjcoef)

if adjcoef~=0
    inflPosX=spatinf.x; inflPosY=spatinf.y;
    dXinf=max(inflPosX)-min(inflPosX);
    dYinf=max(inflPosY)-min(inflPosY);
    if dXinf<dYinf
        LadjLR=(round(adjcoef*lambda_p/dom.dx)+1)*dom.dx;
        inflPosX=[inflPosX(1);inflPosX;inflPosX(end)];
        inflPosY=[inflPosY(1)-LadjLR;inflPosY;inflPosY(end)+LadjLR];
        
        inflPosX1=inflPosX-LadjLR;
        inflPosX2=inflPosX+LadjLR;
        inflPosXN=[inflPosX1;flipud(inflPosX2)];
        inflPosYN=[inflPosY;flipud(inflPosY)];
       
        mappedindex=funC_mappedmatrixindexes(dom.X,dom.Y,inflPosXN,inflPosYN);
        ChiAdjBox=1-inpolygon(dom.XX,dom.YY,dom.XX(mappedindex),dom.YY(mappedindex));
       
        ChiAdjtemp(1,:,:)=ChiAdjBox;
        Nsmooth=max(1,floor(LadjLR/dom.dx));
        if mod(Nsmooth,2)==0, Nsmooth=Nsmooth-1;end
         if Nsmooth>29
            Nsmooth=max(29,floor(LadjLR/dom.dx/4));
         end
        ChiAdjtemp=smooth3(ChiAdjtemp,'box',Nsmooth);
        ChiAdj=squeeze(ChiAdjtemp(1,:,:));
       
    else
        LadjBT=(round(adjcoef*lambda_p/dom.dy)+1)*dom.dy;
        inflPosY=[inflPosY(1);inflPosY;inflPosY(end)];
        inflPosX=[inflPosX(1)-LadjBT;inflPosX;inflPosX(end)+LadjBT];
        inflPosY1=inflPosY-LadjBT;
        inflPosY2=inflPosY+LadjBT;
        inflPosXN=[inflPosX;flipud(inflPosX)];
        inflPosYN=[inflPosY1;flipud(inflPosY2)];
        mappedindex=funC_mappedmatrixindexes(dom.X,dom.Y,inflPosXN,inflPosYN);
        ChiAdjBox=1-inpolygon(dom.XX,dom.YY,dom.XX(mappedindex),dom.YY(mappedindex));
        ChiAdjtemp(1,:,:)=ChiAdjBox;
        Nsmooth=max(1,floor(LadjBT/dom.dy));
        if mod(Nsmooth,2)==0, Nsmooth=Nsmooth-1;end
        if Nsmooth>29
            Nsmooth=max(29,floor(LadjBT/dom.dy/4));
         end
        ChiAdjtemp=smooth3(ChiAdjtemp,'box',Nsmooth);
        ChiAdj=squeeze(ChiAdjtemp(1,:,:));
    end
    
    
    %     for jj=1:length(inflPosY)
    %         indyinfl=funC_closest(dom.Y,inflPosY(jj));
    %             ChiAdjL    =funC_cf(dom.X,inflPosX(jj),LadjLR);
    %             ChiAdjR    =1-funC_cf(dom.X,inflPosX(jj)-LadjLR,LadjLR);
    %             ChiAdj(indyinfl,:)    = ChiAdj(indyinfl,:).*(ChiAdjR+ChiAdjL);
    %     end
    %     LadjBT=(round(adjcoef*lambda_p/dom.dy)+1)*dom.dy;
    %     for jj=1:length(inflPosX)
    %            indxinfl=funC_closest(dom.X,inflPosX(jj));
    %             ChiAdjB    =funC_cf(dom.Y,inflPosY(jj),LadjBT)';
    %             ChiAdjT    =1-funC_cf(dom.Y,inflPosY(jj)-LadjBT,LadjBT)';
    %             ChiAdj(:,indxinfl)=ChiAdj(:,indxinfl).*(ChiAdjB+ChiAdjT);
    %     end
    %
    %     ChiAdjtemp(1,:,:)=ChiAdj;
    %     ChiAdjtemp=smooth3(ChiAdjtemp,'box',5);
    %     ChiAdj=squeeze(ChiAdjtemp(1,:,:));
    %  pause
else
    ChiAdj=ones(dom.Ny,dom.Nx);
end
%
%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(2)
%    pause;
%     figure
%     [X,Y]=meshgrid(dom.X,dom.Y);
%     mesh(X,Y,ChiAdj);
%     view(3)
