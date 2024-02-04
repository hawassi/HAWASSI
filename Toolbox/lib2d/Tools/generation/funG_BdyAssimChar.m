function [charupdate,charnonlinAdj,dir]=funG_BdyAssimChar(dom,bdyassim)
    
    if bdyassim.shapeOpt==1 %halfcircle
        r1=bdyassim.halfcirc_R1;
        xc=bdyassim.halfcirc_xc;yc=bdyassim.halfcirc_yc;      
        smoothfact=bdyassim.smoothfact;
        distNonlinAdj=bdyassim.nonlinAdj_distance;
        smoothNonlinAdj=bdyassim.nonlinAdj_smooth;
        if bdyassim.propdir==1 %South
            dir='South';distNonlinAdj=-distNonlinAdj;
        elseif bdyassim.propdir==2 %North
            dir='North';
        elseif  bdyassim.propdir==3 %West
             dir='West';distNonlinAdj=-distNonlinAdj;
        elseif bdyassim.propdir==4 %East
            dir='East';
        end
        X1=dom.X(1)+dom.fbl.l;X2=dom.X(end)-dom.fbl.r;
        Y1=dom.Y(1)+dom.fbl.b;Y2=dom.Y(end)-dom.fbl.t;
        dx=dom.X(2)-dom.X(1); dy=dom.Y(2)-dom.Y(1);
        del.l=smoothfact*dx;del.r=smoothfact*dx;
        del.b=smoothfact*dy;del.t=smoothfact*dy;
        
        
        charupdate=funC_char_circle_smooth(dom.X,dom.Y,xc,yc,r1,smoothfact,dir).*...
           funC_SmoothCharSquare(dom.X,dom.Y,X1,X2,Y1,Y2,del);
        
        
        if bdyassim.cb_nonlinAdj==1
        charnonlinAdj=1-funC_char_circle_smooth(dom.X,dom.Y,xc,yc,r1+distNonlinAdj,smoothNonlinAdj,dir);
        else
         charnonlinAdj=ones(size(dom.XX));   
        end
    else
    end
    
    
%                        figure(111)
%                set(gcf,'Renderer','zbuffer'); %due to graphics driver
%                 surf(dom.XX,dom.YY,charupdate,'edgecolor','none')
%         

% figure(111)
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% subplot(2,1,1)
% surf(dom.XX,dom.YY,charupdate,'edgecolor','none');
% subplot(2,1,2)
% surf(dom.XX,dom.YY,charnonlinAdj,'edgecolor','none');

        
end

