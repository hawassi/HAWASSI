function shipprof=funSA_shiplib(ship)
B=ship.width;
b=B/2;
draft=ship.draft;

if ship.shape==1
    Np=159;
    dX=2*b/Np;
    dZ=draft/Np;
    ZSR=linspace(-draft-dZ,0,Np);
    XSR=zeros(size(ZSR));
    IndR=closest(ZSR,-draft);
    XSR(1:IndR)=(ZSR(1:IndR)-ZSR(1))*(b-dX)/dZ;
    XSR(IndR:end)=(ZSR(IndR:end)-ZSR(IndR))*dX/(draft)+XSR(IndR);
    ZSL=linspace(-draft-dZ,0,Np);
    XSL=zeros(size(ZSL));
    IndL=closest(ZSL,-draft);
    XSL(1:IndL)=-(ZSL(1:IndL)-ZSL(1))*(b-dX)/dZ;
    XSL(IndL:end)=-((ZSL(IndL:end)-ZSL(IndL))*dX/(draft)-XSL(IndL));
    shipprof=[sort(XSL.') sort(ZSL.','descend' ); XSR(2:end).' ZSR(2:end).'];
    
elseif ship.shape==2
    zT0min=linspace(-draft,0,200);
    theta=acos(zT0min/draft);
    X_min=-b*sqrt(1-cos(theta).^2);
    zT0plus=linspace(-draft,0,200);
    theta=acos(zT0plus/draft);
    X_plus=b*sqrt(1-cos(theta).^2);
    %      figure;
    %     plot(X_min,zT0min,'r',X_plus,zT0plus,'--b')
    
    shipprof=[sort(X_min.') sort(zT0min.','descend' ); X_plus(2:end).' zT0plus(2:end).'];
    % title(['Area=',num2str(Area),' B*T=',num2str(b*draft*2)])
    
elseif ship.shape==3
    zT0min=linspace(-draft,0,100);
    X_min=-b.*(zT0min./draft+1);
    zT0plus=linspace(-draft,0,100);
    X_plus=b.*(zT0plus./draft+1);
    shipprof=[sort(X_min.') sort(zT0min.','descend' ); X_plus(2:end).' zT0plus(2:end).'];

else
    shipprof=ship.shapedat;
end
% datshipprof=shipprof;
% save datshipprof datshipprof

end