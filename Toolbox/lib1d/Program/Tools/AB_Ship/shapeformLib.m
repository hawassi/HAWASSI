function [shape,dxshape,shapedat]=shapeformLib(x,Sxc,Szc,Sthetac,T,ell,betaL,shipform)
draft=T-Szc;
shapedat=0;
switch shipform
    case 'cfGauss'
        %%
        shape   = - 2*draft*(exp(-(2*x/ell).^2)-1/2);
        % dxshape = 4/(ell)^2*x.*(2*draft*(exp(-(2*x/ell).^2)));
       % dxshape = gradient(shape,x(2)-x(1));
        %         figure(1)
        %         plot(x,shape,'b',x,dxshape,'r:')
    case 'U4'
        %%
        shape = min(draft/ell^4*(x.^2-ell^2).*(x.^2+ell^2),4);
        % dxshape = draft/ell^4*x.^3;%.*chiship;
     %   dxshape = gradient(shape,x(2)-x(1));
        %         figure(1)
        %         plot(x,shape,'b',x,dxshape,'r:')
        %%
    case 'Barge'
        %%
        
        %         dx=x(2)-x(1);
        %         ell=(ell*2+dx)/2;
        %         ind0=closest(x,Sxc);
        %          shape=zeros(size(x));
        %          indl=closest(x,x(ind0)-ell)+1;
        %          indr=ind0+ind0-indl;
        %          shape(indl:indr)=-draft;
        b=ell;
        Np=159;
        dX=2*b/Np;
        dZ=draft/Np;
        ZSR=[linspace(-draft-dZ,-draft,50) linspace(-draft+dZ,0,Np)];
        
        XSR=zeros(size(ZSR));
        IndR=closest(ZSR,-draft);
        XSR(1:IndR)=(ZSR(1:IndR)-ZSR(1))*(b-dX)/dZ;
        XSR(IndR:end)=(ZSR(IndR:end)-ZSR(IndR))*dX/(draft)+XSR(IndR);
        ZSL=[linspace(-draft-dZ,-draft,50) linspace(-draft+dZ,0,Np)];
       
        XSL=zeros(size(ZSL));
        IndL=closest(ZSL,-draft);
        XSL(1:IndL)=-(ZSL(1:IndL)-ZSL(1))*(b-dX)/dZ;
        XSL(IndL:end)=-((ZSL(IndL:end)-ZSL(IndL))*dX/(draft)-XSL(IndL));
        shapedat=[Sxc+sort(XSL.') sort(ZSL.','descend' );Sxc+XSR(2:end).' ZSR(2:end).'];
        X=shapedat(:,1);Z=shapedat(:,2);
        [shapedat(:,1), shapedat(:,2)]=fun_rotation(X,Z,Sthetac,Sxc);
       
       
%         figure;
%         plot(shapedat0(:,1),shapedat0(:,2),'b',shapedat(:,1),shapedat(:,2),'r')
        
        shape=zeros(size(x));
        indx1=closest(x,min(shapedat(:,1)));
        indx2=closest(x,max(shapedat(:,1)));
       
        if x(indx1)<shapedat(1,1)
            indx1=indx1+1;
        end
        
        if x(indx2)>shapedat(end,1)
            indx2=indx2-1;
        end
        
        shape(indx1:indx2)=interp1(shapedat(:,1),shapedat(:,2),x(indx1:indx2));
        dx=x(2)-x(1);indl=indx1;indr=indx2;
        DDL =zeros(1,length(x)); DDL(indl-1)=draft/dx;
        DDR =zeros(1,length(x));  DDR(indr+1) =draft/dx;
        dxshape0 = DDR-DDL;dxshapeGrad=gradient(shape,dx);
        
%                 figure
%                 plot(shapedat(:,1),shapedat(:,2))
        
        %         dxshape = dxshape';
        %         figure(1)
        %         plot(x,shape,'r',x(ind0),0,'+r',x(indl),0,'+r',x(indr),0,'+r');%,x,chiship,'k:',xL,0,'c*',xR,0,'c*',x,dxshape,'r:')
        %  xlabel('x [m]') ; ylabel('\eta [m]');
        %                 grid on;
        % title(['Shipform ',shipform]);
        %            saveas(gcf,[sf_savename,'Neumanncond-',RBdyn,'.fig'])
        %%
    case 'Wedge'
        zT0min=linspace(-draft,0,100);
        X_min=-ell.*(zT0min./draft+1);
        zT0plus=linspace(-draft,0,100);
        X_plus=ell.*(zT0plus./draft+1);
        shapedat=[Sxc+sort(X_min.') sort(zT0min.','descend' ); Sxc+X_plus(2:end).' zT0plus(2:end).'];
        shape=zeros(size(x));
        indx1=closest(x,shapedat(1,1));
        indx2=closest(x,shapedat(end,1));
        
        if x(indx1)<shapedat(1,1)
            indx1=indx1+1;
        end
        
        if x(indx2)>shapedat(end,1)
            indx2=indx2-1;
        end
        
        shape(indx1:indx2)=interp1(shapedat(:,1),shapedat(:,2),x(indx1:indx2));
       % dxshape=gradient(shape,x(2)-x(1));
        figure;
        plot(x,shape,'r',shapedat(:,1),shapedat(:,2),'--r');xlim([-200; 200])
        
        
        %%
    case 'Half circle'
        zT0min=linspace(-draft,0,200);
        theta=acos(zT0min/draft);
        X_min=-ell*sqrt(1-cos(theta).^2);
        zT0plus=linspace(-draft,0,200);
        theta=acos(zT0plus/draft);
        X_plus=ell*sqrt(1-cos(theta).^2);
        shapedat=[Sxc+sort(X_min.') sort(zT0min.','descend' ); Sxc+X_plus(2:end).' zT0plus(2:end).'];
        shape=zeros(size(x));
        indx1=closest(x,shapedat(1,1));
        indx2=closest(x,shapedat(end,1));
        
        if x(indx1)<shapedat(1,1)
            indx1=indx1+1;
        end
        
        if x(indx2)>shapedat(end,1)
            indx2=indx2-1;
        end
        
        
        shape(indx1:indx2)=interp1(shapedat(:,1),shapedat(:,2),x(indx1:indx2));
   %     dxshape=gradient(shape,x(2)-x(1));
        %     figure;
        %     plot(x,shape,'r',shapedat(:,1),shapedat(:,2),'--r');xlim([-200; 200])
        
    case 'smsquare'
        %%
        del = ell/10;
        shape = -2*draft*(cf(x,-ell-del/2,del)-cf(x,ell-del/2,del))+draft;
   %     dxshape = gradient(shape,x(2)-x(1));
        %         shape = shape.*(cf(x,-ell,del)).*(1-cf(x,ell-del,del));
        %         dxshape = Ifft(1i*k.*fft(shape)')';
        figure(1)
        plot(x,shape,'b',x,dxshape,'r:')
        %%
    case 'Lewis'
        %%
        B=ell*2;
        T=draft;
        theta=[0.5*pi:0.1:1.5*pi];
        H=B/(2*T);
        %betaL=A/(B*T);
        p=betaL-pi/4;
        q=(H-1)/(H+1);
        b=(0.75*pi+sqrt((0.25*pi)^2-0.5*pi.*p*(1-q^2)))./(pi+p.*(1-q^2))-1;
        a=(b+1)*q;
        figure(1)
        for i=1:length(betaL)
            y=((1+a(i))*sin(theta)-b(i)*sin(3.*theta))*(B/(2.*(1+a(i)+b(i))));
            z=((1-a(i))*cos(theta)+b(i)*cos(3*theta))*(B/(2*(1+a(i)+b(i))));
            plot(y,z,'-b');hold on;
            xlim([-ell ell]); ylim([-draft 0])
        end
        shape=zeros(size(x));
        indx1=closest(x,-ell);indx2=closest(x,ell);
        shape(indx1:indx2)=interp1(y,z,x(indx1:indx2),'spline');
     %   dxshape=gradient(shape,x(2)-x(1));
        figure;
        plot(x,shape)
        xlim([x(1) x(end)]); ylim([-draft 0])
        %pause
        
end

dxshape=funOprt_FDGradient1dShip(shape,x,Sxc);

% figure;
% subplot(2,1,1)
% plot(x,shape);
% subplot(2,1,2)
% plot(x,dxshape0,'b',x,dxshape,'--r')


if shapedat==0
    shapedat=zeros(length(x),2);
    shapedat(:,1)=x;
    shapedat(:,2)=shape;
end