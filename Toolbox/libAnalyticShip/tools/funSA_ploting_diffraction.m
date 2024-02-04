    g=par.g;AreaO=par.AreaO;B=ship.width;T=ship.draft;Ainc=par.Ainc;kAinc=par.kAinc;D=bottom.depth;b=B/2;
    if calc.normVartype==1
    factx=g*AreaO.*kAinc;factz=g*B.*Ainc;facttheta=g*B^3.*kAinc/12;
    elseif calc.normVartype==2
    factx=g*T.*Ainc;factz=g*B.*Ainc;facttheta=g*B.*T*Ainc;
    elseif calc.normVartype==3
    factx=g*T.*Ainc;factz=g*b.*Ainc;facttheta=g*b.*T*Ainc;
    elseif calc.normVartype==4
    factx=g*B.*Ainc;factz=g*B.*Ainc;facttheta=g*b.*B*Ainc;
    elseif calc.normVartype==5
    factx=g*D.*Ainc;factz=g*D.*Ainc;facttheta=g*D^2*Ainc;
    end
    

    datCalForces= OutdiffRadd.diff.Forces;
    datCalReflTrans= OutdiffRadd.diff.ReflTrans;
    save([proj.workdir,'\datCalForces.mat'],'datCalForces');
    save([proj.workdir,'\datCalReflTrans.mat'],'datCalReflTrans');
   
    
    wwnorm=wave.input;
    Fx=datCalForces(:,2)./factx;
    Fz=datCalForces(:,3)./factz;
    Ftheta=datCalForces(:,4)./facttheta;    
    Rc=datCalReflTrans(:,2);Tc=datCalReflTrans(:,3);
    
    if wave.inputVar==1
      xlabstr=['$\omega \sqrt{B/2g}$'];  
    elseif wave.inputVar==2
       xlabstr=['$\omega^2 B/2g$'];
    elseif wave.inputVar==3
       xlabstr=['$\lambda/D$']; 
    elseif wave.inputVar==4
       xlabstr=['$\lambda/B$'];    
    elseif wave.inputVar==5
       xlabstr=['$kD$'];     
    elseif wave.inputVar==6
       xlabstr=['$kB$'];   
     elseif wave.inputVar==7
       xlabstr=['$kb$'];   
    elseif wave.inputVar==8
         xlabstr=['$kT$']; 
    end
   
   if calc.normVartype==1
      ylabFx=['$|F_x|/\rho gAka$'];
      ylabFz='$|F_z|/\rho gBa$';
      ylabFtheta=['$|M_{\theta}|/(\rho gB^3ka/12)$'];
   elseif calc.normVartype==2
        ylabFx=['$|F_x|/\rho gTA$'];
       ylabFz='$|F_z|/\rho gBA$';
       ylabFtheta=['$|M_{\theta}|/(\rho gBTa)$'];
   elseif calc.normVartype==3
       ylabFx=['$|F_x|/\rho gTA$'];
       ylabFz='$|F_z|/\rho gba$';
       ylabFtheta=['$|M_{\theta}|/(\rho gbTa)$'];   
    elseif calc.normVartype==4
       ylabFx=['$|F_x|/\rho gBa$'];
       ylabFz='$|F_z|/\rho gBa$';
       ylabFtheta=['$2|M_{\theta}|/(\rho gB^2a)$'];   
   elseif calc.normVartype==5
       ylabFx=['$|F_x|/\rho gDa$'];
       ylabFz='$|F_z|/\rho gDa$';
       ylabFtheta=['$|M_{\theta}|/(\rho gD^2a)$'];        
   end

    if isempty(calc.valdata)
        figure;
        subplot(3,1,1)
        plot(wwnorm,Fx,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFx,'interpreter','latex');
        title('Horizontal Excitation force');
        plot_properties;
        subplot(3,1,2)
        plot(wwnorm,Fz,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFz,'interpreter','latex');
        title('Vertical Excitation force');
        % xlim([0 2]);ylim([0 1.5]);
        plot_properties;
        
        subplot(3,1,3)
        plot(wwnorm,Ftheta,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFtheta,'interpreter','latex');
        title('Excitation moment');
        % xlim([0 2]);%ylim([0 0.2]);
        plot_properties;
        figure;
        plot(wwnorm,Rc,'r',wwnorm,Tc,'--r',wwnorm,Rc.^2+Tc.^2,':r')
        title('Reflection & Transmission coefficients');
        legend('R','T','R^2+T^2')
        xlabel(xlabstr,'Interpreter','latex');
        plot_properties;
        
    else
        datval=calc.valdata;
        datval(datval(:,1)==0,1)=NaN;
        datval(datval(:,3)==0,3)=NaN;
        datval(datval(:,5)==0,5)=NaN;
        RTData=1;
        if size(datval,2)==10
        datval(datval(:,7)==0,7)=NaN;
        datval(datval(:,9)==0,9)=NaN;
        else
        datval(:,7:10)=NaN;
        RTData=0;
        end
     
%         disp('here')
        
        figure;
        subplot(3,1,1)
        plot(datval(:,1),datval(:,2),'ob',wwnorm,Fx,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFx,'interpreter','latex');
        title('Horizontal Excitation force');
        legend('data','analytic')
        xlim([wwnorm(1) wwnorm(end)])
        plot_properties;
        subplot(3,1,2)
        plot(datval(:,3),datval(:,4),'ob',wwnorm,Fz,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFz,'interpreter','latex');
        title('Vertical Excitation force');
        % xlim([0 2]);ylim([0 1.5]);
          xlim([wwnorm(1) wwnorm(end)])
        plot_properties;
        
        subplot(3,1,3)
        plot(datval(:,5),datval(:,6),'ob',wwnorm,Ftheta,'r');
        xlabel(xlabstr,'Interpreter','latex');
        ylabel(ylabFtheta,'interpreter','latex');
        title('Excitation moment');
        % xlim([0 2]);%ylim([0 0.2]);
          xlim([wwnorm(1) wwnorm(end)])
        plot_properties;
        
        figure;
        title('Reflection & Transmission coefficients');
        if RTData==1
        plot(datval(:,7),datval(:,8),'*b',datval(:,9),datval(:,10),'ob',wwnorm,Rc,'r',wwnorm,Tc,'--r',wwnorm,Rc.^2+Tc.^2,':r')
        legend('R data','T data','R analytic','T analytic','R^2+T^2')
        else
         plot(wwnorm,Rc,'r',wwnorm,Tc,'--r',wwnorm,Rc.^2+Tc.^2,':r')   
          legend('R analytic','T analytic','R^2+T^2')   
        end
        xlabel(xlabstr,'Interpreter','latex');
          xlim([wwnorm(1) wwnorm(end)])
        plot_properties;
    end