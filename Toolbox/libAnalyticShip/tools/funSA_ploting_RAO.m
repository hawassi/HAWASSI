
  datCalRAO=OutRAO.notcoupled;
  datCalRAOCM=OutRAO.coupled;
  g=par.g;B=ship.width;
  wwnorm=wave.input;
  
  
save([proj.workdir,'\datCalRAO.mat'],'datCalRAO');
save([proj.workdir,'\datCalRAOCM.mat'],'datCalRAOCM');
  
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
  
  if isempty(calc.valdata)
    flagVdat=0;
else
    flagVdat=1;
    datval=calc.valdata;
    datval(datval(:,1)==0,1)=NaN;
    datval(datval(:,3)==0,3)=NaN;
    datval(datval(:,5)==0,5)=NaN;
  end

 figure;
 subplot(2,1,1)
 if flagVdat==0
 plot(wwnorm,datCalRAO(:,2),'r',wwnorm,datCalRAO(:,3),'g',wwnorm,datCalRAO(:,4),'c');
 legend({'$|\xi_2/a|$','$|\xi_3/a| $','$|\xi_4/ka|$'},'Interpreter','latex');
 else
 plot(datval(:,1),datval(:,2),'-ob',datval(:,3),datval(:,4),'-db',datval(:,5),datval(:,6),'-*b',wwnorm,datCalRAO(:,2),'r',wwnorm,datCalRAO(:,3),'g',wwnorm,datCalRAO(:,4),'c');    
legend({'$|\xi_2/a|$ data','$|\xi_3/a|$ data','$|\xi_4/ka|$ data','$|\xi_1/a|$ Analytic','$|\xi_3/a|$ Analytic','$|\xi_5/ka|$ Analytic'},'Interpreter','latex');
 
 end
 title('RAO for uncoupled motions, no visc damping applied')
 xlabel(xlabstr,'Interpreter','latex');
 ylabel('RAO','Interpreter','latex');
 %ylim([0 4]);
 plot_properties;
 subplot(2,1,2)
 if flagVdat==0
     plot(wwnorm,datCalRAOCM(:,2),'r',wwnorm,datCalRAOCM(:,3),'g',wwnorm,datCalRAOCM(:,4),'c');

 else
     plot(datval(:,1),datval(:,2),'-ob',datval(:,3),datval(:,4),'-db',datval(:,5),datval(:,6),'-*b',wwnorm,datCalRAOCM(:,2),'r',wwnorm,datCalRAOCM(:,3),'g',wwnorm,datCalRAOCM(:,4),'c');
 end
      title('RAO for coupled motions')
 xlabel(xlabstr,'Interpreter','latex');
 ylabel('RAO','Interpreter','latex');
 %ylim([0 4]);
 plot_properties;
 