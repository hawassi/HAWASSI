datCalAddedMass=OutdiffRadd.rad.addedmass;
datCalDampCoef=OutdiffRadd.rad.dampCoef;%
datCalRadWaveAmpl=OutdiffRadd.rad.waveAmpl;%

save([proj.workdir,'\datCalAddedMass.mat'],'datCalAddedMass');
save([proj.workdir,'\datCalDampCoef.mat'],'datCalDampCoef');
save([proj.workdir,'\datCalRadWaveAmpl.mat'],'datCalRadWaveAmpl');
  
B=ship.width;g=par.g;AreaO=par.AreaO;D=bottom.depth;w0=par.w0;b=B/2;
wwnorm=wave.input;T=ship.draft;

 if calc.normVartype==3
    factA=D^2;factAB=D^3;factAB2=D^4;  
    factB=D^2/sqrt(D/g);factBB=D^3/sqrt(D/g);factBB2=D^4/sqrt(D/g); 
 elseif calc.normVartype==2
    factA=AreaO;factAB=AreaO*B;factAB2=AreaO*B*B;  
    factB=AreaO*w0;factBB=AreaO*B*w0;factBB2=AreaO*B*B*w0;
 else
    factA=AreaO;factAB=AreaO*B;factAB2=AreaO*B*B;  
    factB=AreaO/sqrt(B/2/g);factBB=AreaO*B/sqrt(B/2/g);factBB2=AreaO*B*B/sqrt(B/2/g);
 end
 
A11=datCalAddedMass(:,2)./factA;A13=datCalAddedMass(:,3)./factA;A15=datCalAddedMass(:,4)./factAB;
A31=datCalAddedMass(:,5)./factA;A33=datCalAddedMass(:,6)./factA;A35=datCalAddedMass(:,7)./factAB;
A51=datCalAddedMass(:,8)./factAB;A53=datCalAddedMass(:,9)./factAB;A55=datCalAddedMass(:,10)./factAB2;

B11=datCalDampCoef(:,2)./factB;B13=datCalDampCoef(:,3)./factB;B15=datCalDampCoef(:,4)./factBB;
B31=datCalDampCoef(:,5)./factB;B33=datCalDampCoef(:,6)./factB;B35=datCalDampCoef(:,7)./factBB;
B51=datCalDampCoef(:,8)./factBB;B53=datCalDampCoef(:,9)./factBB;B55=datCalDampCoef(:,10)./factBB2;

eta1=datCalRadWaveAmpl(:,2:3);
eta3=datCalRadWaveAmpl(:,4:5);
eta5=datCalRadWaveAmpl(:,6:7)*2/B;

if calc.normVartype==3
ylabA11=['$\mathcal{A}_{22}/\rho D^2$'];
ylabA13=['$\mathcal{A}_{23}/\rho D^2 $'];
ylabA15=['$\mathcal{A}_{24}/\rho D^3$'];
ylabA31=['$\mathcal{A}_{32}/\rho D^2$'];
ylabA33=['$\mathcal{A}_{33}/\rho D^2 $'];
ylabA35=['$\mathcal{A}_{34}/\rho D^3$'];
ylabA51=['$\mathcal{A}_{42}/\rho D^3$'];
ylabA53=['$\mathcal{A}_{43}/\rho D^3$'];
ylabA55=['$\mathcal{A}_{44}/\rho D^4$'];
ylabB11=['$(\mathcal{B}_{22}/\rho D^2)\sqrt{D/g}$'];
ylabB13=['$(\mathcal{B}_{23}/\rho D^2)\sqrt{D/g}$'];
ylabB15=['$(\mathcal{B}_{24}/\rho D^3)\sqrt{D/g}$'];
ylabB31=['$(\mathcal{B}_{32}/\rho D^2)\sqrt{D/g}$'];
ylabB33=['$(\mathcal{B}_{33}/\rho D^2)\sqrt{D/g}$'];
ylabB35=['$(\mathcal{B}_{34}/\rho D^3)\sqrt{D/g}$'];
ylabB51=['$(\mathcal{B}_{42}/\rho D^3)\sqrt{D/g}$'];
ylabB53=['$(\mathcal{B}_{43}/\rho D^3)\sqrt{D/g}$'];
ylabB55=['$(\mathcal{B}_{44}/\rho D^4)\sqrt{D/g}$'];
elseif calc.normVartype==2
ylabA11=['$\mathcal{A}_{22}/\rho A$'];
ylabA13=['$\mathcal{A}_{23}/\rho A$'];
ylabA15=['$\mathcal{A}_{24}/\rho AB$'];
ylabA31=['$\mathcal{A}_{32}/\rho A$'];
ylabA33=['$\mathcal{A}_{33}/\rho A $'];
ylabA35=['$\mathcal{A}_{34}/\rho AB$'];
ylabA51=['$\mathcal{A}_{42}/\rho AB$'];
ylabA53=['$\mathcal{A}_{43}/\rho AB$'];
ylabA55=['$\mathcal{A}_{44}/\rho AB^2$'];
ylabB11=['$\mathcal{B}_{22}/\omega\rho A$'];
ylabB13=['$\mathcal{B}_{23}/\omega\rho A$'];
ylabB15=['$\mathcal{B}_{24}/\omega\rho AB$'];
ylabB31=['$\mathcal{B}_{32}/\omega\rho A$'];
ylabB33=['$\mathcal{B}_{33}/\omega\rho A$'];
ylabB35=['$\mathcal{B}_{34}/\omega\rho AB$'];
ylabB51=['$\mathcal{B}_{42}/\omega\rho AB$'];
ylabB53=['$\mathcal{B}_{43}/\omega\rho AB$'];
ylabB55=['$\mathcal{B}_{44}/\omega\rho AB^2$'];
else
ylabA11=['$\mathcal{A}_{22}/\rho A$'];
ylabA13=['$\mathcal{A}_{23}/\rho A$'];
ylabA15=['$\mathcal{A}_{24}/\rho AB$'];
ylabA31=['$\mathcal{A}_{32}/\rho A$'];
ylabA33=['$\mathcal{A}_{33}/\rho A $'];
ylabA35=['$\mathcal{A}_{34}/\rho AB$'];
ylabA51=['$\mathcal{A}_{42}/\rho AB$'];
ylabA53=['$\mathcal{A}_{43}/\rho AB$'];
ylabA55=['$\mathcal{A}_{44}/\rho AB^2$'];
ylabB11=['$(\mathcal{B}_{22}/\rho A)\sqrt{B/2g}$'];
ylabB13=['$(\mathcal{B}_{23}/\rho A)\sqrt{B/2g} $'];
ylabB15=['$(\mathcal{B}_{23}/\rho AB)\sqrt{B/2g}$'];
ylabB31=['$(\mathcal{B}_{32}/\rho A)\sqrt{B/2g}$'];
ylabB33=['$(\mathcal{B}_{33}/\rho A)\sqrt{B/2g} $'];
ylabB35=['$(\mathcal{B}_{34}/\rho AB)\sqrt{B/2g}$'];
ylabB51=['$(\mathcal{B}_{42}/\rho AB)\sqrt{B/2g}$'];
ylabB53=['$(\mathcal{B}_{43}/\rho AB)\sqrt{B/2g} $'];
ylabB55=['$(\mathcal{B}_{44}/\rho AB^2)\sqrt{B/2g} $'];
end
ylabeta1=['$|\eta/\xi_2|$'];
ylabeta3=['$|\eta/\xi_3|$'];
ylabeta5=['$2|\eta/\xi_4|/B$'];

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
end

figure('name','Added mass');
subplot(3,3,1)
if flagVdat==0
    plot(wwnorm,A11,'r')
else
    if datval(1,1)~=0
        datA11=datval(:,1:2);
        datA11(datA11==0)=NaN;
        plot(datA11(:,1),datA11(:,2),'ob',wwnorm,A11,'r');
    else
        plot(wwnorm,A11,'r')
    end
end
ylabel(ylabA11,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;
       
subplot(3,3,2)
if flagVdat==0
    plot(wwnorm,A13,'r')
else
    if datval(1,3)~=0
        datA13=datval(:,3:4);
        datA13(datA13==0)=NaN;
        plot(datA13(:,1),datA13(:,2),'ob',wwnorm,A13,'r');
    else
        plot(wwnorm,A13,'r')
    end
end

ylabel(ylabA13,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,3)
if flagVdat==0
    plot(wwnorm,A15,'r')
else
    if datval(1,5)~=0
        datA15=datval(:,5:6);
        datA15(datA15==0)=NaN;
        plot(datA15(:,1),datA15(:,2),'ob',wwnorm,A15,'r');
    else
        plot(wwnorm,A15,'r')
    end
end

ylabel(ylabA15,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,4)
if flagVdat==0
    plot(wwnorm,A31,'r')
else
    if datval(1,7)~=0
        datA31=datval(:,7:8);
        datA31(datA31==0)=NaN;
        plot(datA31(:,1),datA31(:,2),'ob',wwnorm,A31,'r');
    else
        plot(wwnorm,A31,'r')
    end
end

ylabel(ylabA31,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,5)
if flagVdat==0
    plot(wwnorm,A33,'r')
else
    if datval(1,9)~=0
        datA33=datval(:,9:10);
        datA33(datA33==0)=NaN;
        plot(datA33(:,1),datA33(:,2),'ob',wwnorm,A33,'r');
    else
        plot(wwnorm,A33,'r')
    end
end

ylabel(ylabA33,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,6)
if flagVdat==0
    plot(wwnorm,A35,'r')
else
    if datval(1,11)~=0
        datA35=datval(:,11:12);
        datA35(datA35==0)=NaN;
        plot(datA35(:,1),datA35(:,2),'ob',wwnorm,A35,'r');
    else
        plot(wwnorm,A35,'r')
    end
end
ylabel(ylabA35,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;


subplot(3,3,7)
if flagVdat==0
    plot(wwnorm,A51,'r')
else
    if datval(1,13)~=0
        datA51=datval(:,13:14);
        datA51(datA51==0)=NaN;
        plot(datA51(:,1),datA51(:,2),'ob',wwnorm,A51,'r');
    else
        plot(wwnorm,A51,'r')
    end
end
ylabel(ylabA51,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,8)
if flagVdat==0
    plot(wwnorm,A53,'r')
else
    if datval(1,15)~=0
        datA53=datval(:,15:16);
        datA53(datA53==0)=NaN;
        plot(datA53(:,1),datA53(:,2),'ob',wwnorm,A53,'r');
    else
        plot(wwnorm,A53,'r')
    end
end

ylabel(ylabA53,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,9)
if flagVdat==0
    plot(wwnorm,A55,'r')
else
    if datval(1,17)~=0
        datA55=datval(:,17:18);
        datA55(datA55==0)=NaN;
        plot(datA55(:,1),datA55(:,2),'ob',wwnorm,A55,'r');
    else
        plot(wwnorm,A55,'r')
    end
end

ylabel(ylabA55,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

figure('name','Damping coefficients');
subplot(3,3,1)
if flagVdat==0
    plot(wwnorm,B11,'r')
else
    if datval(1,19)~=0
        datB11=datval(:,19:20);
        datB11(datB11==0)=NaN;
        plot(datB11(:,1),datB11(:,2),'ob',wwnorm,B11,'r');
    else
        plot(wwnorm,B11,'r')
    end
end

ylabel(ylabB11,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;
       
subplot(3,3,2)
if flagVdat==0
    plot(wwnorm,B13,'r')
else
    if datval(1,21)~=0
        datB13=datval(:,21:22);
        datB13(datB13==0)=NaN;
        plot(datB13(:,1),datB13(:,2),'ob',wwnorm,B13,'r');
    else
        plot(wwnorm,B13,'r')
    end
end
ylabel(ylabB13,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,3)
if flagVdat==0
    plot(wwnorm,B15,'r')
else
    if datval(1,23)~=0
        datB15=datval(:,23:24);
        datB15(datB15==0)=NaN;
        plot(datB15(:,1),datB15(:,2),'ob',wwnorm,B15,'r');
    else
        plot(wwnorm,B15,'r')
    end
end
ylabel(ylabB15,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,4)
if flagVdat==0
    plot(wwnorm,B31,'r')
else
    if datval(1,25)~=0
        datB31=datval(:,25:26);
        datB31(datB31==0)=NaN;
        plot(datB31(:,1),datB31(:,2),'ob',wwnorm,B31,'r');
    else
        plot(wwnorm,B31,'r')
    end
end
ylabel(ylabB31,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,5)
if flagVdat==0
    plot(wwnorm,B33,'r')
else
    if datval(1,27)~=0
        datB33=datval(:,27:28);
        datB33(datB33==0)=NaN;
        plot(datB33(:,1),datB33(:,2),'ob',wwnorm,B33,'r');
    else
        plot(wwnorm,B33,'r')
    end
end
ylabel(ylabB33,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,6)
if flagVdat==0
    plot(wwnorm,B35,'r')
else
    if datval(1,27)~=0
        datB35=datval(:,29:30);
        datB35(datB35==0)=NaN;
        plot(datB35(:,1),datB35(:,2),'ob',wwnorm,B35,'r');
    else
        plot(wwnorm,B35,'r')
    end
end
ylabel(ylabB35,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;


subplot(3,3,7)
if flagVdat==0
    plot(wwnorm,B51,'r')
else
    if datval(1,31)~=0
        datB51=datval(:,31:32);
        datB51(datB51==0)=NaN;
        plot(datB51(:,1),datB51(:,2),'ob',wwnorm,B51,'r');
    else
        plot(wwnorm,B51,'r')
    end
end

ylabel(ylabB51,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,8)
if flagVdat==0
    plot(wwnorm,B53,'r')
else
    if datval(1,33)~=0
        datB53=datval(:,33:34);
        datB53(datB53==0)=NaN;
        plot(datB53(:,1),datB53(:,2),'ob',wwnorm,B53,'r');
    else
        plot(wwnorm,B53,'r')
    end
end
ylabel(ylabB53,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

subplot(3,3,9)
if flagVdat==0
    plot(wwnorm,B55,'r')
else
    if datval(1,35)~=0
        datB55=datval(:,35:36);
        datB55(datB55==0)=NaN;
        plot(datB55(:,1),datB55(:,2),'ob',wwnorm,B55,'r');
    else
        plot(wwnorm,B55,'r')
    end
end
ylabel(ylabB55,'Interpreter','latex')
xlabel(xlabstr,'Interpreter','latex');
plot_properties;

figure('name','Radiated wave amplitude');
subplot(3,1,1)
if flagVdat==0
plot(wwnorm,eta1(:,1),'r')
else
     if datval(1,37)~=0
      dateta1=datval(:,37:38);
      dateta1(dateta1==0)=NaN;
      plot(dateta1(:,1),dateta1(:,2),'ob',wwnorm,eta1(:,1),'r')
     else
      plot(wwnorm,eta1(:,1),'r')
     end
end
xlabel(xlabstr,'Interpreter','latex');
ylabel(ylabeta1,'Interpreter','latex')
plot_properties;

subplot(3,1,2)
if flagVdat==0
plot(wwnorm,eta3(:,1),'r')
else
     if datval(1,39)~=0
      dateta3=datval(:,39:40);
      dateta3(dateta3==0)=NaN;
      plot(dateta3(:,1),dateta3(:,2),'ob',wwnorm,eta3(:,1),'r')
     else
      plot(wwnorm,eta3(:,1),'r')
     end
end

xlabel(xlabstr,'Interpreter','latex');
ylabel(ylabeta3,'Interpreter','latex')
plot_properties;

subplot(3,1,3)
if flagVdat==0
plot(wwnorm,eta5(:,1),'r')
else
     if datval(1,41)~=0
      dateta5=datval(:,41:42);
      dateta5(dateta5==0)=NaN;
      plot(dateta5(:,1),dateta5(:,2),'ob',wwnorm,eta5(:,1),'r')
     else
      plot(wwnorm,eta5(:,1),'r')
     end
end

xlabel(xlabstr,'Interpreter','latex');
ylabel(ylabeta5,'Interpreter','latex')
plot_properties;
