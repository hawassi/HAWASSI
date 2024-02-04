function plot_PP2_3fig(CheckId,X,Z,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,a_t,ax_t,az_t,StringVar,stringlabel,stringtitle,SetAx,Id)
Check_Vz=CheckId.Vz;
Check_Vx=CheckId.Vx;
Check_V=CheckId.V;
Check_a=CheckId.a;
Check_ax=CheckId.ax;
Check_az=CheckId.az;
Check_P=CheckId.P;
Check_Pdyn=CheckId.Pdyn;
Check_Pnonlindyn=CheckId.Pnonlindyn;


set(gcf,'Renderer','zbuffer'); %due to graphics driver
if  Check_P==1 && Check_Pdyn==0 &&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0 && Check_Pdyn==1 &&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0 && Check_Pdyn==0 &&Check_Pnonlindyn==1
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_P==1 && Check_Pdyn==1 &&Check_Pnonlindyn==0
     subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
elseif  Check_P==1 && Check_Pdyn==0 &&Check_Pnonlindyn==1
     subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
elseif  Check_P==0 && Check_Pdyn==1 &&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
elseif  Check_P==1 && Check_Pdyn==1 &&Check_Pnonlindyn==1
     subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
end


if  Check_V==1 && Check_Vx==0 &&Check_Vz==0
      plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_V==0 && Check_Vx==1 &&Check_Vz==0
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_V==0 && Check_Vx==0 &&Check_Vz==1
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_V==1 && Check_Vx==1 &&Check_Vz==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
elseif  Check_V==1 && Check_Vx==0 &&Check_Vz==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
elseif  Check_V==0 && Check_Vx==1 &&Check_Vz==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
elseif  Check_V==1 && Check_Vx==1 &&Check_Vz==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
end

if  Check_a==1 && Check_ax==0 &&Check_az==0
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_a==0 && Check_ax==1 &&Check_az==0
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_a==0 && Check_ax==0 &&Check_az==1
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_a==1 && Check_ax==1 &&Check_az==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
elseif  Check_a==1 && Check_ax==0 &&Check_az==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
elseif  Check_a==0 && Check_ax==1 &&Check_az==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
elseif  Check_a==1 && Check_ax==1 &&Check_az==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
end


 
