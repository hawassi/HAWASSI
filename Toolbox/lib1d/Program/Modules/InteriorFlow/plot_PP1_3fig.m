function plot_PP1_3fig(sf_savename,CheckId,X,Z,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,a_t,ax_t,az_t,StringVar,stringlabel,stringtitle,SetAx,Id,statusbarObj)
Check_Vz=CheckId.Vz;
Check_Vx=CheckId.Vx;
Check_V=CheckId.V;
Check_a=CheckId.a;
Check_ax=CheckId.ax;
Check_az=CheckId.az;
Check_P=CheckId.P;
Check_Pdyn=CheckId.Pdyn;
Check_Pnonlindyn=CheckId.Pnonlindyn;

FigP=[];FigV=[];FigA=[];


if  Check_P==1 && Check_Pdyn==0 &&Check_Pnonlindyn==0
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0 && Check_Pdyn==1 &&Check_Pnonlindyn==0
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0 && Check_Pdyn==0 &&Check_Pnonlindyn==1
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_P==1 && Check_Pdyn==1 &&Check_Pnonlindyn==0
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
elseif  Check_P==1 && Check_Pdyn==0 &&Check_Pnonlindyn==1
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
elseif  Check_P==0 && Check_Pdyn==1 &&Check_Pnonlindyn==1
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
elseif  Check_P==1 && Check_Pdyn==1 &&Check_Pnonlindyn==1
    FigP=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.05, 0.6,0.7]);
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{3},stringlabel,SetAx,Id)
end
set(FigP,'Renderer','zbuffer'); %due to graphics driver

if  Check_V==1 && Check_Vx==0 &&Check_Vz==0
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_V==0 && Check_Vx==1 &&Check_Vz==0
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_V==0 && Check_Vx==0 &&Check_Vz==1
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_V==1 && Check_Vx==1 &&Check_Vz==0
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
elseif  Check_V==1 && Check_Vx==0 &&Check_Vz==1
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
elseif  Check_V==0 && Check_Vx==1 &&Check_Vz==1
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
elseif  Check_V==1 && Check_Vx==1 &&Check_Vz==1
    FigV=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.05, 0.6,0.7]);
    subplot(3,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{6},stringlabel,SetAx,Id)
end
set(FigV,'Renderer','zbuffer'); %due to graphics driver

if  Check_a==1 && Check_ax==0 &&Check_az==0
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_a==0 && Check_ax==1 &&Check_az==0
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_a==0 && Check_ax==0 &&Check_az==1
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
    title(stringtitle);
elseif  Check_a==1 && Check_ax==1 &&Check_az==0
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
elseif  Check_a==1 && Check_ax==0 &&Check_az==1
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
elseif  Check_a==0 && Check_ax==1 &&Check_az==1
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    subplot(2,1,1)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
elseif  Check_a==1 && Check_ax==1 &&Check_az==1
    FigA=figure;
    set(gcf,'visible','off');
    set(gcf, 'units','normalized','Position', [0.2, 0.05, 0.6,0.7]);
    subplot(3,1,1)
    plot_PP1_Var(X,Z,a_t,StringVar{7},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,ax_t,StringVar{8},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,az_t,StringVar{9},stringlabel,SetAx,Id)
end
set(FigA,'Renderer','zbuffer'); %due to graphics driver

set(gcf,'Renderer','zbuffer'); %due to graphics driver


if CheckId.SaveFig==1
    statusbarObj.setText('Saving figures...');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.blue);
    if Id==1
        if ~isempty(FigP)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigP,[sf_savename,'PP1_Interior_plot_XZ_Pressure',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigP,[sf_savename,'PP1_Interior_plot_XZ_Pressure',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigV)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigV,[sf_savename,'PP1_Interior_plot_XZ_Velocity',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigV,[sf_savename,'PP1_Interior_plot_XZ_Velocity',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigA)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigA,[sf_savename,'PP1_Interior_plot_XZ_Acceleration',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigA,[sf_savename,'PP1_Interior_plot_XZ_Acceleration',CheckId.SaveFig_type])
            end
        end
    elseif Id==2
        if ~isempty(FigP)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigP,[sf_savename,'PP1_Interior_plot_ZT_Pressure',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigP,[sf_savename,'PP1_Interior_plot_ZT_Pressure',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigV)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigV,[sf_savename,'PP1_Interior_plot_ZT_Velocity',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigV,[sf_savename,'PP1_Interior_plot_ZT_Velocity',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigA)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigA,[sf_savename,'PP1_Interior_plot_ZT_Acceleration',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigA,[sf_savename,'PP1_Interior_plot_ZT_Acceleration',CheckId.SaveFig_type])
            end
        end
    elseif Id==3
        if ~isempty(FigP)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigP,[sf_savename,'PP1_Interior_plot_XT_Pressure',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigP,[sf_savename,'PP1_Interior_plot_XT_Pressure',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigV)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigV,[sf_savename,'PP1_Interior_plot_XT_Velocity',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigV,[sf_savename,'PP1_Interior_plot_XT_Velocity',CheckId.SaveFig_type])
            end
        end
        if ~isempty(FigA)
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigA,[sf_savename,'PP1_Interior_plot_XT_Acceleration',CheckId.SaveFig_type],'epsc')
            else
                saveas(FigA,[sf_savename,'PP1_Interior_plot_XT_Acceleration',CheckId.SaveFig_type])
            end
        end
    end
    statusbarObj.setText('');
end


if ~isempty(FigP) && isempty(FigV)  && isempty(FigA)
    set(FigP,'visible','on')
elseif isempty(FigP) && ~isempty(FigV)  && isempty(FigA)
    set(FigV,'visible','on')
elseif isempty(FigP) && isempty(FigV)  && ~isempty(FigA)
    set(FigA,'visible','on')
elseif ~isempty(FigP) && ~isempty(FigV)  && isempty(FigA)
    FigXZ=[FigP FigV]; tabname={'Pressure', 'Velocity'};
    figs2tabs([FigXZ(:)],tabname)
elseif ~isempty(FigP) && isempty(FigV)  && ~isempty(FigA)
    FigXZ=[FigP FigA]; tabname={'Pressure', 'Acceleration'};
    figs2tabs([FigXZ(:)],tabname)
elseif isempty(FigP) && ~isempty(FigV)  && ~isempty(FigA)
    FigXZ=[FigV FigA]; tabname={'Velocity','Acceleration'};
    figs2tabs([FigXZ(:)],tabname)
elseif ~isempty(FigP) && ~isempty(FigV)  && ~isempty(FigA)
    FigXZ=[FigP FigV FigA]; tabname={'Pressure','Velocity','Acceleration'};
    figs2tabs([FigXZ(:)],tabname)
end





