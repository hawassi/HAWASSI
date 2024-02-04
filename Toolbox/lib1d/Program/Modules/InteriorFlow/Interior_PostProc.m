% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

function Interior_PostProc(IP,PP,handles,statusbarObj)

statusbarObj.setText('Processing...');

CurDir=cellstr(get(handles.popupmenu_CurDir, 'string'));
projdir=CurDir{get(handles.popupmenu_CurDir,'Value')};
savename=IP.savename;
global OS
if strcmp(OS,'Windows')
    sf_savename = [projdir,'\',savename,'\'];
else
    sf_savename = [projdir,'/',savename,'/'];
end

if ~isdir(sf_savename)
    mkdir(sf_savename);
end

StringVar={'Pressure (Tot.) [N/m^2]', 'Pressure (dyn.) [N/m^2]', 'Pressure (Non-lin. dyn.) [N/m^2]',...
    'Velocity [m/s]', 'Horizontal Velocity [m/s]','Vertical Velocity [m/s]', ...
    'Acceleration [m/s^2]','Horizontal Acceleration [m/s^2]', 'Vertical Acceleration [m/s^2]'};


CheckId.P=get(handles.checkbox_P,'value');
CheckId.Pdyn=get(handles.checkbox_P_dyn,'value');
CheckId.Pnonlindyn=get(handles.checkbox_P_NonlinDyn,'value');
CheckId.V=get(handles.checkbox_V,'value');
CheckId.Vx=get(handles.checkbox_Vx,'value');
CheckId.Vz=get(handles.checkbox_Vz,'value');
CheckId.a=get(handles.checkbox_A,'value');
CheckId.ax=get(handles.checkbox_ax,'value');
CheckId.az=get(handles.checkbox_az,'value');
CheckId.SaveFig=get(handles.checkbox_savefigure,'value');
filetype=cellstr(get(handles.popupmenu_saveFig, 'string'));
CheckId.SaveFig_type=filetype{get(handles.popupmenu_saveFig,'value')};
CheckId.Save=get(handles.checkbox_savedata_PP,'value');
IdSavedata=get(handles.edit_saveId_PP,'string');

SetAx.Clim.Id=get(handles.checkbox22,'value');
SetAx.Clim.val=get(handles.edit_colorbar,'userdata');
SetAx.Xlim.Id=get(handles.checkbox25,'value');
SetAx.Xlim.val=get(handles.edit_xlim,'userdata');
SetAx.Ylim.Id=get(handles.checkbox_setting_ylim,'value');
SetAx.Ylim.val=get(handles.edit34,'userdata');
SetAx.tcoarse.Id=get(handles.checkbox_tcoarse,'value');
SetAx.tcoarse.val=get(handles.edit_PP_tcoarse,'userdata');
SetAx.GIF.Id=get(handles.checkbox_PP2_GIF,'value');
SetAx.GIF.val=get(handles.edit_PP2_GIF,'userdata');

if CheckId.P+CheckId.V+CheckId.Vx+CheckId.Vz+CheckId.Pdyn+CheckId.Pnonlindyn+...
        CheckId.a+CheckId.ax+CheckId.az==0
    statusbarObj.setText('Choose an option (Pressure or Velocities or Acceleration)!');
    statusbarTxt = statusbarObj.getComponent(0);
    statusbarTxt.setForeground(java.awt.Color.red);
    return
end;

P_tot=IP.Ptot;
Vx=IP.dxPHI;
Vz=IP.dzPHI;
V=sqrt(Vx.^2+Vz.^2);
ax=IP.dt_dxPHI;
az=IP.dt_dzPHI;
a=sqrt(ax.^2+az.^2);
P_dyn=-IP.dtPHI;
P_nonlindyn=-(IP.dtPHI+0.5.*(Vx.^2+Vz.^2));
X=IP.X;Z=IP.Z;
T=IP.timeIP;


if PP==1
    Check_XZ=get(handles.checkbox_XZ,'value');
    Check_ZT=get(handles.checkbox_ZT,'value');
    Check_XT=get(handles.checkbox_XT,'value');
    Check_X=get(handles.checkbox_PP1_line_X,'value');
    Check_Z=get(handles.checkbox_PP1_line_Z,'value');
    Check_T=get(handles.checkbox_PP1_line_T,'value');
    
    if Check_XZ+Check_ZT+Check_XT+Check_X+Check_Z+Check_T==0
        statusbarObj.setText('Choose an axes for plotting!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    end
    
    
    XZ_T   =get(handles.edit_PP1_XZ_time,'Userdata');
    ZT_X   =get(handles.edit_PP1_ZT_X,'Userdata');
    XT_Z   =get(handles.edit_PP1_XT_Z,'Userdata');
    
    line_X_TZ=get(handles.edit_PP1_line_X,'Userdata');
    line_Z_TX=get(handles.edit_PP1_line_Z,'Userdata');
    line_T_XZ=get(handles.edit_PP1_line_T,'Userdata');
    
    
    if Check_XZ==1
        if isempty(XZ_T)||length(XZ_T)~=1
            statusbarObj.setText('Input a time!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if XZ_T<T(1)
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(T(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif XZ_T>T(end)
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(T(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        indt=closest(T,XZ_T);
        P_t=squeeze(P_tot(indt,:,:))';
        Pdyn_t=squeeze(P_dyn(indt,:,:))';
        Pnonlindyn_t=squeeze(P_nonlindyn(indt,:,:))';
        V_t=squeeze(V(indt,:,:))';
        Vx_t=squeeze(Vx(indt,:,:))';
        Vz_t=squeeze(Vz(indt,:,:))';
        a_t=squeeze(a(indt,:,:))';
        ax_t=squeeze(ax(indt,:,:))';
        az_t=squeeze(az(indt,:,:))';
        
        if CheckId.Save==1
            data_interior_PP.X=X;
            data_interior_PP.Z=Z;
            data_interior_PP.T=T(indt);
            data_interior_PP.P=P_t;
            data_interior_PP.Pdyn=Pdyn_t;
            data_interior_PP.Pnonlindyn=Pnonlindyn_t;
            data_interior_PP.V=V_t;
            data_interior_PP.Vx=Vx_t;
            data_interior_PP.Vz=Vz_t;
            data_interior_PP.a=a_t;
            data_interior_PP.ax=ax_t;
            data_interior_PP.az=az_t;
            save ('-v7.3',[sf_savename,'data_interior_PP_T','_',IdSavedata,'.mat'],'data_interior_PP');
        end
    end
    
    if Check_ZT==1
        
        
        if isempty(ZT_X)||length(ZT_X)~=1
            statusbarObj.setText('Input a horizontal point (x)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if ZT_X<X(1)
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(X(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif ZT_X>X(end)
            statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(X(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        indx=closest(X,ZT_X);
        P_x=squeeze(P_tot(:,indx,:))';
        Pdyn_x=squeeze(P_dyn(:,indx,:))';
        Pnonlindyn_x=squeeze(P_nonlindyn(:,indx,:))';
        V_x=squeeze(V(:,indx,:))';
        Vx_x=squeeze(Vx(:,indx,:))';
        Vz_x=squeeze(Vz(:,indx,:))';
        a_x=squeeze(a(:,indx,:))';
        ax_x=squeeze(ax(:,indx,:))';
        az_x=squeeze(az(:,indx,:))';
        
        if CheckId.Save==1
            data_interior_PP.X=X(indx);
            data_interior_PP.Z=Z;
            data_interior_PP.T=T;
            data_interior_PP.P=P_x;
            data_interior_PP.Pdyn=Pdyn_x;
            data_interior_PP.Pnonlindyn=Pnonlindyn_x;
            data_interior_PP.V=V_x;
            data_interior_PP.Vx=Vx_x;
            data_interior_PP.Vz=Vz_x;
            data_interior_PP.a=a_x;
            data_interior_PP.ax=ax_x;
            data_interior_PP.az=az_x;
            save ('-v7.3',[sf_savename,'data_interior_PP_X','_',IdSavedata,'.mat'],'data_interior_PP');
        end
        
    end
    
    if Check_XT==1
        
        
        if isempty(XT_Z)||length(XT_Z)~=1
            statusbarObj.setText('Input a vertical point (z)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if XT_Z<Z(1)
            statusbarObj.setText(['Error:Input value z< minimum value (z=',num2str(Z(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif XT_Z>Z(end)
            statusbarObj.setText(['Error:Input value z> maximum value (z=',num2str(Z(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        indz=closest(Z,XT_Z);
        P_z=squeeze(P_tot(:,:,indz));
        Pdyn_z=squeeze(P_dyn(:,:,indz));
        Pnonlindyn_z=squeeze(P_nonlindyn(:,:,indz));
        V_z=squeeze(V(:,:,indz));
        Vx_z=squeeze(Vx(:,:,indz));
        Vz_z=squeeze(Vz(:,:,indz));
        a_z=squeeze(a(:,:,indz));
        ax_z=squeeze(ax(:,:,indz));
        az_z=squeeze(az(:,:,indz));
        
        if CheckId.Save==1
            data_interior_PP.X=X;
            data_interior_PP.Z=Z(indz);
            data_interior_PP.T=T;
            data_interior_PP.P=P_z;
            data_interior_PP.Pdyn=Pdyn_z;
            data_interior_PP.Pnonlindyn=Pnonlindyn_z;
            data_interior_PP.V=V_z;
            data_interior_PP.Vx=Vx_z;
            data_interior_PP.Vz=Vz_z;
            data_interior_PP.a=a_z;
            data_interior_PP.ax=ax_z;
            data_interior_PP.az=az_z;
            save ('-v7.3',[sf_savename,'data_interior_PP_Z','_',IdSavedata,'.mat'],'data_interior_PP');
        end
    end
    
    if SetAx.Clim.Id==1
        if isempty(SetAx.Clim.val)|| length(SetAx.Clim.val)~=2
            statusbarObj.setText('Input colorbar limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if SetAx.Clim.val(1)>SetAx.Clim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [min;max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
    end
    
    if SetAx.Xlim.Id==1
        if isempty(SetAx.Xlim.val)|| length(SetAx.Xlim.val)~=2
            statusbarObj.setText('Input horizontal axes limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if SetAx.Xlim.val(1)>SetAx.Xlim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [x_min;x_max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
    end
    
    
    
    if SetAx.Ylim.Id==1
        if isempty(SetAx.Ylim.val)|| length(SetAx.Ylim.val)~=2
            statusbarObj.setText('Input vertical axes limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if SetAx.Ylim.val(1)>SetAx.Ylim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [y_min;y_max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
    end
    
    if Check_X==1
        
        if isempty(line_X_TZ)||length(line_X_TZ)~=2
            statusbarObj.setText('Input a time and a point at vertical axes (z)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if line_X_TZ(1)<T(1)
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(T(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_X_TZ(1)>T(end)
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(T(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if line_X_TZ(2)<Z(1)
            statusbarObj.setText(['Error:Input value z< minimum value (z=',num2str(Z(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_X_TZ(2)>Z(end)
            statusbarObj.setText(['Error:Input value z> maximum value (z=',num2str(Z(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        
        indt=closest(T,line_X_TZ(1));
        indz=closest(Z,line_X_TZ(2));
        
        P_tz=squeeze(P_tot(indt,:,indz))';
        Pdyn_tz=squeeze(P_dyn(indt,:,indz))';
        Pnonlindyn_tz=squeeze(P_nonlindyn(indt,:,indz))';
        V_tz=squeeze(V(indt,:,indz))';
        Vx_tz=squeeze(Vx(indt,:,indz))';
        Vz_tz=squeeze(Vz(indt,:,indz))';
        a_tz=squeeze(a(indt,:,indz))';
        ax_tz=squeeze(ax(indt,:,indz))';
        az_tz=squeeze(az(indt,:,indz))';
        
        if CheckId.Save==1
            data_interior_PP.X=X;
            data_interior_PP.Z=Z(indz);
            data_interior_PP.T=T(indt);
            data_interior_PP.P=P_tz;
            data_interior_PP.Pdyn=Pdyn_tz;
            data_interior_PP.Pnonlindyn=Pnonlindyn_tz;
            data_interior_PP.V=V_tz;
            data_interior_PP.Vx=Vx_tz;
            data_interior_PP.Vz=Vz_tz;
            data_interior_PP.a=a_tz;
            data_interior_PP.ax=ax_tz;
            data_interior_PP.az=az_tz;
            save ('-v7.3',[sf_savename,'data_interior_PP_','TZ','_',IdSavedata,'.mat'],'data_interior_PP');
        end
    end
    
    if Check_Z==1
        
        if isempty(line_Z_TX)||length(line_Z_TX)~=2
            statusbarObj.setText('Input a time and a point at horizonta axes (x)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        
        if line_Z_TX(1)<T(1)
            statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(T(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_Z_TX(1)>T(end)
            statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(T(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if line_Z_TX(2)<X(1)
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(X(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_Z_TX(2)>X(end)
            statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(X(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        
        indt=closest(T,line_Z_TX(1));
        indx=closest(X,line_Z_TX(2));
        
        P_tx=squeeze(P_tot(indt,indx,:))';
        Pdyn_tx=squeeze(P_dyn(indt,indx,:))';
        Pnonlindyn_tx=squeeze(P_nonlindyn(indt,indx,:))';
        V_tx=squeeze(V(indt,indx,:))';
        Vx_tx=squeeze(Vx(indt,indx,:))';
        Vz_tx=squeeze(Vz(indt,indx,:))';
        a_tx=squeeze(a(indt,indx,:))';
        ax_tx=squeeze(ax(indt,indx,:))';
        az_tx=squeeze(az(indt,indx,:))';
        
        if CheckId.Save==1
            data_interior_PP.X=X(indx);
            data_interior_PP.Z=Z;
            data_interior_PP.T=T(indt);
            data_interior_PP.P=P_tx;
            data_interior_PP.Pdyn=Pdyn_tx;
            data_interior_PP.Pnonlindyn=Pnonlindyn_tx;
            data_interior_PP.V=V_tx;
            data_interior_PP.Vx=Vx_tx;
            data_interior_PP.Vz=Vz_tx;
            data_interior_PP.a=a_tx;
            data_interior_PP.ax=ax_tx;
            data_interior_PP.az=az_tx;
            save ('-v7.3',[sf_savename,'data_interior_PP_','TX','_',IdSavedata,'.mat'],'data_interior_PP');
        end
    end
    
    if Check_T==1
        
        if isempty(line_T_XZ)||length(line_T_XZ)~=2
            statusbarObj.setText('Input a point at horizontal axes (x) and at vertical axes (z)!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        
        if line_T_XZ(1)<X(1)
            statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(X(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_T_XZ(1)>X(end)
            statusbarObj.setText(['Error:Input value x> maximum value (t=',num2str(X(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if line_T_XZ(2)<Z(1)
            statusbarObj.setText(['Error:Input value z< minimum value (z=',num2str(Z(1)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        elseif line_T_XZ(2)>Z(end)
            statusbarObj.setText(['Error:Input value z> maximum value (z=',num2str(Z(end)),')']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        
        indx=closest(X,line_T_XZ(1));
        indz=closest(Z,line_T_XZ(2));
        
        P_xz=squeeze(P_tot(:,indx,indz))';
        Pdyn_xz=squeeze(P_dyn(:,indx,indz))';
        Pnonlindyn_xz=squeeze(P_nonlindyn(:,indx,indz))';
        V_xz=squeeze(V(:,indx,indz))';
        Vx_xz=squeeze(Vx(:,indx,indz))';
        Vz_xz=squeeze(Vz(:,indx,indz))';
        a_xz=squeeze(a(:,indx,indz))';
        ax_xz=squeeze(ax(:,indx,indz))';
        az_xz=squeeze(az(:,indx,indz))';
        
        if CheckId.Save==1
            data_interior_PP.X=X(indx);
            data_interior_PP.Z=Z(indz);
            data_interior_PP.T=T;
            data_interior_PP.P=P_xz;
            data_interior_PP.Pdyn=Pdyn_xz;
            data_interior_PP.Pnonlindyn=Pnonlindyn_xz;
            data_interior_PP.V=V_xz;
            data_interior_PP.Vx=Vx_xz;
            data_interior_PP.Vz=Vz_xz;
            data_interior_PP.a=a_xz;
            data_interior_PP.ax=ax_xz;
            data_interior_PP.az=az_xz;
            save ('-v7.3',[sf_savename,'data_interior_PP_','XZ','_',IdSavedata,'.mat'],'data_interior_PP');
        end
    end
    
    
    
    
    
    if Check_XZ==1
        
        stringlabel={'x[m]','z[m]'};
        stringtitle=['Mesh plot @ time: ',num2str(roundn(XZ_T,-2)), ' [s]'];
        
        plot_PP1_3fig(sf_savename,CheckId,X,Z,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,a_t,ax_t,az_t,StringVar,stringlabel,stringtitle,SetAx,1,statusbarObj);
        
    end
    
    if Check_ZT==1
        
        stringlabel={'time[s]','z[m]'};
        stringtitle=['Mesh plot @ x: ',num2str(roundn(ZT_X,-2)), ' [m]'];
        plot_PP1_3fig(sf_savename,CheckId,T,Z,P_x,Pdyn_x,Pnonlindyn_x,V_x,Vx_x,Vz_x,a_x,ax_x,az_x,StringVar,stringlabel,stringtitle,SetAx,2,statusbarObj);
        
    end
    
    
    if Check_XT==1
        stringlabel={'x[m]','time[s]'};
        stringtitle=['Mesh plot @ z: ',num2str(roundn(XT_Z,-2)), ' [m]'];
        plot_PP1_3fig(sf_savename,CheckId,X,T,P_z,Pdyn_z,Pnonlindyn_z,V_z,Vx_z,Vz_z,a_z,ax_z,az_z,StringVar,stringlabel,stringtitle,SetAx,3,statusbarObj);
    end
    
    if Check_X==1
        stringlabel={'x[m]'};
        stringtitle=['plot @ time: ',num2str(roundn(line_X_TZ(1),-2)), ' [s]',...
            ' @ z: ',num2str(roundn(line_X_TZ(2),-2)), ' [m]'];
        plot_PP1_3fig_line(sf_savename,CheckId,X,P_tz,Pdyn_tz,Pnonlindyn_tz,V_tz,Vx_tz,Vz_tz,a_tz,ax_tz,az_tz,StringVar,stringlabel,stringtitle,SetAx,1,statusbarObj)
        
    end
    
    if Check_Z==1
        FigZ=figure;
        set(FigZ,'visible','off');
        stringlabel={'z[m]'};
        stringtitle=['plot @ time: ',num2str(roundn(line_Z_TX(1),-2)), ' [s]',...
            ' @ x: ',num2str(roundn(line_Z_TX(2),-2)), ' [m]'];
        plot_PP1_3fig_line(sf_savename,CheckId,Z,P_tx,Pdyn_tx,Pnonlindyn_tx,V_tx,Vx_tx,Vz_tx,a_tx,ax_tx,az_tx,StringVar,stringlabel,stringtitle,SetAx,2,statusbarObj)
    end
    
    
    if Check_T==1
        FigT=figure;
        set(FigT,'visible','off');
        stringlabel={'time [s]'};
        stringtitle=['plot @ x: ',num2str(roundn(line_T_XZ(1),-2)), ' [m]',...
            ' @ z: ',num2str(roundn(line_T_XZ(2),-2)), ' [m]'];
        plot_PP1_3fig_line(sf_savename,CheckId,T,P_xz,Pdyn_xz,Pnonlindyn_xz,V_xz,Vx_xz,Vz_xz,a_xz,ax_xz,az_xz,StringVar,stringlabel,stringtitle,SetAx,3,statusbarObj)
        
    end
    
    
    
    
elseif PP==2
    Tinterv =get(handles.edit_PP2_T,'Userdata');
    Xinterv =get(handles.edit_PP2_X,'Userdata');
    Zinterv =get(handles.edit_PP2_Z,'Userdata');
    
    
    if isempty(Tinterv)||length(Tinterv)~=2
        statusbarObj.setText('Input time interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Tinterv(1)<T(1)
        statusbarObj.setText(['Error:Input value t< minimum value (t=',num2str(T(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Tinterv(2)>T(end)
        statusbarObj.setText(['Error:Input value t> maximum value (t=',num2str(T(end)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Tinterv(1)>Tinterv(2)
        statusbarObj.setText(['Error:Wrong input format! [t_start;t_end]']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    
    if isempty(Xinterv)||length(Xinterv)~=2
        statusbarObj.setText('Input horizontal interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Xinterv(1)<X(1)
        statusbarObj.setText(['Error:Input value x< minimum value (x=',num2str(X(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Xinterv(2)>X(end)
        statusbarObj.setText(['Error:Input value x> maximum value (x=',num2str(X(end)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Xinterv(1)>Xinterv(2)
        statusbarObj.setText(['Error:Wrong input format! [x_start;x_end]']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if isempty(Zinterv)||length(Zinterv)~=2
        statusbarObj.setText('Input vertical interval!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if  Zinterv(1)<Z(1)
        statusbarObj.setText(['Error:Input value z< minimum value (z=',num2str(Z(1)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    elseif Zinterv(2)>Z(end)
        statusbarObj.setText(['Error:Input value z> maximum value (z=',num2str(Z(end)),')']);
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    
    if CheckId.P+CheckId.Pdyn+CheckId.Pnonlindyn>0 && CheckId.V+CheckId.Vx+CheckId.Vz>0 ...
            && CheckId.a+CheckId.ax+CheckId.az>0
        statusbarObj.setText('Choose either Pressure, velocity or acceleration !');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    elseif CheckId.P+CheckId.Pdyn+CheckId.Pnonlindyn>0 && CheckId.V+CheckId.Vx+CheckId.Vz>0 ...
            && CheckId.a+CheckId.ax+CheckId.az==0
        statusbarObj.setText('Choose either Pressure, Velocity or Acceleration !');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    elseif CheckId.P+CheckId.Pdyn+CheckId.Pnonlindyn>0 && CheckId.V+CheckId.Vx+CheckId.Vz==0 ...
            && CheckId.a+CheckId.ax+CheckId.az>0
        statusbarObj.setText('Choose either Pressure, Velocity or Acceleration !');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    elseif CheckId.P+CheckId.Pdyn+CheckId.Pnonlindyn==0 && CheckId.V+CheckId.Vx+CheckId.Vz>0 ...
            && CheckId.a+CheckId.ax+CheckId.az>0
        statusbarObj.setText('Choose either Pressure, Velocity or Acceleration !');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    else
        statusbarObj.setText('');
    end
    
    if SetAx.Clim.Id==1
        if isempty(SetAx.Clim.val)|| length(SetAx.Clim.val)~=2
            statusbarObj.setText('Input colorbar limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
        if SetAx.Clim.val(1)>SetAx.Clim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [min;max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
    end
    
    if SetAx.Xlim.Id==1
        if isempty(SetAx.Xlim.val)|| length(SetAx.Xlim.val)~=2
            statusbarObj.setText('Input horizontal axes limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if SetAx.Xlim.val(1)>SetAx.Xlim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [x_min;x_max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        
    end
    
    
    
    if SetAx.Ylim.Id==1
        if isempty(SetAx.Ylim.val)|| length(SetAx.Ylim.val)~=2
            statusbarObj.setText('Input vertical axes limit!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if SetAx.Ylim.val(1)>SetAx.Ylim.val(2)
            statusbarObj.setText(['Error:Wrong input format! [y_min;y_max]']);
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
    end
    
    if SetAx.tcoarse.Id==1
        if isempty(SetAx.tcoarse.val)|| length(SetAx.tcoarse.val)~=1
            statusbarObj.setText('Input tcoarse!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if fix(SetAx.tcoarse.val)~=SetAx.tcoarse.val
            statusbarObj.setText('Input only an integer value of tcoarse!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
        if SetAx.tcoarse.val<0
            statusbarObj.setText('Error: tcoarse is negative.');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return;
        end
        
        if SetAx.tcoarse.val==0
            statusbarObj.setText('Error: tcoarse is zero.');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return;
        end
        
    end
    
    if SetAx.GIF.Id==1
        if isempty(SetAx.GIF.val)|| length(SetAx.GIF.val)~=2
            statusbarObj.setText('Input GIF setting!');
            statusbarTxt = statusbarObj.getComponent(0);
            statusbarTxt.setForeground(java.awt.Color.red);
            return
        else
            statusbarObj.setText('');
        end
    end
    
    
    
    
    GIF_delaytime=SetAx.GIF.val(1);
    GIF_loopcount=SetAx.GIF.val(2);
    
    indti=closest(T,Tinterv(1));
    indtf=closest(T,Tinterv(2));
    
    if SetAx.tcoarse.Id==1
        stepdt=SetAx.tcoarse.val;
    else
        stepdt=1;
    end
    
    xinit=Xinterv(1);xend=Xinterv(2);
    stepdx=1;
    indxI=closest(X,xinit);indxE=closest(X,xend);
    zinit=Zinterv(1);zend=Zinterv(2);
    stepdz=1;
    indzI=closest(Z,zinit);indzE=closest(Z,zend);
    
    
    SetAx.Xlim.Id=1;
    SetAx.Xlim.val=[Xinterv(1) Xinterv(2)];
    
    
    filename = [sf_savename,'PP2_Interior_animation.gif'];
    hf2=figure;
    numframes=length(T(indti:stepdt:indtf));
    stringlabel={'x[m]','z[m]'};
    PP2_X=X(indxI:stepdx:indxE);
    PP2_Z=Z(indzI:stepdz:indzE);
    
    if CheckId.P+CheckId.Pdyn+CheckId.Pnonlindyn==3 ||CheckId.V+CheckId.Vx+CheckId.Vz==3 ...
            ||+CheckId.a+CheckId.ax+CheckId.az==3
        set(hf2, 'units','normalized','Position', [0.2, 0.05, 0.6,0.7]);
    else
        set(hf2, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
    end
    
    
    
    for i=1:numframes
        figure(hf2);
        indt=indti+(i-1)*stepdt;
        if indt>=indtf, break; end;
        
        P_t=squeeze(P_tot(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        Pdyn_t=squeeze(P_dyn(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        Pnonlindyn_t=squeeze(P_nonlindyn(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        V_t=squeeze(V(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        Vx_t=squeeze(Vx(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        Vz_t=squeeze(Vz(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        a_t=squeeze(a(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        ax_t=squeeze(ax(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        az_t=squeeze(az(indt,indxI:stepdx:indxE,indzI:stepdz:indzE))';
        
        stringtitle=['Mesh plot @ time: ',num2str(roundn(T(indt),-2)), ' [s]'];
        plot_PP2_3fig(CheckId,PP2_X,PP2_Z,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,a_t,ax_t,az_t,StringVar,stringlabel,stringtitle,SetAx,1);
        
        
        pause(GIF_delaytime)
        drawnow
        if ~ishandle(hf2)
            break;
        end
        
        try
            frame = getframe(hf2);
        catch
            
        end
        
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1;
            imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime,'loopcount',GIF_loopcount);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',GIF_delaytime, 'Writemode', 'append');
        end
    end
    
elseif PP==3
    MeasData=get(handles.pushbutton_measdata,'Userdata');
    Xmeas =MeasData(1,1);
    Zmeas =MeasData(2,2:end);
    Tmeas =MeasData(3:end,1);
    Depth= MeasData(1,2);
    PhysVal_m=MeasData(3:end,2:end);
    IdDensityPlot=get(handles.checkbox_PP3_densityplot,'value');
    IdErrorPlot=get(handles.checkbox_PP3_errorplot,'value');
    IdTimemeanPlot=get(handles.checkbox_timeMeanVariation,'value');
    IdLinePlot=get(handles.checkbox_PP3_lineplot,'value');
    if CheckId.P+CheckId.V+CheckId.Vx+CheckId.Vz+CheckId.Pdyn+CheckId.Pnonlindyn+CheckId.ax+CheckId.az~=1
        statusbarObj.setText('Choose only one physical quantity (pressure, velocity or acceleration)!');
        statusbarTxt = statusbarObj.getComponent(0);
        statusbarTxt.setForeground(java.awt.Color.red);
        return
    end;
    
    indx=closest(X,Xmeas);
    if CheckId.P==1
        PhysVal_s= squeeze(P_tot(:,indx,:));
        varname=StringVar{1};
        varnameS='P';
    elseif CheckId.Pdyn==1
        PhysVal_s= squeeze(P_dyn(:,indx,:));
        varname=StringVar{2};
        varnameS='Pdyn';
    elseif CheckId.Pnonlindyn==1
        PhysVal_s= squeeze(P_nonlindyn(:,indx,:));
        varname=StringVar{3};
        varnameS='Pnonlindyn';
    elseif CheckId.V==1
        PhysVal_s= squeeze(V(:,indx,:));
        varname=StringVar{4};
        varnameS='U';
    elseif CheckId.Vx==1
        PhysVal_s= squeeze(Vx(:,indx,:));
        varname=StringVar{5};
        varnameS='u';
    elseif CheckId.Vz==1
        PhysVal_s= squeeze(Vz(:,indx,:));
        varname=StringVar{6};
        varnameS='v';
    elseif CheckId.a==1
        PhysVal_s= squeeze(a(:,indx,:));
        varname=StringVar{7};
        varnameS='a';
    elseif CheckId.ax==1
        PhysVal_s= squeeze(ax(:,indx,:));
        varname=StringVar{8};
        varnameS='ax';
    elseif CheckId.az==1
        PhysVal_s= squeeze(az(:,indx,:));
        varname=StringVar{9};
        varnameS='az';
    end
    
    
    
    stringlabel={'time[s]','z[m]'};
    stringtitle1=['Plot @ x: ',num2str(roundn(Xmeas,-2)), ' [m]', ' (Simulation)'];
    stringtitle2=['Plot @ x: ',num2str(roundn(Xmeas,-2)), ' [m]', ' (Measurement)'];
    stringtitle3=['Plot @ x: ',num2str(roundn(Xmeas,-2)), ' [m]', ' (Error (%))'];
    
    
    if IdDensityPlot+IdErrorPlot>0 && IdTimemeanPlot==0 &&IdLinePlot==0
        if SetAx.Xlim.Id==1
            indt1=closest(T,SetAx.Xlim.val(1));
            indt2=closest(T,SetAx.Xlim.val(2));
            Tsimul=T(indt1:indt2);
        else
            Tsimul=T;
        end
        
        indZsimul=closest(Z,Zmeas);
        indTmeas=closest(Tmeas,Tsimul); %follows time of simulation
        PhysVal_meas=PhysVal_m(indTmeas,:);
        
        PhysVal_simul=PhysVal_s(:,indZsimul);
        PhysVal_meas=Adjust_time_simul_measPP3(PhysVal_simul,PhysVal_meas,Zmeas);
        
        PhysVal_simul=PhysVal_simul';
        PhysVal_meas=PhysVal_meas';
        
        
        Error=((PhysVal_simul-PhysVal_meas)...
            ./max(max(max(PhysVal_simul)),max(max(PhysVal_meas)))).*100;
        
        if IdDensityPlot==1 && IdErrorPlot==0
            figure
            set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
            subplot(2,1,1)
            plot_PP1_Var(Tsimul,Zmeas,PhysVal_simul,varname,stringlabel,SetAx,2)
            title(stringtitle1);
            subplot(2,1,2)
            plot_PP1_Var(Tsimul,Zmeas,PhysVal_meas,varname,stringlabel,SetAx,2)
            title(stringtitle2);
            set(gcf,'Renderer','zbuffer'); %due to graphics driver
            if CheckId.SaveFig==1
                if strcmp(CheckId.SaveFig_type,'.eps')
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_density_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type],'epsc')
                else
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_density_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type])
                end
            end
        elseif IdDensityPlot==1 &&IdErrorPlot==1
            figure
            set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
            subplot(3,1,1)
            plot_PP1_Var(Tsimul,Zmeas,PhysVal_simul,varname,stringlabel,SetAx,2)
            title(stringtitle1);
            subplot(3,1,2)
            plot_PP1_Var(Tsimul,Zmeas,PhysVal_meas,varname,stringlabel,SetAx,2)
            title(stringtitle2);
            subplot(3,1,3)
            plot_PP1_Var(Tsimul,Zmeas,Error,'Error(%)',stringlabel,SetAx,2)
            title(stringtitle3);
            set(gcf,'Renderer','zbuffer'); %due to graphics driver
            if CheckId.SaveFig==1
                if strcmp(CheckId.SaveFig_type,'.eps')
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_density_error_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type],'epsc')
                else
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_density_error_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type])
                end
            end
        elseif IdDensityPlot==0 &&IdErrorPlot==1
            figure
            set(gcf, 'units','normalized','Position', [0.2, 0.2, 0.6,0.5]);
            plot_PP1_Var(Tsimul,Zmeas,Error,'Error(%)',stringlabel,SetAx,2);
            title(stringtitle3);
            set(gcf,'Renderer','zbuffer'); %due to graphics driver
            if CheckId.SaveFig==1
                if strcmp(CheckId.SaveFig_type,'.eps')
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_error_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type],'epsc')
                else
                    saveas(gcf,[sf_savename,'PP3_',varnameS,'_error_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type])
                end
            end
        end
        
        
        
    elseif IdLinePlot==1 && IdTimemeanPlot==0 && IdDensityPlot+IdErrorPlot==0
        if SetAx.Xlim.Id==1
            indt1=closest(T,SetAx.Xlim.val(1));
            indt2=closest(T,SetAx.Xlim.val(2));
            Tsimul=T(indt1:indt2);
        else
            Tsimul=T;
        end
        
        indTmeas=closest(Tmeas,Tsimul); %follows time of simulation
        PhysVal_meas=PhysVal_m(indTmeas,:);
        PhysVal_simul=PhysVal_s;
        
        if get(handles.checkbox39,'value')==1
            dZstep=get(handles.edit_dZmeas,'Userdata');
            Zlineplot=Zmeas(1:dZstep:end);
        else
            Zlineplot=get(handles.edit_PP3_Z,'userdata');
        end
        
        Nfig=length(Zlineplot);
        FigPP3_line=figure;
        
        FlagYlabel=0;
        for i=1:Nfig
            indZsimul=closest(Z,Zlineplot(i));
            indZmeas =closest(Zmeas,Zlineplot(i));
            tempsimul=PhysVal_simul(:,indZsimul);
            tempsimul(isnan(tempsimul))=zeros;
            [Xcorel,lags]=xcorr(tempsimul,PhysVal_meas(:,indZmeas),'coeff');
            OptCorr1=max(Xcorel);
            indMaxCorel=closest(Xcorel,OptCorr1);
            OptCorr2=lags(indMaxCorel);
            PhysVal_measOpt = circshift(PhysVal_meas(:,indZmeas),OptCorr2);
            
            subaxis(Nfig,1,Nfig-i+1, 'Spacing', 0, 'Padding', 0, 'Margin', 0.13);
            plot(T,PhysVal_measOpt,'-b',T,PhysVal_simul(:,indZsimul),'--r');
            
            MM1 =max([max(PhysVal_measOpt) max(PhysVal_simul(:,indZsimul))]);
            MM =MM1+MM1/2;
            if i==1
                set(gca, 'box','off');
                xlabel('time [s]');
            else
                set(gca,'xtick',[],'xcolor','w');
                set(gca,'box','off');
            end
            
            xlim([Tsimul(1),Tsimul(end)]);
            ylim([-MM MM])
            if mod(i,2)==0
                set(gca,'ytick',[]);
            else
                if MM1*10>10
                    set(gca,'ytick',[-roundn(MM1,0) roundn(MM1,0)]);
                elseif MM1*10>1
                    set(gca,'ytick',[-roundn(MM1,-1) roundn(MM1,-1)]);
                else
                    set(gca,'ytick',[-roundn(MM1,-2) roundn(MM1,-2)]);
                end
            end
            th1=text(Tsimul(1)+1, MM/2,['Z=',num2str(roundn(Z(indZsimul),-2))]);
            
            
            if i==round(Nfig/2) &&  mod(i,2)~=0
                ylabel(varname);
                FlagYlabel=1;
            end
            if FlagYlabel==0 && i==round(Nfig/2)+1
                ylabel(varname);
            end
            
            plot_properties;
        end
        
        
        
        if Nfig>15
            set(FigPP3_line, 'units','normalized','Position', [0.1, 0.1, 0.3,0.8]);
        else
            set(FigPP3_line, 'units','normalized','Position', [0.1, 0.1, 0.6,0.5]);
        end
        if CheckId.SaveFig==1
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(FigPP3_line,[sf_savename,'PP3_',varnameS,'_line_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type],'epsc')
            else
                saveas(FigPP3_line,[sf_savename,'PP3_',varnameS,'_line_plot_at_x=',num2str(roundn(Xmeas,-2)),CheckId.SaveFig_type])
            end
        end
    elseif IdTimemeanPlot==1 && IdDensityPlot+IdErrorPlot==0 &&IdLinePlot==0
        Tsimul=T;
        indZsimul=closest(Z,Zmeas);
        indTmeas=closest(Tmeas,Tsimul); %follows time of simulation
        PhysVal_meas=PhysVal_m(indTmeas,:);
        PhysVal_simul=PhysVal_s(:,indZsimul);
        
        [TMPhysVal_simul,TMPhysVal_meas]=timemeanPP3(PhysVal_simul,PhysVal_meas,Zmeas,Depth,Tsimul);
        %         assignin('base','TMPhysVal_simul',TMPhysVal_simul)
        %         assignin('base','TMPhysVal_meas',TMPhysVal_meas)
        
        indx=closest(X,Xmeas);
        Zrel=(Zmeas-IP.meanEta(indx))./Depth;
        figure
        plot(TMPhysVal_meas,Zrel,'*b',TMPhysVal_simul,Zrel,'or')
        
        legend('Meas','Simul')
        if SetAx.Xlim.Id==1
            xlim([SetAx.Xlim.val(1),SetAx.Xlim.val(2)]);
        end
        
         if SetAx.Ylim.Id==1
            ylim([SetAx.Ylim.val(1),SetAx.Ylim.val(2)]);
        end
        
        xlabel(['$\bar{',varnameS,'}/\sqrt{gD}$'],'interpreter','latex')
        ylabel('$(z-\bar{\eta})/D$','interpreter','latex');
        plot_properties;
        if CheckId.SaveFig==1
            if strcmp(CheckId.SaveFig_type,'.eps')
                saveas(gcf,[sf_savename,'PP3_timemean_',varnameS,'_@x=',num2str(Xmeas),CheckId.SaveFig_type],'epsc')
            else
                saveas(gcf,[sf_savename,'PP3_timemean_',varnameS,'_@x=',num2str(Xmeas),CheckId.SaveFig_type])
            end
        end
    end
end
statusbarObj.setText('');

