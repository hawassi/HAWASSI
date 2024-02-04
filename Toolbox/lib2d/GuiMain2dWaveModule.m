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

function GuiMain2dWaveModule(inputdata)

%%%% Creating figure and panels
h.fig = figure('Name', strcat('AB2D-Wave',' [',inputdata.projectname,']'),'unit','normalized', 'Position',...
    [0.2 0.2 0.7 0.6],'NumberTitle', 'Off','MenuBar', 'None', 'Toolbar', 'None',...
    'Color', [.94 .94 .94]);
% [pathstr,~,~] = fileparts(mfilename('fullpath'));
% cd (pathstr);
% addpath(genpath(pathstr))
%
% h.pathnow=pathstr;%nputdata.pathnow;
% h.projechist=[];%inputdata.projecthist;
% h.module=[];%inputdata.module;
% h.projectname='';%inputdata.projectname;
% h.usernote='';%inputdata.usernote;
% h.projectdirectory=[];%inputdata.projectdirectory;

h.pathnow=inputdata.pathnow;
h.projechist=inputdata.projecthist;
h.module=inputdata.module;
h.moduleRestrict=inputdata.restrict;
h.projectname=inputdata.projectname;
h.usernote=inputdata.usernote;
h.projectdirectory=inputdata.projectdirectory;
h.flagOpenProj=inputdata.flagOpenProj;
if inputdata.flagOpenProj~=0
    h.GUIinput=inputdata.GUIinput;
end

ID_model_phiform=0; %% 1 formulation in phi otherwise in velocity

flag_warn_workdir=0;
%%%% define demo or full version
h.fullversion=1-inputdata.IdDemo;

%l for layout panel
try
    logo_axes = axes('parent', h.fig,'unit','normalized',...
        'Position', [.8 .915 0.2 0.08]);
    
    if ~isdeployed
        imshow([h.pathnow,'\Toolbox\lib2d\Tools\gui\images\logo\hawassi_black.jpg']); %
    else
        imshow('hawassi_black.jpg');
    end
catch
end


lmenu.main       = uix.HBoxFlex; %this is using toolbox gui layout
lmenu.menupanel  = uipanel('Parent',lmenu.main);
warning('off');
lmenu.tabgroup   = uitabgroup('Parent',lmenu.menupanel);
w = warning('query','last');  %turn off warning
warning('off',w.identifier);
lmenu.tabinput    = uitab(lmenu.tabgroup,'title','Input');
lmenu.tabpreview  = uitab(lmenu.tabgroup,'title','Preview');
lmenu.tabpostproc = uitab(lmenu.tabgroup,'title','Post-processing');
lmenu.tabIF = uitab(lmenu.tabgroup,'title','Internal flow');
%p for panel
%menu panel

menupanel = uiextras.CardPanel( 'Parent',lmenu.main);
% menupanel = uiextras.CardPanel( 'Parent', Figure3, 'Position', [.15 .15 .85 .85], 'Padding', 5 );
panel.project           = uipanel( 'Parent', menupanel);
panel.model             = uipanel( 'Parent', menupanel);
panel.spatial           = uipanel( 'Parent', menupanel);
panel.spatialdomain     = uipanel( 'Parent', menupanel);
panel.spatialbdy        = uipanel( 'Parent', menupanel);
panel.spatialbdyWall    = uipanel( 'Parent', menupanel);
panel.spatialbdyDamping = uipanel( 'Parent', menupanel);
panel.spatialbathy      = uipanel( 'Parent', menupanel);
panel.waveinput         = uipanel( 'Parent', menupanel);
panel.waveinput_ivp     = uipanel( 'Parent', menupanel);
panel.waveinput_influx  = uipanel( 'Parent', menupanel);
panel.waveinput_BdyAssim  = uipanel( 'Parent', menupanel);
panel.time              = uipanel( 'Parent', menupanel);
panel.options           = uipanel( 'Parent', menupanel);
menupanel.SelectedChild = 1;

%%%%%%Creating uitree%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmenu.project = uitreenode('v0', 'Project', 'Project', [], true);
tmenu.model = uitreenode('v0', 'Model', 'Model', [], true);
tmenu.spatial = uitreenode('v0', 'Spatial', 'Spatial', [], false);
tmenu.spatial.add(uitreenode('v0', 'Domain Area', 'Domain Area', [], true));
boundaries=uitreenode('v0', 'Boundaries', 'Boundaries', [], false);
boundaries.add(uitreenode('v0', 'Wall', 'Wall', [], true));
boundaries.add(uitreenode('v0', 'Fourier Bdry.', 'Fourier Bdry.', [], true));
tmenu.spatial.add(boundaries);
tmenu.spatial.add(uitreenode('v0', 'Bottom', 'Bottom', [], true));
tmenu.waveinput = uitreenode('v0', 'Wave Input', 'Wave Input', [], false);
tmenu.waveinput.add(uitreenode('v0', 'Initial Condition', 'Initial Condition',...
    [], true));
tmenu.waveinput.add(uitreenode('v0', 'Influx', 'Influx', [], true));
tmenu.waveinput.add(uitreenode('v0', 'Boundary Assimilation', 'Boundary Assimilation', [], true));
tmenu.time = uitreenode('v0', 'Time', 'Time', [], true);
tmenu.options = uitreenode('v0', 'Options', 'Options', [], true);

tmenu.root = uitreenode('v0', strcat('Input'), strcat('Input'), [], false);
tmenu.root.add(tmenu.project);
tmenu.root.add(tmenu.model);
tmenu.root.add(tmenu.spatial);
tmenu.root.add(tmenu.waveinput);
tmenu.root.add(tmenu.time);
tmenu.root.add(tmenu.options);


%%%%%panel tab preview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panel.prev_dispersion = uipanel( 'Parent', menupanel);
panel.prev_spatial    = uipanel( 'Parent', menupanel);
panel.prev_wave       = uipanel( 'Parent', menupanel);
panel.prev_log        = uipanel( 'Parent', menupanel);
menupanel.SelectedChild = 1;

tprev.dispersion = uitreenode('v0', 'Dispersion', 'Dispersion', [], true);
tprev.domain = uitreenode('v0', 'Spatial', 'Spatial', [], true);
tprev.wave = uitreenode('v0', 'Wave-input', 'Wave-input', [], true);
tprev.log = uitreenode('v0', 'Log file', 'Log file', [], true);

tprev.root = uitreenode('v0', strcat('Preview'), strcat('Preview'), [], false);
tprev.root.add(tprev.dispersion);
tprev.root.add(tprev.domain);
tprev.root.add(tprev.wave);
tprev.root.add(tprev.log);

%%%%%post proc tab preview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panel.pp_project   = uipanel( 'Parent', menupanel);
panel.pp_plot    = uipanel( 'Parent', menupanel);
panel.pp_plotdensity= uipanel( 'Parent', menupanel);
panel.pp_plotdensity_profile= uipanel( 'Parent', menupanel);
panel.pp_plotdensity_signal= uipanel( 'Parent', menupanel);
panel.pp_plotdensity_statistic= uipanel( 'Parent', menupanel);
panel.pp_plotline= uipanel( 'Parent', menupanel);
panel.pp_plotline_profile= uipanel( 'Parent', menupanel);
panel.pp_plotline_buoy= uipanel( 'Parent', menupanel);
panel.pp_plotline_statistic= uipanel( 'Parent', menupanel);
panel.pp_plotline_ham_mom= uipanel( 'Parent', menupanel);
panel.pp_plotline_breaking= uipanel( 'Parent', menupanel);
panel.pp_plotline_MTAA= uipanel( 'Parent', menupanel);
panel.pp_plotline_extreme= uipanel( 'Parent', menupanel);
panel.pp_quantitative_buoy = uipanel( 'Parent', menupanel);
panel.pp_anim      = uipanel( 'Parent', menupanel);
panel.pp_anim_density  = uipanel( 'Parent', menupanel);
panel.pp_anim_line     = uipanel( 'Parent', menupanel);
panel.pp_validation      = uipanel( 'Parent', menupanel);
panel.pp_validation_buoy = uipanel( 'Parent', menupanel);
panel.pp_validation_buoy_quant = uipanel( 'Parent', menupanel);

menupanel.SelectedChild = 1;

tpostproc.project = uitreenode('v0', 'Project', 'Project', [], true);
tpostproc.plotting = uitreenode('v0', 'Plotting', 'Plotting', [], false);
densityplot=uitreenode('v0', 'Density plots', 'Density plots', [], false);
densityplot.add(uitreenode('v0', 'Profile', 'Profile', [], true))
densityplot.add(uitreenode('v0', 'Signal', 'Signal', [], true))
densityplot.add(uitreenode('v0', 'Statistic', 'Statistic', [], true))
tpostproc.plotting.add(densityplot)
lineplot=uitreenode('v0', 'Line plots ', 'Line plots ', [], false);
lineplot.add(uitreenode('v0', 'Profile ', 'Profile ', [], true));
lineplot.add(uitreenode('v0', 'Buoy ', 'Buoy ', [], true));
lineplot.add(uitreenode('v0', 'Statistic1', 'Statistic ', [], true));
lineplot.add(uitreenode('v0', 'Hamiltonian & Momentum', 'Hamiltonian & Momentum', [], true));
lineplot.add(uitreenode('v0', 'Breaking events', 'Breaking events', [], true));
lineplot.add(uitreenode('v0', 'Max. Avg. Temp. Area Amplitude',  'Max. Avg. Temp. Area Amplitude', [], true));
lineplot.add(uitreenode('v0', 'Extreme events',  'Extreme events', [], true));
tpostproc.plotting.add(lineplot);
tpostproc.quantbuoy=uitreenode('v0', 'Quantitative (Buoy)', 'Quantitative (Buoy)', [], true);
tpostproc.animation= uitreenode('v0', 'Animation', 'Animation', [], false);
tpostproc.animation.add(uitreenode('v0', 'Density ', 'Density', [], true))
tpostproc.animation.add(uitreenode('v0', 'Line ', 'Line', [], true))
tpostproc.validation = uitreenode('v0', 'Validation', 'Validation', [],false);
tpostproc.validation.add(uitreenode('v0', 'Buoy', 'Buoy', [], true))
tpostproc.validation.add(uitreenode('v0', 'Quantitative', 'Quantitative', [], true))

tpostproc.root = uitreenode('v0', strcat('Post Processing'), strcat('Post Processing'), [], false);
tpostproc.root.add(tpostproc.project);
tpostproc.root.add(tpostproc.plotting);
tpostproc.root.add(tpostproc.quantbuoy);
tpostproc.root.add(tpostproc.animation);
tpostproc.root.add(tpostproc.validation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%panel Internal Flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panel.IF_project   = uipanel( 'Parent', menupanel);
panel.IF_calc      = uipanel( 'Parent', menupanel);
panel.IF_plotting   = uipanel( 'Parent', menupanel);
panel.IF_plotdensity   = uipanel( 'Parent', menupanel);
panel.IF_plotanimation = uipanel( 'Parent', menupanel);

menupanel.SelectedChild = 1;
tIF.project = uitreenode('v0', 'Project', 'Project', [], true);
tIF.calculate = uitreenode('v0', 'Calculation', 'Calculation', [], true);
tIF.plotting = uitreenode('v0', 'Plotting', 'Plotting', [], false);
IFdensityplot=uitreenode('v0', 'Density plots', 'Density plots', [], true);
IFanimplot=uitreenode('v0', 'Animation', 'Animation', [], true);
tIF.plotting.add(IFdensityplot)
tIF.plotting.add(IFanimplot)

tIF.root = uitreenode('v0', strcat('Post Processing'), strcat('Post Processing'), [], false);
tIF.root.add(tIF.project);
tIF.root.add(tIF.calculate);
tIF.root.add(tIF.plotting);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepanel = uix.VBoxFlex('Parent', h.fig, 'Padding', 5);
set(lmenu.main, 'width', [-1 -3], 'Spacing', 5, 'Parent', basepanel);
panel.monitor = uipanel( 'Parent', basepanel,'backgroundcolor',[0.8 0.8 0.8] );
set( basepanel, 'height',[-4 -1],'Spacing', 5,...
    'Units', 'Normalized', 'Position', [0 0 1 .9]);

set(lmenu.tabgroup ,'selectedtab',lmenu.tabIF);
set(lmenu.tabgroup ,'selectedtab',lmenu.tabpostproc);
set(lmenu.tabgroup ,'selectedtab',lmenu.tabpreview);
set(lmenu.tabgroup ,'selectedtab',lmenu.tabinput);

h.panel=panel;

%%%%%% UI controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%monitor box panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.monitorbox =uicontrol('Style','text','String','>> Initializing GUI, please wait ...', ...
    'HorizontalAlignment', 'left', 'Parent', panel.monitor,...
    'Unit', 'Normalized', 'Position', [0 0.15 0.8 0.8],...
    'backgroundcolor',[0.8 0.8 0.8],'Fontunits','normalized','FontSize', 0.15,...
    'fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%menu bar%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.menu_file = uimenu(h.fig,'Label','File');
h.menu_file_newproj=uimenu(h.menu_file,'Label','New project');
h.menu_file_saveproj=uimenu(h.menu_file,'Label','Save project');
h.menu_file_clearproj=uimenu(h.menu_file,'Label','Clear project');
h.menu_file_quitproj=uimenu(h.menu_file,'Label','Quit');


h.menu_tool = uimenu(h.fig,'Label','Tool');
h.menu_tool_calc = uimenu(h.menu_tool,'Label','Calculator');

h.menu_run = uimenu(h.fig,'Label','Run');
h.menu_run_preproc = uimenu(h.menu_run,'Label','Pre-processing');
h.menu_run_simul = uimenu(h.menu_run,'Label','Start simulation ...');


h.menu_help = uimenu(h.fig,'Label','Help');
h.menu_help_about = uimenu(h.menu_help,'Label','About');
h.menu_help_doc = uimenu(h.menu_help,'Label','Documentation');

[h,defdatwall]=handles_maingui(panel,h,ID_model_phiform);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%set selection function of menu list %%%%%%%%%%%%%%%%%%%%%%%%
[mtree, container] = uitree('v0', 'Root', tmenu.root, 'Parent',  lmenu.tabinput...
    ,'SelectionChangeFcn',{@menuSelectFcn,h});
mtree.expand(mtree.getRoot);
mtree.expand(tmenu.spatial);
mtree.expand(tmenu.waveinput);
mtree.expand(boundaries);
mtree.setSelectedNode(tmenu.root)
set(container, 'Parent', lmenu.tabinput);
set(container, 'Unit', 'normalized', 'Position', [0 0 1 1]);

[mtree_prev, container_prev] = uitree('v0', 'Root', tprev.root, 'Parent', ...
    lmenu.tabpreview, 'SelectionChangeFcn',@mySelectFcn_preview);
mtree_prev.expand(mtree_prev.getRoot);
mtree_prev.setSelectedNode(tprev.root);
set(container_prev, 'Parent', lmenu.tabpreview);
set(container_prev, 'Unit', 'normalized', 'Position', [0 0 1 1]);

[mtree_postproc, container_postproc] = uitree('v0', 'Root', tpostproc.root, 'Parent', ...
    lmenu.tabpostproc, 'SelectionChangeFcn',{@mySelectFcn_postproc,h});
mtree_postproc.expand(mtree_postproc.getRoot);
mtree_postproc.expand(tpostproc.plotting);
mtree_postproc.expand(densityplot);
mtree_postproc.expand(lineplot);
mtree_postproc.expand(tpostproc.animation);
mtree_postproc.expand(tpostproc.validation);
mtree_postproc.setSelectedNode(tpostproc.root);
set(container_postproc, 'Parent', lmenu.tabpostproc);
set(container_postproc, 'Unit', 'normalized', 'Position', [0 0 1 1]);

[mtree_IF, container_IF] = uitree('v0', 'Root', tIF.root, 'Parent', ...
    lmenu.tabIF, 'SelectionChangeFcn',{@mySelectFcn_IF,h});
mtree_IF.expand(mtree_IF.getRoot);
mtree_IF.expand(tIF.plotting);

set(container_IF, 'Parent', lmenu.tabIF);
set(container_IF, 'Unit', 'normalized', 'Position', [0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%set callback handles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h.fig,'ResizeFcn',{@guiAB2dwaveResizeFcn,h});
set(h.fig,'CloseRequestFcn',{@callback_menu_file_quitproj,h});
set(h.menu_file_newproj, 'callback',{@callback_menu_file_newproj,h});
set(h.menu_file_saveproj, 'callback',{@callback_menu_file_saveproj,h});
set(h.menu_file_clearproj, 'callback',{@callback_menu_file_clearproj,h});
set(h.menu_file_quitproj, 'callback',{@callback_menu_file_quitproj,h});
set(h.menu_run_preproc,'callback',{@callback_preproc,h});
set(h.menu_run_simul,'callback',{@callback_run_odesolver,h});
set(h.menu_tool_calc,'callback',{@callback_calculator})
set(h.menu_help_about,'callback',{@callback_about_AB2})
set(h.menu_help_doc,'callback',{@callback_doc})
set(h.project_edit_name,'callback',{@callback_projectname,h})
set(h.project_button_projdir,'callback',{@callback_project_button_projdir,h})
set(h.model_popup_dynamic,'callback',{@callback_modeldyn,h});
set(h.model_popup_dispersion,'Callback',{@callback_dispersionmenu,h});
set(h.model_popup_breaking,'Callback',{@callback_breaking,h});
set(h.model_break_edit_initiation,'Callback',{@callback_breaking_initiation,h});
set(h.model_break_edit_termination,'Callback',{@callback_breaking_termination,h});
set(h.model_break_edit_Tstar,'Callback',{@callback_breaking_Tstar,h});
set(h.model_current_cb,'Callback',{@callback_current_cb,h});
set(h.model_current_ux_edit,'Callback',{@callback_current_ux,h});
set(h.model_current_uy_edit,'Callback',{@callback_current_uy,h})
set(h.spatial_edit_xmin,'Callback',{@callback_spatial_xmin,h});
set(h.spatial_edit_xmax,'Callback',{@callback_spatial_xmax,h});
set(h.spatial_edit_ymin,'Callback',{@callback_spatial_ymin,h});
set(h.spatial_edit_ymax,'Callback',{@callback_spatial_ymax,h});
set(h.spatial_edit_dx,'Callback',{@callback_spatial_dx,h});
set(h.spatial_edit_dy,'Callback',{@callback_spatial_dy,h});
set(h.cutfrac_k,'Callback',{@callback_cutfrac_k,h});
set(h.wall_popup,'callback',{@callback_wall_popup,h});
set(h.wall_table,'CellEditCallback',{@callback_wall_table,h},...
    'CellSelectionCallback',{@tableDelSelection,h,h.wall_button_deleterow});
set(h.wall_button_addrow,'callback',{@callback_addrow_table,h,h.wall_table,defdatwall});
set(h.wall_button_deleterow,'callback',{@callback_deleterow_table,h,h.wall_table})


set(h.damping_popup,'callback',{@callback_damping_popup,h});
set(h.damping_table,'CellEditCallback',{@callback_damping_table,h},...
    'CellSelectionCallback',{@tableDelSelection,h,h.damping_button_deleterow});
set(h.damping_button_addrow,'callback',{@callback_damping_addrow,h});
set(h.damping_button_deleterow,'callback',{@callback_damping_deleterow,h})

set(h.bathymetry_popup_type,'callback',{@callback_bathy_popup,h})
set(h.bathymetry_edit_depth,'callback',{@callback_bathy_depth,h})
set(h.bathymetry_edit_mindepth,'callback',{@callback_bathy_mindepth,h});
set(h.bathymetry_edit_middepth,'callback',{@callback_bathy_middepth,h}); % nunu
set(h.bathymetry_edit_maxdepth,'callback',{@callback_bathy_maxdepth,h});
set(h.bathymetry_edit_slope,'callback',{@callback_bathy_slope,h});
set(h.bathymetry_edit_startslope,'callback',{@callback_bathy_startslope,h});
set(h.bathymetry_edit_mindepthshore,'callback',{@callback_bathy_mindepthshore,h});
set(h.bathymetry_edit_maxdepthshore,'callback',{@callback_bathy_maxdepthshore,h});
set(h.bathymetry_edit_slopeshore,'callback',{@callback_bathy_slopeshore,h})
set(h.bathymetry_edit_shoreposition,'callback',{@callback_bathy_shoreposition,h})
set(h.bathymetry_button_load,'callback',{@callback_bathy_userdefined,h});
set(h.bathymetry_popup_interpdepth,'callback',{@callback_bathy_popup_interpdepth,h}); % nunu
set(h.bathymetry_checkbox_friction,'callback',{@callback_friction_check,h});
set(h.bathymetry_table_friction,'CellEditCallback',{@callback_friction_table,h},...
    'CellSelectionCallback',{@tableDelSelection,h,h.bathymetry_button_deleterow});
set(h.bathymetry_button_addrow,'callback',{@callback_friction_addrow,h});
set(h.bathymetry_button_deleterow,'callback',{@callback_deleterow_table,h,h.bathymetry_table_friction});
set(h.waveinput_ivp_popup_type,'callback',{@callback_ivp_popup,h});
set(h.waveinput_ivp_edit_A,'callback',{@callback_ivp_amplitude,h});
set(h.waveinput_ivp_edit_centre_position,'callback',{@callback_ivp_centreposition,h});
set(h.waveinput_ivp_edit_sd,'callback',{@callback_ivp_stdev,h});
set(h.waveinput_ivp_button_load,'callback',{@callback_ivp_loaddata,h});
set(h.waveinput_influx_popup,'callback',{@callback_influx_popup,h});
set(h.waveinput_influx_button_addrow,'callback',{@callback_influx_addrow,h});
set(h.waveinput_influx_button_deleterow,'callback',{@callback_influx_deleterow,h});
set(h.waveinput_influx_checkbox_ramp,'callback',{@callback_influx_check_ramp,h})
set(h.waveinput_influx_edit_ramp,'callback',{@callback_influx_ramp,h});
set(h.waveinput_influxline_checkbox_ramp,'callback',{@callback_influxline_check_ramp,h})
set(h.waveinput_influxline_edit_ramp,'callback',{@callback_influxline_ramp,h});
set(h.waveinput_influx_checkbox_nonlinadj,'callback',{@callback_influx_check_nonlinadj,h})
set(h.waveinput_influx_edit_nonlinadj,'callback',{@callback_influx_nonlinadj,h});
set(h.waveinput_influx_proptable,'CellEditCallback',{@callback_influx_prop_table,h},...
    'CellSelectionCallback',{@tableDelSelection,h,h.waveinput_influx_button_deleterow});
set(h.waveinput_influx_methodtable,'CellEditCallback',{@callback_influx_method_table,h},...
     'CellSelectionCallback',{@tableDelSelection,h,h.waveinput_influx_button_deleterow});

set(h.waveinput_bdy_assim_popup,'callback',{@callback_bdyassim_popup,h});
set(h.waveinput_bdy_assim_shape_popup,'callback',{@callback_bdyassim_shape,h})
set(h.waveinput_bdy_assim_edit_R1,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_edit_xc,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_edit_yc,'callback',{@callback_edit_param,h})

set(h.waveinput_bdy_assim_edit_smooth,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_time_edit_interval_init,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_time_edit_interval_end,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_time_edit_step,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_button_load,'callback',{@callback_bdyaasim_loaddata,h,...
    h.waveinput_bdy_assim_edit_data})
if ID_model_phiform==1
set(h.waveinput_bdy_assim_cb_phi,'callback',{@callback_bdyaasim_cb_phi,h})
set(h.waveinput_bdy_assim_button_load_phi,'callback',{@callback_bdyaasim_loaddata,h,...
    h.waveinput_bdy_assim_edit_data_phi})
else
set(h.waveinput_bdy_assim_cb_vel,'callback',{@callback_bdyaasim_cb_vel,h})
set(h.waveinput_bdy_assim_button_load_u,'callback',{@callback_bdyaasim_loaddata_vel,h,...
    h.waveinput_bdy_assim_edit_data_vel})    
set(h.waveinput_bdy_assim_button_load_v,'callback',{@callback_bdyaasim_loaddata_vel,h,...
    h.waveinput_bdy_assim_edit_data_vel})    
end

set(h.waveinput_bdy_assim_cb_nonlinAdj,'callback',{@callback_bdyassim_nonlinadj,h})
set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'callback',{@callback_edit_param,h})
set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'callback',{@callback_edit_param,h})

set(h.time_edit_interval_init,'callback',{@callback_time_tinit,h});
set(h.time_edit_interval_end,'callback',{@callback_time_tend,h});
set(h.time_edit_step,'callback',{@callback_time_dt,h})
set(h.options_checkbox_intflow,'callback',{@callback_option_inflow_check,h})
set(h.options_edit_interval_init,'callback',{@callback_option_inflow_tinit,h});
set(h.options_edit_interval_end,'callback',{@callback_option_inflow_tend,h});
set(h.options_edit_step,'callback',{@callback_option_inflow_dtstep,h});
set(h.options_edit_partition,'callback',{@callback_option_partition,h});
set(h.options_checkbox_default,'callback',{@callback_ode_partition,h})
set(h.options_edit_ode_tol,'callback',{@callback_ode_tol,h});

set(h.options_montecarlo_cb_check,'callback',{@callback_mc_cb_check,h});
set(h.options_montecarlo_cb_numset,'callback',{@callback_mc_cb_numset,h});
set(h.options_montecarlo_cb_Nsimul,'callback',{@callback_mc_N_run,h});
set(h.options_montecarlo_cb_waveinput,'callback',{@callback_mc_cb_waveinput,h})
set(h.options_montecarlo_popup_waveinput,'callback',{@callback_mc_popup_waveinput,h})
set(h.options_montecarlo_cb_A,'callback',{@callback_mc_cb_param,h.options_montecarlo_edit_A});
set(h.options_montecarlo_cb_Tp,'callback',{@callback_mc_cb_param,h.options_montecarlo_edit_Tp});
set(h.options_montecarlo_cb_gamma,'callback',{@callback_mc_cb_param,h.options_montecarlo_edit_gamma});
set(h.options_montecarlo_cb_s,'callback',{@callback_mc_cb_param,h.options_montecarlo_edit_s});
set(h.options_montecarlo_edit_A,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_Tp,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_gamma,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_s,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_Nsimul,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_px,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_edit_py,'callback',{@callback_mc_edit_param});
set(h.options_montecarlo_pb_waveinput,'callback',{@callback_mc_pb_storedata,h})

set(h.preview.dispersion_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.preview.dispersion_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.preview.dispersion_button_datacursor,'callback',{@callback_datacursor,h})
set(h.preview.spatial_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.preview.spatial_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.preview.spatial_button_datacursor,'callback',{@callback_datacursor,h})
set(h.preview.spatial_popup_var,'callback',{@callback_preview_spatial,h})

set(h.preview.wave_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.preview.wave_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.preview.wave_button_datacursor,'callback',{@callback_datacursor,h})
set(h.preview.wave_popup_var,'callback',{@callback_preview_wave,h})
set(h.preview.wave_popup_spect,'callback',{@callback_preview_wave_spectrum,h})

set(h.pp_project_edit_name,'callback',{@callback_pp_projectname,h})
set(h.pp_project_button_projdir,'callback',{@callback_pp_project_button_projdir,h})
set(h.pp_project_load_data,'callback',{@callback_pp_loaddata,h});
set(h.pp_project_toggle_customize_data,'callback',{@callback_pp_customize_data,h});

set(h.pp_density_profile_plot_button,'callback',{@callback_pp_plot_density_profile,h})
set(h.pp_density_profile_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.pp_density_profile_setting_edit_view})
set(h.pp_density_profile_setting_edit_view,'callback',{@callback_pp_setting_edit_view,h})
set(h.pp_density_profile_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.pp_density_profile_setting_edit_clim})
set(h.pp_density_profile_setting_edit_clim,'callback',{@callback_pp_setting_edit_clim,h})
set(h.pp_density_profile_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_density_profile_setting_edit_xlim})
set(h.pp_density_profile_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h})
set(h.pp_density_profile_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_density_profile_setting_edit_ylim})
set(h.pp_density_profile_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h})
set(h.pp_density_profile_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_density_profile_setting_edit_coarse})
set(h.pp_density_profile_setting_edit_coarse,'callback',{@callback_pp_setting_edit_coarse,h})
set(h.pp_density_profile_setting_cb_ampliref,'callback',{@callback_pp_setting_cb_ampliref,...
    h.pp_density_profile_setting_edit_ampliref})
set(h.pp_density_profile_setting_edit_ampliref,'callback',{@callback_pp_setting_edit_ampliref,h});

set(h.pp_density_profile_setting_cb_level,'callback',...
    {@callback_pp_setting_cb_level,h,h.pp_density_profile_setting_edit_level,...
    h.pp_density_profile_cb_level_eta,h.pp_density_profile_cb_level_phi,...
    h.pp_density_profile_cb_level_quiver,h.pp_density_profile_cb_level_extremeCrest,...
    h.pp_density_profile_cb_level_extremeTrough})
set(h.pp_density_profile_setting_edit_level,'callback',{@callback_pp_setting_edit_level,h})
set(h.pp_density_profile_edit_bathyScale,'callback',{@callback_anim_edit_bathy_scale})
set(h.pp_density_profile_cb_bathy,'callback',{@callback_anim_cb_bathy,h.pp_density_profile_edit_bathyScale})



set(h.pp_density_profile_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_density_profile_setting_popup_savefig})
set(h.pp_density_profile_time_edit,'callback',{@callback_density_profile_at_time,h})
set(h.pp_density_profile_cb_level_eta,'callback',{@callback_density_profile_level_eta,...
    h,h.pp_density_profile_cb_level_phi,h.pp_density_profile_cb_level_quiver,...
    h.pp_density_profile_setting_cb_level,h.pp_density_profile_setting_edit_level})
set(h.pp_density_profile_cb_level_phi,'callback',{@callback_density_profile_level_phi,...
    h,h.pp_density_profile_cb_level_eta,h.pp_density_profile_cb_level_quiver,...
    h.pp_density_profile_setting_cb_level,h.pp_density_profile_setting_edit_level})
set(h.pp_density_profile_cb_level_quiver,'callback',{@callback_density_profile_quiver,...
    h.pp_density_profile_cb_level_eta,h.pp_density_profile_cb_level_phi,...
    h.pp_density_profile_setting_cb_level,h.pp_density_profile_setting_edit_level})
set(h.pp_density_profile_cb_level_extremeCrest,'callback',...
    {@callback_density_profile_level_extremecrest,...
    h.pp_density_profile_setting_cb_level,h.pp_density_profile_setting_edit_level,...
    h.pp_density_profile_setting_cb_ampliref,h.pp_density_profile_setting_edit_ampliref})
set(h.pp_density_profile_cb_level_extremeTrough,'callback',...
    {@callback_density_profile_level_extremecrest,....
    h.pp_density_profile_setting_cb_level,h.pp_density_profile_setting_edit_level,...
    h.pp_density_profile_setting_cb_ampliref,h.pp_density_profile_setting_edit_ampliref})
set(h.pp_density_profile_popup_var,'callback',{@callback_density_profile_popup_var,h,...
    h.pp_density_profile_cb_level_eta,h.pp_density_profile_cb_level_phi,...
    h.pp_density_profile_cb_level_quiver,h.pp_density_profile_cb_level_extremeCrest,...
    h.pp_density_profile_setting_cb_level})
set(h.pp_density_button_zoomin,'callback',{@callback_zoom_in,h});
set(h.pp_density_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_density_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_density_button_rotate,'callback',{@callback_rotate,h})
set(h.pp_density_button_pan,'callback',{@callback_pan,h})


set(h.pp_density_signal_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_density_signal_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_density_signal_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_density_signal_button_rotate,'callback',{@callback_rotate,h})
set(h.pp_density_signal_button_pan,'callback',{@callback_pan,h})

set(h.pp_density_signal_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_density_signal_x_cb,'callback',...
    {@callback_pp_cb_x,h,h.pp_density_signal_x_edit,h.pp_density_signal_y_cb});
set(h.pp_density_signal_x_edit,'callback',{@callback_pp_edit_x,h});
set(h.pp_density_signal_y_cb,'callback',...
    {@callback_pp_cb_y,h,h.pp_density_signal_y_edit,h.pp_density_signal_x_cb});
set(h.pp_density_signal_y_edit,'callback',{@callback_pp_edit_y,h});
set(h.pp_density_signal_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.pp_density_signal_setting_edit_view});
set(h.pp_density_signal_setting_edit_view,'callback',...
    {@callback_pp_setting_edit_view,h});
set(h.pp_density_signal_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.pp_density_signal_setting_edit_clim});
set(h.pp_density_signal_setting_edit_clim,'callback',...
    {@callback_pp_setting_edit_clim,h});
set(h.pp_density_signal_setting_cb_spatlim,'callback',...
    {@callback_pp_setting_cb_spatlim,h.pp_density_signal_setting_edit_spatlim});
set(h.pp_density_signal_setting_edit_spatlim,'callback',{@callback_pp_setting_edit_spatlim,h});
set(h.pp_density_signal_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_density_signal_setting_edit_tlim});
set(h.pp_density_signal_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_density_signal_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_density_signal_setting_edit_coarse});
set(h.pp_density_signal_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_density_signal_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_density_signal_setting_popup_savefig})
set(h.pp_density_signal_plot_button,'callback',{@callback_pp_plot_density_signal,h})

set(h.pp_density_statistic_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_density_statistic_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_density_statistic_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_density_statistic_button_rotate,'callback',{@callback_rotate,h})
set(h.pp_density_statistic_button_pan,'callback',{@callback_pan,h})

set(h.pp_density_statistic_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.pp_density_statistic_setting_edit_view});
set(h.pp_density_statistic_setting_edit_view,'callback',...
    {@callback_pp_setting_edit_view,h});
set(h.pp_density_statistic_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.pp_density_statistic_setting_edit_clim});
set(h.pp_density_statistic_setting_edit_clim,'callback',...
    {@callback_pp_setting_edit_clim,h});
set(h.pp_density_statistic_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_density_statistic_setting_edit_xlim});
set(h.pp_density_statistic_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h});
set(h.pp_density_statistic_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_density_statistic_setting_edit_ylim});
set(h.pp_density_statistic_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h});
set(h.pp_density_statistic_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_density_statistic_setting_edit_tlim});
set(h.pp_density_statistic_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_density_statistic_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_density_statistic_setting_edit_coarse});
set(h.pp_density_statistic_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_density_statistic_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_density_statistic_setting_popup_savefig})
% set(h.pp_density_statistic_hs_cb,'callback',{@callback_cb_hs,h});
% set(h.pp_density_statistic_mtc_cb,'callback',{@callback_cb_mtc,h});
% set(h.pp_density_statistic_mtt_cb,'callback',{@callback_cb_mtt,h});
% set(h.pp_density_statistic_mwl_cb,'callback',{@callback_cb_mwl,h});
set(h.pp_density_statistic_plot_button,'callback',{@callback_pp_plot_density_statistic,h})

set(h.pp_line_profile_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_profile_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_profile_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_profile_button_pan,'callback',{@callback_pan,h})

set(h.pp_line_profile_edit_bathyScale,'callback',{@callback_anim_edit_bathy_scale})
set(h.pp_line_profile_cb_bathy,'callback',{@callback_anim_cb_bathy,h.pp_line_profile_edit_bathyScale})

set(h.pp_line_profile_time_edit,'callback',{@callback_line_profile_at_time,h})
set(h.pp_line_profile_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_line_profile_x_cb,'callback',...
    {@callback_pp_cb_x,h,h.pp_line_profile_x_edit,h.pp_line_profile_y_cb});
set(h.pp_line_profile_x_edit,'callback',{@callback_pp_edit_x,h});
set(h.pp_line_profile_y_cb,'callback',...
    {@callback_pp_cb_y,h,h.pp_line_profile_y_edit,h.pp_line_profile_x_cb});
set(h.pp_line_profile_y_edit,'callback',{@callback_pp_edit_y,h});

set(h.pp_line_profile_MTC_cb,'callback',{@callback_pp_MTA_check,h})
set(h.pp_line_profile_MTT_cb,'callback',{@callback_pp_MTA_check,h})
set(h.pp_line_profile_setting_cb_MTAtlim,'callback',{@callback_pp_MTA_tlim_check,h})
set(h.pp_line_profile_setting_edit_MTAtlim,'callback',{@callback_pp_setting_edit_tlim,h})
set(h.pp_line_profile_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_profile_setting_edit_zlim})
set(h.pp_line_profile_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_profile_setting_cb_spatlim,'callback',...
    {@callback_pp_setting_cb_spatlim,h.pp_line_profile_setting_edit_spatlim});
set(h.pp_line_profile_setting_edit_spatlim,'callback',{@callback_pp_setting_edit_spatlim,h});
set(h.pp_line_profile_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_profile_setting_edit_coarse});
set(h.pp_line_profile_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_profile_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_profile_setting_popup_savefig})
set(h.pp_line_profile_plot_button,'callback',{@callback_pp_plot_line_profile,h})


set(h.pp_line_buoy_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_buoy_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_buoy_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_buoy_button_pan,'callback',{@callback_pan,h})
set(h.pp_line_buoy_x_edit,'callback',{@callback_pp_line_buoy_edit_x,h});
set(h.pp_line_buoy_y_edit,'callback',{@callback_pp_line_buoy_edit_y,h});
set(h.pp_line_buoy_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_line_buoy_combine_tinterv_edit,'callback',{@callback_edit_param,h})
set(h.pp_line_buoy_cb_combine_buoys,'callback',{@callback_cb_combine_buoys,h})

set(h.pp_line_buoy_cb_signal,'callback',...
    {@callback_pp_cb_signal,h.pp_line_buoy_cb_spectrum,...
    h.pp_line_buoy_popup_spectrum_var,h.pp_line_buoy_setting_cb_spsmooth,...
    h.pp_line_buoy_setting_edit_spsmooth})
set(h.pp_line_buoy_cb_spectrum,'callback',{@callback_pp_cb_spectrum,...
    h.pp_line_buoy_popup_spectrum_var,h.pp_line_buoy_cb_signal,...
    h.pp_line_buoy_setting_cb_spsmooth,h.pp_line_buoy_setting_edit_spsmooth})

set(h.pp_line_buoy_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_buoy_setting_edit_zlim})
set(h.pp_line_buoy_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_buoy_setting_cb_horznlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_buoy_setting_edit_horznlim});
set(h.pp_line_buoy_setting_edit_horznlim,'callback',{@callback_pp_setting_edit_horzlim,h});
set(h.pp_line_buoy_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_buoy_setting_edit_tlim});
set(h.pp_line_buoy_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_buoy_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_buoy_setting_edit_coarse});
set(h.pp_line_buoy_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_buoy_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_buoy_setting_popup_savefig})
set(h.pp_line_buoy_setting_cb_spsmooth,'callback',{@callback_pp_setting_cb_coarse,...
    h.pp_line_buoy_setting_edit_spsmooth})
set(h.pp_line_buoy_setting_edit_spsmooth,'callback',{@callback_edit_param,h})
set(h.pp_line_buoy_plot_button,'callback',{@callback_pp_plot_line_buoy,h})


set(h.pp_line_statistic_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_statistic_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_statistic_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_statistic_button_pan,'callback',{@callback_pan,h})

set(h.pp_line_statistic_along_axes_cb,'callback',{@callback_pp_statistic_cb_along_axes,h})
set(h.pp_line_statistic_x_cb,'callback',...
    {@callback_pp_cb_x,h,h.pp_line_statistic_x_edit,h.pp_line_statistic_y_cb});
set(h.pp_line_statistic_x_edit,'callback',{@callback_pp_edit_x,h});
set(h.pp_line_statistic_y_cb,'callback',...
    {@callback_pp_cb_y,h,h.pp_line_statistic_y_edit,h.pp_line_statistic_x_cb});
set(h.pp_line_statistic_y_edit,'callback',{@callback_pp_edit_y,h});
set(h.pp_line_statistic_buoys_cb,'callback',{@callback_pp_statistic_cb_at_buoys,h})
set(h.pp_line_statistic_buoy_x_edit,'callback',{@callback_edit_param,h})
set(h.pp_line_statistic_buoy_y_edit,'callback',{@callback_edit_param,h})

set(h.pp_line_statistic_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_statistic_setting_edit_zlim})
set(h.pp_line_statistic_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_statistic_setting_cb_spatlim,'callback',...
    {@callback_pp_setting_cb_spatlim,h.pp_line_statistic_setting_edit_spatlim});
set(h.pp_line_statistic_setting_edit_spatlim,'callback',{@callback_pp_setting_edit_spatlim,h});
set(h.pp_line_statistic_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_statistic_setting_edit_tlim});
set(h.pp_line_statistic_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_statistic_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_statistic_setting_edit_coarse});
set(h.pp_line_statistic_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_statistic_setting_cb_Hs_ref,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_statistic_setting_edit_Hs_ref})
set(h.pp_line_statistic_setting_edit_Hs_ref,'callback',{@callback_edit_param,h});
set(h.pp_line_statistic_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_statistic_setting_popup_savefig})
set(h.pp_line_statistic_plot_button,'callback',{@callback_pp_plot_line_statistic,h})
set(h.pp_line_statistic_popup_var,'callback',{@callback_pp_plot_line_statistic_var,h})


set(h.pp_line_ham_mom_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_ham_mom_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_ham_mom_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_ham_mom_button_pan,'callback',{@callback_pan,h})


set(h.pp_line_ham_mom_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_ham_mom_setting_edit_zlim})
set(h.pp_line_ham_mom_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_ham_mom_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_line_ham_mom_setting_edit_xlim});
set(h.pp_line_ham_mom_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h});
set(h.pp_line_ham_mom_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_line_ham_mom_setting_edit_ylim});
set(h.pp_line_ham_mom_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h});
set(h.pp_line_ham_mom_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_ham_mom_setting_edit_tlim});
set(h.pp_line_ham_mom_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_ham_mom_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_ham_mom_setting_edit_coarse});
set(h.pp_line_ham_mom_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_ham_mom_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_ham_mom_setting_popup_savefig})
set(h.pp_line_ham_mom_plot_button,'callback',{@callback_pp_plot_line_ham_mom,h})

set(h.pp_line_breaking_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_breaking_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_breaking_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_breaking_button_pan,'callback',{@callback_pan,h})


set(h.pp_line_breaking_popup_var,'callback',{@callback_line_breaking_popup_var,h})
set(h.pp_line_breaking_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_breaking_setting_edit_zlim})
set(h.pp_line_breaking_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_breaking_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_line_breaking_setting_edit_xlim});
set(h.pp_line_breaking_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h});
set(h.pp_line_breaking_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_line_breaking_setting_edit_ylim});
set(h.pp_line_breaking_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h});
set(h.pp_line_breaking_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_breaking_setting_edit_tlim});
set(h.pp_line_breaking_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_breaking_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_breaking_setting_edit_coarse});
set(h.pp_line_breaking_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_breaking_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_breaking_setting_popup_savefig})
set(h.pp_line_breaking_plot_button,'callback',{@callback_pp_plot_line_breaking,h})

set(h.pp_line_MTAA_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_line_MTAA_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_MTAA_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_MTAA_button_pan,'callback',{@callback_pan,h})


set(h.pp_line_MTAA_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_MTAA_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_MTAA_setting_edit_zlim})
set(h.pp_line_MTAA_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_MTAA_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_line_MTAA_setting_edit_xlim});
set(h.pp_line_MTAA_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h});
set(h.pp_line_MTAA_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_line_MTAA_setting_edit_ylim});
set(h.pp_line_MTAA_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h});
set(h.pp_line_MTAA_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_MTAA_setting_edit_tlim});
set(h.pp_line_MTAA_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_MTAA_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_MTAA_setting_edit_coarse});
set(h.pp_line_MTAA_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_MTAA_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_MTAA_setting_popup_savefig})
set(h.pp_line_MTAA_plot_button,'callback',{@callback_pp_plot_line_MTAA,h})

set(h.pp_line_Extreme_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_line_Extreme_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_line_Extreme_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_line_Extreme_button_rotate,'callback',{@callback_rotate,h})
set(h.pp_line_Extreme_button_pan,'callback',{@callback_pan,h})

set(h.pp_line_Extreme_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_line_Extreme_setting_edit_zlim})
set(h.pp_line_Extreme_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_line_Extreme_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_line_Extreme_setting_edit_xlim});
set(h.pp_line_Extreme_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h});
set(h.pp_line_Extreme_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_line_Extreme_setting_edit_ylim});
set(h.pp_line_Extreme_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h});
set(h.pp_line_Extreme_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_line_Extreme_setting_edit_tlim});
set(h.pp_line_Extreme_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h});
set(h.pp_line_Extreme_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_line_Extreme_setting_edit_coarse});
set(h.pp_line_Extreme_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_line_Extreme_setting_cb_ampliref,'callback',{@callback_pp_setting_cb_ampliref,...
    h.pp_line_Extreme_setting_edit_ampliref})
set(h.pp_line_Extreme_setting_edit_ampliref,'callback',{@callback_pp_setting_edit_ampliref,h});
set(h.pp_line_Extreme_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_line_Extreme_setting_popup_savefig})
set(h.pp_line_Extreme_plot_button,'callback',{@callback_pp_plot_line_Extreme,h})


set(h.pp_quant_buoy_x_edit,'callback',{@callback_edit_param,h});
set(h.pp_quant_buoy_y_edit,'callback',{@callback_edit_param,h});
set(h.pp_quant_buoy_tinterv_edit,'callback',{@callback_edit_param,h})
set(h.pp_quant_buoy_button,'callback',{@callback_quant_buoy_calculate,h});

set(h.pp_anim_density_plot_button,'callback',{@callback_pp_plot_anim_density,h})
set(h.pp_anim_density_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.pp_anim_density_setting_edit_view})
set(h.pp_anim_density_setting_edit_view,'callback',{@callback_pp_setting_edit_view,h})
set(h.pp_anim_density_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.pp_anim_density_setting_edit_clim})
set(h.pp_anim_density_setting_edit_clim,'callback',{@callback_pp_setting_edit_clim,h})
set(h.pp_anim_density_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.pp_anim_density_setting_edit_xlim})
set(h.pp_anim_density_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h})
set(h.pp_anim_density_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.pp_anim_density_setting_edit_ylim})
set(h.pp_anim_density_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h})
set(h.pp_anim_density_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_anim_density_setting_edit_tlim})
set(h.pp_anim_density_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h})
set(h.pp_anim_density_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_anim_density_setting_edit_coarse})
set(h.pp_anim_density_setting_edit_coarse,'callback',{@callback_pp_setting_edit_coarse,h})
set(h.pp_anim_density_setting_cb_ampliref,'callback',{@callback_pp_setting_cb_ampliref,...
    h.pp_anim_density_setting_edit_ampliref})
set(h.pp_anim_density_setting_edit_ampliref,'callback',{@callback_pp_setting_edit_ampliref,h});
set(h.pp_anim_density_setting_cb_level,'callback',...
    {@callback_pp_setting_cb_level,h,h.pp_anim_density_setting_edit_level,...
    h.pp_anim_density_cb_level_eta,h.pp_anim_density_cb_level_phi,...
    h.pp_anim_density_cb_level_quiver,h.pp_anim_density_cb_level_extremeCrest,...
    h.pp_anim_density_cb_level_extremeTrough})
set(h.pp_anim_density_setting_edit_level,'callback',{@callback_pp_setting_edit_level,h})
set(h.pp_anim_density_setting_cb_saveanim,'callback',...
    {@callback_pp_setting_cb_saveanim,h.pp_anim_density_setting_edit_gifset,h.pp_anim_density_setting_text_gifset})
set(h.pp_anim_density_setting_edit_gifset,'callback',{@callback_pp_setting_edit_gifset,h})
set(h.pp_anim_density_cb_level_eta,'callback',{@callback_density_profile_level_eta,...
    h,h.pp_anim_density_cb_level_phi,h.pp_anim_density_cb_level_quiver,...
    h.pp_anim_density_setting_cb_level,h.pp_anim_density_setting_edit_level})
set(h.pp_anim_density_cb_level_phi,'callback',{@callback_density_profile_level_phi,...
    h,h.pp_anim_density_cb_level_eta,h.pp_anim_density_cb_level_quiver,...
    h.pp_anim_density_setting_cb_level,h.pp_anim_density_setting_edit_level})
set(h.pp_anim_density_cb_level_quiver,'callback',{@callback_density_profile_quiver,...
    h.pp_anim_density_cb_level_eta,h.pp_anim_density_cb_level_phi,...
    h.pp_anim_density_setting_cb_level,h.pp_anim_density_setting_edit_level})
set(h.pp_anim_density_cb_level_extremeCrest,'callback',...
    {@callback_density_profile_level_extremecrest,h.pp_anim_density_setting_cb_level,...
    h.pp_anim_density_setting_edit_level,h.pp_anim_density_setting_cb_ampliref,...
    h.pp_anim_density_setting_edit_ampliref})
set(h.pp_anim_density_cb_level_extremeTrough,'callback',...
    {@callback_density_profile_level_extremecrest,h.pp_anim_density_setting_cb_level,...
    h.pp_anim_density_setting_edit_level,h.pp_anim_density_setting_cb_ampliref,...
    h.pp_anim_density_setting_edit_ampliref})
set(h.pp_anim_density_popup_var,'callback',{@callback_density_profile_popup_var,h,...
    h.pp_anim_density_cb_level_eta,h.pp_anim_density_cb_level_phi,...
    h.pp_anim_density_cb_level_quiver,h.pp_anim_density_cb_level_extremeCrest,...
    h.pp_anim_density_setting_cb_level})
set(h.pp_anim_density_edit_bathyScale,'callback',{@callback_anim_edit_bathy_scale})
set(h.pp_anim_density_cb_bathy,'callback',{@callback_anim_cb_bathy,h.pp_anim_density_edit_bathyScale})


set(h.pp_anim_play,'Callback',{@callback_resume,h.pp_anim_density_axes,h.pp_anim_pause})
set(h.pp_anim_pause,'Callback',{@callback_pause});
set(h.pp_anim_stop,'Callback',{@callback_stop,h.pp_anim_pause});


set(h.pp_anim_line_play,'Callback',{@callback_resume,h.pp_anim_line_axes,h.pp_anim_line_pause})
set(h.pp_anim_line_pause,'Callback',{@callback_pause});
set(h.pp_anim_line_stop,'Callback',{@callback_stop,h.pp_anim_line_pause});
set(h.pp_anim_line_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_anim_line_x_cb,'callback',...
    {@callback_pp_cb_x,h,h.pp_anim_line_x_edit,h.pp_anim_line_y_cb});
set(h.pp_anim_line_x_edit,'callback',{@callback_pp_edit_x,h});
set(h.pp_anim_line_y_cb,'callback',...
    {@callback_pp_cb_y,h,h.pp_anim_line_y_edit,h.pp_anim_line_x_cb});
set(h.pp_anim_line_y_edit,'callback',{@callback_pp_edit_y,h});
set(h.pp_anim_line_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_anim_line_setting_edit_zlim})
set(h.pp_anim_line_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_anim_line_setting_cb_spatlim,'callback',...
    {@callback_pp_setting_cb_spatlim,h.pp_anim_line_setting_edit_spatlim});
set(h.pp_anim_line_setting_edit_spatlim,'callback',{@callback_pp_setting_edit_spatlim,h});
set(h.pp_anim_line_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_anim_line_setting_edit_tlim})
set(h.pp_anim_line_setting_edit_tlim,'callback',...
    {@callback_pp_setting_edit_tlim,h});
set(h.pp_anim_line_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_anim_line_setting_edit_coarse});
set(h.pp_anim_line_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});

set(h.pp_anim_line_setting_cb_saveanim,'callback',...
    {@callback_pp_setting_cb_saveanim,h.pp_anim_line_setting_edit_gifset,h.pp_anim_line_setting_text_gifset})
set(h.pp_anim_line_setting_edit_gifset,'callback',{@callback_pp_setting_edit_gifset,h})

set(h.pp_anim_line_edit_bathyScale,'callback',{@callback_anim_edit_bathy_scale})
set(h.pp_anim_line_cb_bathy,'callback',{@callback_anim_cb_bathy,h.pp_anim_line_edit_bathyScale})


set(h.pp_anim_line_plot_button,'callback',{@callback_pp_plot_anim_line,h})


set(h.pp_validation_buoy_button_zoomin,'callback',{@callback_zoom_in,h})
set(h.pp_validation_buoy_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.pp_validation_buoy_button_datacursor,'callback',{@callback_datacursor,h})
set(h.pp_validation_buoy_button_pan,'callback',{@callback_pan,h})

set(h.pp_validation_buoy_x_edit,'callback',{@callback_pp_validation_buoy_edit_xy,h});
set(h.pp_validation_buoy_y_edit,'callback',{@callback_pp_validation_buoy_edit_xy,h});
set(h.pp_validation_buoy_time_edit,'callback',{@callback_pp_setting_edit_tlim,h})
set(h.pp_validation_buoy_popup_var,'callback',{@callback_pp_density_signal_popup_var,h})
set(h.pp_validation_buoy_measdata_button,'callback',{@callback_pp_validation_load_meas_data,h})

set(h.pp_validation_buoy_cb_timeshift,'callback',{@callback_pp_validation_buoy_cb_timeshift,h})
set(h.pp_validation_buoy_cb_timeshift_dt,'callback',{@callback_pp_validation_buoy_cb_timeshift_dt,h})
set(h.pp_validation_buoy_edit_timeshift_dt,'callback',{@callback_pp_validation_buoy_edit_timeshift_dt,h})
set(h.pp_validation_buoy_cb_timeshift_def,'callback',{@callback_pp_validation_buoy_cb_timeshift_def,h})

set(h.pp_validation_buoy_cb_signal,'callback',...
    {@callback_pp_cb_signal,h.pp_validation_buoy_cb_spectrum,...
    h.pp_validation_buoy_popup_spectrum_var,h.pp_validation_buoy_setting_cb_spsmooth,...
    h.pp_validation_buoy_setting_edit_spsmooth})
set(h.pp_validation_buoy_cb_spectrum,'callback',{@callback_pp_cb_spectrum,...
    h.pp_validation_buoy_popup_spectrum_var,h.pp_validation_buoy_cb_signal,...
    h.pp_validation_buoy_setting_cb_spsmooth,h.pp_validation_buoy_setting_edit_spsmooth})

set(h.pp_validation_buoy_setting_cb_zlim,'callback',...
    {@callback_pp_setting_cb_zlim,h.pp_validation_buoy_setting_edit_zlim})
set(h.pp_validation_buoy_setting_edit_zlim,'callback',...
    {@callback_pp_setting_edit_zlim,h});
set(h.pp_validation_buoy_setting_cb_horznlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.pp_validation_buoy_setting_edit_horznlim});
set(h.pp_validation_buoy_setting_edit_horznlim,'callback',{@callback_pp_setting_edit_horzlim,h});
set(h.pp_validation_buoy_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.pp_validation_buoy_setting_edit_coarse});
set(h.pp_validation_buoy_setting_edit_coarse,'callback',...
    {@callback_pp_setting_edit_coarse,h});
set(h.pp_validation_buoy_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.pp_validation_buoy_setting_popup_savefig})
set(h.pp_validation_buoy_setting_cb_spsmooth,'callback',{@callback_pp_setting_cb_coarse,...
    h.pp_validation_buoy_setting_edit_spsmooth})
set(h.pp_validation_buoy_setting_edit_spsmooth,'callback',{@callback_edit_param,h})
set(h.pp_validation_buoy_plot_button,'callback',{@callback_pp_validation_buoy_plot,h})
set(h.pp_validation_buoy_number_popup,'callback',{@callback_pp_validation_buoy_plot_popup,h})

set(h.pp_validation_quant_buoy_x_edit,'callback',{@callback_pp_validation_buoy_edit_xy,h});
set(h.pp_validation_quant_buoy_y_edit,'callback',{@callback_pp_validation_buoy_edit_xy,h});
set(h.pp_validation_quant_buoy_time_edit,'callback',{@callback_pp_setting_edit_tlim,h})
set(h.pp_validation_quant_buoy_cb_timeshift,'callback',{@callback_pp_validation_quant_buoy_cb_timeshift,h})
set(h.pp_validation_quant_buoy_cb_timeshift_dt,'callback',{@callback_pp_validation_quant_buoy_cb_timeshift_dt,h})
set(h.pp_validation_quant_buoy_edit_timeshift_dt,'callback',{@callback_pp_validation_quant_buoy_edit_timeshift_dt,h})
set(h.pp_validation_quant_buoy_cb_timeshift_def,'callback',{@callback_pp_validation_quant_buoy_cb_timeshift_def,h})
set(h.pp_validation_quant_buoy_measdata_button,'callback',{@callback_pp_validation_quant_load_meas_data,h})
set(h.pp_validation_quant_buoy_button,'callback',{@callback_pp_validation_quant_buoy_calc,h});
 
set(h.IF_project_edit_name,'callback',{@callback_IF_projectname,h})
set(h.IF_project_button_projdir,'callback',{@callback_IF_project_button_projdir,h})
set(h.IF_project_load_data,'callback',{@callback_IF_loaddata,h});
set(h.IF_calc_time_edit_interval_init,'callback',{@callback_IF_edit,h})
set(h.IF_calc_time_edit_interval_end,'callback',{@callback_IF_edit,h})
set(h.IF_calc_time_edit_step,'callback',{@callback_IF_edit,h})

set(h.IF_calc_x_edit_interval_init,'callback',{@callback_IF_edit,h})
set(h.IF_calc_x_edit_interval_end,'callback',{@callback_IF_edit,h})
set(h.IF_calc_x_edit_step,'callback',{@callback_IF_edit,h})
set(h.IF_calc_y_edit_interval_init,'callback',{@callback_IF_edit,h})
set(h.IF_calc_y_edit_interval_end,'callback',{@callback_IF_edit,h})
set(h.IF_calc_y_edit_step,'callback',{@callback_IF_edit,h})
set(h.IF_calc_z_edit_interval_init,'callback',{@callback_IF_edit,h})
set(h.IF_calc_z_edit_interval_end,'callback',{@callback_IF_edit,h})
set(h.IF_calc_z_edit_step,'callback',{@callback_IF_edit,h})
set(h.IF_calc_z_edit_input,'callback',{@callback_IF_edit,h})
set(h.IF_cal_z_checkbox_equidistant,'callback',{@callback_IF_cb_Zequidistant,h});
set(h.IF_cal_z_checkbox_nonequidistant,'callback',{@callback_IF_cb_Znonequidistant,h});
set(h.IF_calc_IDdata_edit,'callback',{@callback_IF_edit,h})
set(h.IF_calc_pushbutton,'callback',{@callback_IF_pb_calculation,h})
set(h.IF_calc_checkbox,'callback',{@callback_IF_cb_calc,h})
set(h.IF_calc_load_checkbox,'callback',{@callback_IF_cb_calc_load,h})
set(h.IF_calc_load_push,'callback',{@callback_IF_pb_calc_load,h})

set(h.IF_plot_density_axesopt_popup,'callback',{@callback_IF_pop_density_axesopt,h})
set(h.IF_plot_density_var1_edit,'callback',{@callback_IF_edit_density,h});
set(h.IF_plot_density_var2_edit,'callback',{@callback_IF_edit_density,h});
set(h.IF_plot_density_profile_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.IF_plot_density_profile_setting_edit_view})
set(h.IF_plot_density_profile_setting_edit_view,'callback',{@callback_pp_setting_edit_view,h})
set(h.IF_plot_density_profile_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.IF_plot_density_profile_setting_edit_clim})
set(h.IF_plot_density_profile_setting_edit_clim,'callback',{@callback_pp_setting_edit_clim,h})
set(h.IF_plot_density_profile_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.IF_plot_density_profile_setting_edit_xlim})
set(h.IF_plot_density_profile_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h})
set(h.IF_plot_density_profile_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.IF_plot_density_profile_setting_edit_ylim})
set(h.IF_plot_density_profile_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h})
set(h.IF_plot_density_profile_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.IF_plot_density_profile_setting_edit_coarse})
set(h.IF_plot_density_profile_setting_edit_coarse,'callback',{@callback_pp_setting_edit_coarse,h})
set(h.IF_plot_density_profile_setting_cb_level,'callback',...
    {@callback_pp_setting_cb_level1,h,h.IF_plot_density_profile_setting_cb_level,h.IF_plot_density_profile_setting_edit_level})
set(h.IF_plot_density_profile_setting_edit_level,'callback',{@callback_pp_setting_edit_level,h})
set(h.IF_plot_density_profile_setting_cb_savefig,'callback',...
    {@callback_pp_setting_cb_savefig,h.IF_plot_density_profile_setting_popup_savefig})
set(h.IF_plot_density_button_zoomin,'callback',{@callback_zoom_in,h});
set(h.IF_plot_density_button_zoomout,'callback',{@callback_zoom_out,h})
set(h.IF_plot_density_button_datacursor,'callback',{@callback_datacursor,h})
set(h.IF_plot_density_button_rotate,'callback',{@callback_rotate,h})
set(h.IF_plot_density_button_pan,'callback',{@callback_pan,h})

set(h.IF_plot_density_profile_plot_button,'callback',{@callback_IF_button_density,h})

set(h.IF_anim_play,'Callback',{@callback_resume,h.IF_anim_density_axes,h.IF_anim_pause})
set(h.IF_anim_pause,'Callback',{@callback_pause});
set(h.IF_anim_stop,'Callback',{@callback_stop,h.IF_anim_pause});

set(h.IF_anim_density_axesopt_popup,'callback',{@callback_IF_pop_animdensity_axesopt,h})
set(h.IF_anim_density_var2_edit,'callback',{@callback_IF_edit_density,h});

set(h.IF_anim_density_setting_cb_view,'callback',...
    {@callback_pp_setting_cb_view,h.IF_anim_density_setting_edit_view})
set(h.IF_anim_density_setting_edit_view,'callback',{@callback_pp_setting_edit_view,h})
set(h.IF_anim_density_setting_cb_clim,'callback',...
    {@callback_pp_setting_cb_clim,h.IF_anim_density_setting_edit_clim})
set(h.IF_anim_density_setting_edit_clim,'callback',{@callback_pp_setting_edit_clim,h})
set(h.IF_anim_density_setting_cb_xlim,'callback',...
    {@callback_pp_setting_cb_xlim,h.IF_anim_density_setting_edit_xlim})
set(h.IF_anim_density_setting_edit_xlim,'callback',{@callback_pp_setting_edit_xlim,h})
set(h.IF_anim_density_setting_cb_ylim,'callback',...
    {@callback_pp_setting_cb_ylim,h.IF_anim_density_setting_edit_ylim})
set(h.IF_anim_density_setting_edit_ylim,'callback',{@callback_pp_setting_edit_ylim,h})
set(h.IF_anim_density_setting_cb_tlim,'callback',...
    {@callback_pp_setting_cb_tlim,h.IF_anim_density_setting_edit_tlim})
set(h.IF_anim_density_setting_edit_tlim,'callback',{@callback_pp_setting_edit_tlim,h})
set(h.IF_anim_density_setting_cb_coarse,'callback',...
    {@callback_pp_setting_cb_coarse,h.IF_anim_density_setting_edit_coarse})
set(h.IF_anim_density_setting_edit_coarse,'callback',{@callback_pp_setting_edit_coarse,h})
set(h.IF_anim_density_setting_cb_level,'callback',...
    {@callback_pp_setting_cb_level1,h,h.IF_anim_density_setting_cb_level,h.IF_anim_density_setting_edit_level})
set(h.IF_anim_density_setting_edit_level,'callback',{@callback_pp_setting_edit_level,h})
set(h.IF_anim_density_setting_cb_saveanim,'callback',...
    {@callback_pp_setting_cb_saveanim,h.IF_anim_density_setting_edit_gifset,h.IF_anim_density_setting_text_gifset})
set(h.IF_anim_density_setting_edit_gifset,'callback',{@callback_pp_setting_edit_gifset,h})
set(h.IF_anim_density_plot_button,'callback',{@calback_IF_anim_plot,h})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Set initialization for gui handles%%%%%%%%%%%%%%%%%%%%%
set(h.monitorbox,'String','>>')
initialization_gui(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function guiAB2dwaveResizeFcn(~,~,h)
        walltabsize = getpixelposition(h.wall_table);
        ts=walltabsize(3)/12.4;
        set(h.wall_table,'ColumnWidth',{ts ts ts ts ts ts ts ts ts*1.2 ts ts ts});
        
        damptabsize = getpixelposition(h.damping_table);
        dts=damptabsize(3)/6.2;
        set(h.damping_table,'ColumnWidth',{dts dts dts dts dts dts});
        
        inflptabsize = getpixelposition(h.waveinput_influx_proptable);
        infp=inflptabsize(3)/7.7;
        set(h.waveinput_influx_proptable,'ColumnWidth',{infp*1.5 infp infp infp infp infp infp})
        
        tsize = getpixelposition(h.waveinput_influx_methodtable);
        tfs=tsize(3)/14.2;
        set(h.waveinput_influx_methodtable,'ColumnWidth',{tfs tfs tfs tfs tfs tfs tfs tfs tfs tfs tfs tfs tfs tfs})
        
        frictsize=getpixelposition(h.bathymetry_table_friction);
        frtss=frictsize(3)/6.2;
        set(h.bathymetry_table_friction, 'ColumnWidth',{frtss frtss frtss frtss frtss frtss})
        
        quant_buoytablesize = getpixelposition(h.pp_quant_buoy_table);
        infpts=quant_buoytablesize(3)/9.3;
        set(h.pp_quant_buoy_table,'ColumnWidth',{infpts infpts infpts infpts infpts infpts infpts infpts infpts})    
   
        val_quant_buoytablesize = getpixelposition(h.pp_validation_quant_buoy_table);
        infpts=val_quant_buoytablesize(3)/11;
        set(h.pp_validation_quant_buoy_table,'ColumnWidth',{infpts infpts infpts infpts infpts infpts infpts infpts infpts infpts})    
     
    end

    function tableDelSelection(hObj,eventdata,h,href)
        index = eventdata.Indices;
        
        if any(index)
            inf_del.Id=1;
            inf_del.row= index(1);
            inf_del.col= index(2);
        else
            inf_del.Id=0;
        end
        set(href,'userdata',inf_del)
    end

    function callback_menu_file_newproj(hObj,eventdata,h)
        ID=close(h.fig);
        if ID==1
        HAWASSI_AB_startpage;
        end
    end

    function callback_menu_file_saveproj(hObj,eventdata,h)
        try
            GUIinput=fun_passing_handles_preproc(h);
            func_save_guistate(GUIinput);
            msgbox(sprintf('Project %s has been saved',GUIinput.proj.name),'Done');
        catch
        end
    end

    function callback_menu_file_clearproj(hObj,eventdata,h)
        reset_handles_gui(h);
        reset_handles_postproc(h);
    end

    function callback_calculator(hObj,eventdata)
        Calculator;
    end

    function callback_about_AB2(hObj,eventdata)
       About_AB2d; 
    end

    function callback_doc(hObj,eventdata)
        if ~isdeployed
            winopen('\Toolbox\lib2d\Tools\misc\Manual_AB2_v1_1.pdf')
        else
            winopen('Manual_AB2_v1_1.pdf')
        end
    end

    function callback_menu_file_quitproj(hObj,eventdata,h)
        projname=get(h.project_edit_name,'string');
        user_response_quit = questdlg(sprintf('Do you want to save changes you made to %s ?',projname),'AB2D-Wave','Save','Don''t Save','Cancel', 'Cancel');
        switch user_response_quit
            case 'Save'
                callback_menu_file_saveproj(hObj,eventdata,h);
                evalin('base','clear all; clc;');
                delete(hObj);
                delete(gcf);
                %delete(findall(0, 'type', 'figure'));
            case 'Don''t Save'
                evalin('base','clear all; clc;');
                delete(hObj);
                delete(gcf)
                %delete(findall(0, 'type', 'figure'));
            case 'Cancel'
        end
    end


    function callback_project_button_projdir(hObj,eventdata,h)
        pathnow=h.projectdirectory;
        
        workingdir = uigetdir(pathnow,'browse a directory');
        if workingdir ~= 0
            set(h.project_popup_projdir,'String',workingdir);
        else
            set(h.monitorbox,'String',['>>Warning: No directory is loaded'],'foregroundcolor','k');
            uicontrol(h.project_button_projdir)
            return;
        end
        
    end

    function callback_projectname(hObj,eventdata,h)
        projname=get(hObj,'String');
        projdir=get(h.project_popup_projdir,'string');
        
        if isempty(projdir)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Specify a working directory');
            uicontrol(h.project_edit_name);
            return;
        end
        
        if exist([projdir,'/',projname],'dir')
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Warning: Project exists already, it will be overwritten');
            uicontrol(h.project_edit_name);
            flag_warn_workdir=1;
        else
            flag_warn_workdir=0;
            set(h.monitorbox,'foregroundcolor','k','string','>>');
        end
        set(h.fig,'Name',['AB2D-Wave [',projname,']']);
    end

    function callback_pp_customize_data(hObj,eventdata,h)
        
        dat=get(h.pp_project_load_data,'userdata');
        if isempty(dat)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> load simulation data');
            return;
        else
            path=get(h.pp_project_popup_projdir,'string');
            dir=get(h.pp_project_edit_name,'string');
            savedir=[path,'\',dir,'\'];
            if ~isdir(savedir)
                mkdir(savedir)
            end
        end
        Gui_customize_data(dat,savedir);
    end



    function callback_modeldyn(hObj,eventdata,h)
        contents = cellstr(get(hObj,'String'));
        dynmodel=contents{get(hObj,'Value')};
        set(h.monitorbox,'String','>>');
        ID=get(hObj,'Value');
        
        if h.fullversion==0
            if ID==2||ID==3
                set(h.model_popup_dynamic,'value',1);
                set(h.monitorbox,'String',['>> HS2 or HS3 is not available'...
                    ' in the demo version'],'foregroundcolor','red');
            else
                set(h.monitorbox,'String','>>','foregroundcolor','k');
            end
        end
        if ID==1
            set(h.model_popup_breaking,'value',2)
            callback_breaking(h.model_popup_breaking,[],h);
            set(h.cutfrac_k,'string','2','userdata',2);
        elseif ID==2
            set(h.cutfrac_k,'string','4','userdata',4);
        elseif ID==3
            set(h.cutfrac_k,'string','6','userdata',6);
        end
        
        if get(h.waveinput_influx_popup,'value')==2
            if ID>1
                set(h.waveinput_influx_checkbox_nonlinadj,'visible','on',...
                    'enable','on')
                if get(h.waveinput_influx_checkbox_nonlinadj,'value')==1
                    set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','on')
                    set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','on');
                else
                    set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','off')
                    set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','off');
                end
            else
                set(h.waveinput_influx_checkbox_nonlinadj,'visible','on',...
                    'enable','off')
                set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','off')
                set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','off');
            end
        end
        
    end

    function callback_dispersionmenu(hObj,eventdata,h)
       
    end

    function callback_breaking(hObj,eventdata,h)
        val = get(hObj,'Value');
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        if val == 2
            set(h.model_break_text4, 'Visible', 'off');
            set(h.model_break_text5, 'Visible', 'off');
            set(h.model_break_edit_initiation, 'Visible', 'off');
            set(h.model_break_text6, 'Visible', 'off');
            set(h.model_break_text7, 'Visible', 'off');
            set(h.model_break_edit_termination, 'Visible', 'off');
            set(h.model_break_text8, 'Visible', 'off');
            set(h.model_break_edit_Tstar, 'Visible', 'off');
        else
            
            dyntype=get(h.model_popup_dynamic,'value');
            if dyntype==1
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Use only a nonlinear model for breaking wave simulation!')
                uicontrol(h.model_popup_dynamic);
                
                set(h.model_break_text4, 'Visible', 'off');
                set(h.model_break_text5, 'Visible', 'off');
                set(h.model_break_edit_initiation, 'Visible', 'off');
                set(h.model_break_text6, 'Visible', 'off');
                set(h.model_break_text7, 'Visible', 'off');
                set(h.model_break_edit_termination, 'Visible', 'off');
                set(h.model_break_text8, 'Visible', 'off');
                set(h.model_break_edit_Tstar, 'Visible', 'off');
                set(h.model_popup_breaking,'value',2);
                
            else
                if h.fullversion==1
                    set(h.model_break_text4, 'Visible', 'on');
                    set(h.model_break_text5, 'Visible', 'on');
                    set(h.model_break_edit_initiation, 'Visible', 'on');
                    set(h.model_break_text6, 'Visible', 'on');
                    set(h.model_break_text7, 'Visible', 'on');
                    set(h.model_break_edit_termination, 'Visible', 'on');
                    set(h.model_break_text8, 'Visible', 'on');
                    set(h.model_break_edit_Tstar, 'Visible', 'on');
                else
                    set(h.model_popup_breaking,'value',2);
                    set(h.monitorbox,'String',['>> Breaking is not available'...
                        ' for demo version'],'foregroundcolor','red');
                end
            end
        end
    end

    function callback_breaking_initiation(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2double(get(hObj,'String'));
        if param<0 || param>1
            set(h.monitorbox,'String','>>the value must be in (0,1)','foregroundcolor','r');
            uicontrol(h.model_break_edit_initiation);
        else
            set(h.model_break_edit_initiation,'Userdata',param);
        end
    end

    function callback_breaking_termination(hObj, eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2double(get(hObj,'String'));
        if param<0 || param>0.5
            set(h.monitorbox,'String','>>the value must be in (0,0.5)','foregroundcolor','r');
            uicontrol(h.model_break_edit_termination);
        else
            set(h.model_break_edit_termination,'Userdata',param);
        end
    end
    function callback_breaking_Tstar(hObj, eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2double(get(hObj,'String'));
        if param<0 || param>0.5
            set(h.monitorbox,'String','>>the value must be in (0.1,0.5)','foregroundcolor','r');
            uicontrol(h.model_break_edit_Tstar);
        else
            set(h.model_break_edit_Tstar,'Userdata',param);
        end
    end

    function callback_current_cb(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.model_current_text1,'enable','on');
            set(h.model_current_text2,'enable','on');
            set(h.model_current_text3,'enable','on');
            set(h.model_current_ux_edit,'enable','on');
            set(h.model_current_uy_edit,'enable','on');
        else
            set(h.model_current_text1,'enable','off');
            set(h.model_current_text2,'enable','off');
            set(h.model_current_text3,'enable','off');
            set(h.model_current_ux_edit,'enable','off');
            set(h.model_current_uy_edit,'enable','off');
        end
        
    end

    function callback_current_ux(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param);
    end

    function callback_current_uy(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param);
    end

    function callback_spatial_xmin(hObj, eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        parammax=get(h.spatial_edit_xmax,'userdata');
        if ~isempty(parammax)
            if param>=parammax
                set(h.monitorbox,'String','>>wrong input: xmin>=xmax','foregroundcolor','r');
                uicontrol(h.spatial_edit_xmin);
            end
        end
        set(h.spatial_edit_xmin,'Userdata',param);
    end
    function callback_spatial_xmax(hObj, eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        parammin=get(h.spatial_edit_xmin,'userdata');
        if ~isempty(parammin)
            if param<=parammin
                set(h.monitorbox,'String','>>wrong input: xmax<=xmin','foregroundcolor','r');
                uicontrol(h.spatial_edit_xmax);
            end
        end
        set(h.spatial_edit_xmax,'Userdata',param);
    end
    function callback_spatial_ymin(hObj, eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        parammax=get(h.spatial_edit_ymax,'userdata');
        if ~isempty(parammax)
            if param>=parammax
                set(h.monitorbox,'String','>>wrong input: ymin>=ymax','foregroundcolor','r');
                uicontrol(h.spatial_edit_ymin);
            end
        end
        set(h.spatial_edit_ymin,'Userdata',param);
    end
    function callback_spatial_ymax(hObj, eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        parammin=get(h.spatial_edit_ymin,'userdata');
        if ~isempty(parammin)
            if param<=parammin
                set(h.monitorbox,'String','>>wrong input: ymax<=ymin','foregroundcolor','r');
                uicontrol(h.spatial_edit_ymax);
            end
        end
        set(h.spatial_edit_ymax,'Userdata',param);
    end

    function callback_spatial_dx(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.spatial_edit_dx,'Userdata',param);
    end

    function callback_spatial_dy(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.spatial_edit_dy,'Userdata',param);
    end

    function callback_cutfrac_k(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.cutfrac_k,'Userdata',param);
    end



    function callback_wall_popup(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        str = get(hObj, 'String');
        val = get(hObj,'Value');
        
        if val==1
            if strcmpi(h.moduleRestrict,'Flat_NoWall')
                    set(hObj,'value',2);val=2;
                    set(h.monitorbox,'String','>> The wall is not available for this licence. ','foregroundcolor','k');
            end
        end
        
        if val == 2
            set(h.wall_button_addrow,'visible','off');
            set(h.wall_button_deleterow,'visible','off');
            set(h.wall_table,'visible','off');
        else
            set(h.wall_button_addrow,'visible','on');
            set(h.wall_button_deleterow,'visible','on');
            set(h.wall_table,'visible','on');
        end
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end


    function callback_bdyassim_popup(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        str = get(hObj, 'String');
        val = get(hObj,'Value');
        
        if strcmp(str{val},'No') == 1
            set(h.waveinput_bdy_assim_text2,'visible','off');
            set(h.waveinput_bdy_assim_shape_popup,'visible','off');
            set(h.waveinput_bdy_assim_text3,'visible','off');
            set(h.waveinput_bdy_assim_text4,'visible','off');
            set(h.waveinput_bdy_assim_edit_R1,'visible','off');
            set(h.waveinput_bdy_assim_text3a,'visible','off');
            set(h.waveinput_bdy_assim_edit_smooth,'visible','off');
            
            set(h.waveinput_bdy_assim_edit_xc,'visible','off');
            set(h.waveinput_bdy_assim_text4a,'visible','off');
            set(h.waveinput_bdy_assim_edit_yc,'visible','off');
            set(h.waveinput_bdy_assim_text5a,'visible','off');
            
            set(h.waveinput_bdy_assim_text7,'visible','off');
            set(h.waveinput_bdy_assim_edit_data,'visible','off');
            set(h.waveinput_bdy_assim_button_load,'visible','off');
            set(h.waveinput_bdy_assim_text8,'visible','off');
            set(h.waveinput_bdy_assim_propdir_popup,'visible','off');
            set(h.waveinput_bdy_assim_time_text1,'visible','off');
            if ID_model_phiform==1
            set(h.waveinput_bdy_assim_button_load_phi,'visible','off');
            set(h.waveinput_bdy_assim_edit_data_phi,'visible','off');
            set(h.waveinput_bdy_assim_cb_phi,'visible','off');
            else
            set(h.waveinput_bdy_assim_button_load_v,'visible','off');    
            set(h.waveinput_bdy_assim_button_load_u,'visible','off');
            set(h.waveinput_bdy_assim_edit_data_vel,'visible','off');
            set(h.waveinput_bdy_assim_cb_vel,'visible','off');    
                
            end
            set(h.waveinput_bdy_assim_cb_nonlinAdj,'visible','off');
            set(h.waveinput_bdy_assim_text5,'visible','off');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'visible','off');
            set(h.waveinput_bdy_assim_text6,'visible','off');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'visible','off')
            
            
            set(h.waveinput_bdy_assim_time_edit_interval_init,'visible','off');
            set(h.waveinput_bdy_assim_time_text2,'visible','off');
            set(h.waveinput_bdy_assim_time_edit_interval_end,'visible','off');
            set(h.waveinput_bdy_assim_time_text3,'visible','off');
            set(h.waveinput_bdy_assim_time_text4,'visible','off');
            set(h.waveinput_bdy_assim_time_edit_step,'visible','off');
            set(h.waveinput_bdy_assim_time_text5,'visible','off');
            
            
        else
            set(h.waveinput_bdy_assim_text2,'visible','on');
            set(h.waveinput_bdy_assim_shape_popup,'visible','on');
            set(h.waveinput_bdy_assim_text3,'visible','on');
            set(h.waveinput_bdy_assim_text4,'visible','on');
            set(h.waveinput_bdy_assim_edit_R1,'visible','on');
            set(h.waveinput_bdy_assim_edit_smooth,'visible','on');
            set(h.waveinput_bdy_assim_text3a,'visible','on');
            set(h.waveinput_bdy_assim_edit_xc,'visible','on');
            set(h.waveinput_bdy_assim_text4a,'visible','on');
            set(h.waveinput_bdy_assim_edit_yc,'visible','on');
            set(h.waveinput_bdy_assim_text5a,'visible','on');
            
            set(h.waveinput_bdy_assim_text7,'visible','on');
            set(h.waveinput_bdy_assim_edit_data,'visible','on');
            set(h.waveinput_bdy_assim_button_load,'visible','on');
            set(h.waveinput_bdy_assim_text8,'visible','on');
            set(h.waveinput_bdy_assim_propdir_popup,'visible','on');
            set(h.waveinput_bdy_assim_time_text1,'visible','on');
            if ID_model_phiform==1
            set(h.waveinput_bdy_assim_button_load_phi,'visible','on');
            set(h.waveinput_bdy_assim_edit_data_phi,'visible','on');
            set(h.waveinput_bdy_assim_cb_phi,'visible','on');
            else
            set(h.waveinput_bdy_assim_button_load_v,'visible','on');    
            set(h.waveinput_bdy_assim_button_load_u,'visible','on');
            set(h.waveinput_bdy_assim_edit_data_vel,'visible','on');
            set(h.waveinput_bdy_assim_cb_vel,'visible','on');    
            end
            set(h.waveinput_bdy_assim_cb_nonlinAdj,'visible','on');
            set(h.waveinput_bdy_assim_text5,'visible','on');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'visible','on');
            set(h.waveinput_bdy_assim_text6,'visible','on');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'visible','on')
            
            
            set(h.waveinput_bdy_assim_time_edit_interval_init,'visible','on');
            set(h.waveinput_bdy_assim_time_text2,'visible','on');
            set(h.waveinput_bdy_assim_time_edit_interval_end,'visible','on');
            set(h.waveinput_bdy_assim_time_text3,'visible','on');
            set(h.waveinput_bdy_assim_time_text4,'visible','on');
            set(h.waveinput_bdy_assim_time_edit_step,'visible','on');
            set(h.waveinput_bdy_assim_time_text5,'visible','on');
        end
    end

    function callback_bdyassim_shape(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.waveinput_bdy_assim_text3,'enable','on')
            set(h.waveinput_bdy_assim_text4,'enable','on')
            set(h.waveinput_bdy_assim_edit_R1,'enable','on')
            set(h.waveinput_bdy_assim_text5,'enable','on')
            set(h.waveinput_bdy_assim_text3a,'enable','on')
            set(h.waveinput_bdy_assim_edit_xc,'enable','on')
            set(h.waveinput_bdy_assim_text4a,'enable','on')
            set(h.waveinput_bdy_assim_edit_yc,'enable','on')
            set(h.waveinput_bdy_assim_text5a,'enable','on')
             set(h.monitorbox,'String','>>','foregroundcolor', 'k');
        else
            set(h.monitorbox,'String','>>the user-defined bdy. shape is not available!','foregroundcolor', 'r');
            set(hObj,'value',1)
%            set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
%            [file_name,directory]=uigetfile([h.pathnow,'\','*.txt; *.dat; *.mat; *.asc'],'Load a boundary domain');
%             if directory~=0
%                 temp=load([directory,file_name]);
%                 if isstruct(temp)
%                     namevar = fieldnames(temp);
%                     my_data=temp.(namevar{1});
%                 else
%                     my_data=temp;
%                 end
%                 clearvars temp;
%                 set(hObj,'userdata',my_data)
%                 set(h.waveinput_bdy_assim_text3,'enable','off')
%                 set(h.waveinput_bdy_assim_text4,'enable','off')
%                 set(h.waveinput_bdy_assim_edit_R1,'enable','off')
%                 set(h.waveinput_bdy_assim_text5,'enable','off')
%                 set(h.waveinput_bdy_assim_text3a,'enable','off')
%                 set(h.waveinput_bdy_assim_edit_xc,'enable','off')
%                 set(h.waveinput_bdy_assim_text4a,'enable','off')
%                 set(h.waveinput_bdy_assim_edit_yc,'enable','off')
%                 set(h.waveinput_bdy_assim_text5a,'enable','off')
%                 set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
%                 
%             else
%                 set(h.waveinput_bdy_assim_text3,'enable','on')
%                 set(h.waveinput_bdy_assim_text4,'enable','on')
%                 set(h.waveinput_bdy_assim_edit_R1,'enable','on')
%                 set(h.waveinput_bdy_assim_text5,'enable','on')
%                 
%                 set(h.waveinput_bdy_assim_text3a,'enable','on')
%                 set(h.waveinput_bdy_assim_edit_xc,'enable','on')
%                 set(h.waveinput_bdy_assim_text4a,'enable','on')
%                 set(h.waveinput_bdy_assim_edit_yc,'enable','on')
%                 set(h.waveinput_bdy_assim_text5a,'enable','on')
%                 
%                 set(hObj,'value',1)
%                 set(h.monitorbox,'String','>>  No loaded data','foregroundcolor', 'k');
%             end
            
        end
    end

    function callback_bdyaasim_loaddata(hObj,eventdata,h,href1)
        set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
        [file_name,directory]=uigetfile([h.projectdirectory,'\','*.txt; *.dat; *.mat; *.asc'],'Load assimilation data');
        if directory~=0
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            FlagDat=1;
            if ~isstruct(my_data)
                Nx=my_data(2,3);Ny=my_data(3,3);Nt=my_data(1,3);
                if length(my_data(4,:))~=Nx
                    FlagDat=0;
                end
                if length(my_data(4:end,:))~=Ny*Nt
                    FlagDat=0;
                end
            else
                FlagDat=0;
            end
            
            if FlagDat==1
                set(hObj,'userdata',my_data);
                set(href1,'string',[file_name])
                set(h.waveinput_bdy_assim_time_edit_interval_init,'userdata',my_data(1,1),...
                    'string',num2str(my_data(1,1)));
                set(h.waveinput_bdy_assim_time_edit_interval_end,'userdata',my_data(1,2),...
                    'string',num2str(my_data(1,2)));
                dt=(my_data(1,2)-my_data(1,1))/(my_data(1,3)-1);
                set(h.waveinput_bdy_assim_time_edit_step,'userdata',dt,...
                    'string',num2str(dt));
                set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
            else
                set(h.monitorbox,'String','>> Wrong input file','foregroundcolor', 'r');
                set(href1,'string','')
                set(hObj,'userdata',[]);
            end
        else
            set(h.monitorbox,'String','>> No loaded data','foregroundcolor', 'k');
            set(href1,'string','')
            set(hObj,'userdata',[]);
        end
    end

    function callback_bdyaasim_cb_phi(hObj,eventdata,h)
        Id =get(hObj,'value');
        if Id==1
            set(h.waveinput_bdy_assim_button_load_phi,'enable','on')
            set(h.waveinput_bdy_assim_propdir_popup,'enable','off');
        else
            set(h.waveinput_bdy_assim_button_load_phi,'enable','off')
            set(h.waveinput_bdy_assim_propdir_popup,'enable','on');
        end
        
    end

   function callback_bdyaasim_cb_vel(hObj,eventdata,h)
        Id =get(hObj,'value');
        if Id==1
            set(h.waveinput_bdy_assim_button_load_u,'enable','on')
            set(h.waveinput_bdy_assim_button_load_v,'enable','on')
            set(h.waveinput_bdy_assim_propdir_popup,'enable','off');
        else
            set(h.waveinput_bdy_assim_button_load_u,'enable','off')
            set(h.waveinput_bdy_assim_button_load_v,'enable','off')
            set(h.waveinput_bdy_assim_propdir_popup,'enable','on');
        end
        
   end

    function callback_bdyaasim_loaddata_vel(hObj,eventdata,h,href1)
        set(h.monitorbox,'String','>>loading data','foregroundcolor', 'k');
        [file_name,directory]=uigetfile([h.projectdirectory,'\','*.txt; *.dat; *.mat; *.asc'],'Load assimilation velocity data');
        if directory~=0
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            FlagDat=1;
            if ~isstruct(my_data)
                Nx=my_data(2,3);Ny=my_data(3,3);Nt=my_data(1,3);
                if length(my_data(4,:))~=Nx
                    FlagDat=0;
                end
                if length(my_data(4:end,:))~=Ny*Nt
                    FlagDat=0;
                end
            else
                FlagDat=0;
            end
            
            if FlagDat==1
                set(hObj,'userdata',my_data);
                set(href1,'string',[file_name])
                set(h.waveinput_bdy_assim_time_edit_interval_init,'userdata',my_data(1,1),...
                    'string',num2str(my_data(1,1)));
                set(h.waveinput_bdy_assim_time_edit_interval_end,'userdata',my_data(1,2),...
                    'string',num2str(my_data(1,2)));
                dt=(my_data(1,2)-my_data(1,1))/(my_data(1,3)-1);
                set(h.waveinput_bdy_assim_time_edit_step,'userdata',dt,...
                    'string',num2str(dt));
                set(h.monitorbox,'String','>> data has been loaded','foregroundcolor', 'k');
            else
                set(h.monitorbox,'String','>> Wrong input file','foregroundcolor', 'r');
                set(href1,'string','')
                set(hObj,'userdata',[]);
            end
        else
            set(h.monitorbox,'String','>> No loaded data','foregroundcolor', 'k');
            set(href1,'string','')
            set(hObj,'userdata',[]);
        end
    end

    function callback_bdyassim_nonlinadj (hObj,eventdata,h)
        Id =get(hObj,'value');
        if Id==1
            set(h.waveinput_bdy_assim_text5,'enable','on')
            set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'enable','on');
            set(h.waveinput_bdy_assim_text6,'enable','on');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'enable','on')
            
        else
            set(h.waveinput_bdy_assim_text5,'enable','off')
            set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'enable','off');
            set(h.waveinput_bdy_assim_text6,'enable','off');
            set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'enable','off')
        end
        
    end

    function callback_addrow_table(hObj,eventdata,h,htable,defdata)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        oldData = get(htable,'Data');
        newData = [oldData; defdata];
        set(htable,'Data',newData);
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_deleterow_table(hObj,eventdata,h,htable)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Ind_del=get(hObj,'userdata');
        if ~isempty(Ind_del)
            if Ind_del.Id==1
                oldData = get(htable,'Data');
                mask = (1:size(oldData,1))';
                mask(Ind_del.row) = [];
                
                newData = oldData(mask ,:);
                set(htable,'Data',newData);
                
                guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            else
                set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
            end
        else
            set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
        end
   end


    function callback_damping_popup(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        str = get(hObj, 'String');
        val = get(hObj,'Value');
        
        if strcmp(str{val},'No') == 1
            set(h.damping_button_addrow,'visible','off')
            set(h.damping_button_deleterow,'visible','off')
            set(h.damping_table,'visible','off')
        else
            set(h.damping_button_addrow,'visible','on')
            set(h.damping_button_deleterow,'visible','on')
            set(h.damping_table,'visible','on')
        end
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_damping_addrow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        oldData = get(h.damping_table,'Data');
        newData = [oldData; {'-','-','-','-','-','-'}];
        set(h.damping_table,'Data',newData);
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_damping_deleterow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Ind_del=get(hObj,'userdata');
        if ~isempty(Ind_del)
            if Ind_del.Id==1
                oldData = get(h.damping_table,'Data');
                if Ind_del.row~=1
                mask = (1:size(oldData,1))';
                mask(Ind_del.row) = [];
                newData = oldData(mask,:);
                set(h.damping_table,'Data',newData);
                else
                set(h.monitorbox,'String','>>The first row in the table cannot be deleted','foregroundcolor','r');
                end
                guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            else
                set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
            end
        else
            set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
        end
        
    end

    function callback_bathy_popup(hObj,eventdata,h)
        str = get(hObj, 'String');
        val = get(hObj,'Value');
%         if val==4 || val==5
%            set(h.monitorbox,'String','>> Shore is not available','foregroundcolor','k'); 
%            val=1; 
%            set(hObj,'Value',1);
%         end
        if strcmpi(h.moduleRestrict,'Flat_NoWall')
            if val>1
             set(hObj,'value',1);val=1;
             set(h.monitorbox,'String','>> Only flat bottom is available for this licence. ','foregroundcolor','k');    
            end
        end
        
        switch str{val}
            case 'Flat'
                set(h.bathymetry_text1, 'Visible', 'on');
                set(h.bathymetry_edit_depth, 'Visible', 'on');
                set(h.bathymetry_text2, 'Visible', 'on');
                set(h.bathymetry_text3,'visible','off');
                set(h.bathymetry_edit_mindepth,'visible','off');
                set(h.bathymetry_text4,'visible','off');
                set(h.bathymetry_text5,'visible','off');
                set(h.bathymetry_edit_maxdepth,'visible','off');
                set(h.bathymetry_text6,'visible','off');
                set(h.bathymetry_text7,'visible','off');
                set(h.bathymetry_edit_slope,'visible','off');
                set(h.bathymetry_text8,'visible','off');
                set(h.bathymetry_text9,'visible','off');
                set(h.bathymetry_edit_startslope,'visible','off');
                set(h.bathymetry_text10,'visible','off');
                 set(h.bathymetry_text11a,'visible','off');
                set(h.bathymetry_edit_mindepthshore,'visible','off');
                set(h.bathymetry_text11b,'visible','off');

                set(h.bathymetry_text11,'visible','off');
                set(h.bathymetry_edit_maxdepthshore,'visible','off');
                set(h.bathymetry_text12,'visible','off');
                set(h.bathymetry_text13,'visible','off');
                set(h.bathymetry_edit_slopeshore,'visible','off');
                set(h.bathymetry_text14,'visible','off');
                set(h.bathymetry_text15,'visible','off');
                set(h.bathymetry_edit_shoreposition,'visible','off');
                set(h.bathymetry_text16,'visible','off');
                set(h.bathymetry_text17,'visible','off');
                set(h.bathymetry_edit_filename,'visible','off');
                set(h.bathymetry_button_load,'visible','off');
                set(h.bathymetry_text18,'visible','off');
                set(h.bathymetry_popup_interpdepth,'visible','off');
                set(h.bathymetry_text19,'visible','off'); % nunu
                set(h.bathymetry_edit_middepth,'visible','off'); % nunu
                set(h.bathymetry_text20,'visible','off'); % nunu
    
            case 'Slope (in x-axis)'
                set(h.bathymetry_edit_mindepth,'tooltipString','Specify depth in xmin');
                set(h.bathymetry_edit_maxdepth,'tooltipString','Specify depth in xmax');
                set(h.bathymetry_edit_startslope,'tooltipString','Specify the x-position of the foot of slope ');
                set(h.bathymetry_edit_middepth,'tooltipString','Specify depth-reference for 3-depth interpolation'); % nunu

                set(h.bathymetry_text1, 'Visible', 'off');
                set(h.bathymetry_edit_depth, 'Visible', 'off');
                set(h.bathymetry_text2, 'Visible', 'off');
                set(h.bathymetry_text3,'visible','on','string','Depth (xmin):');
                set(h.bathymetry_edit_mindepth,'visible','on');
                set(h.bathymetry_text4,'visible','on');
                set(h.bathymetry_text5,'visible','on','string','(xmax):');
                set(h.bathymetry_edit_maxdepth,'visible','on');
                set(h.bathymetry_text6,'visible','on');
                set(h.bathymetry_text7,'visible','on');
                set(h.bathymetry_edit_slope,'visible','on');
                set(h.bathymetry_text8,'visible','on');
                set(h.bathymetry_text9,'visible','off');
                set(h.bathymetry_edit_startslope,'visible','on');
                set(h.bathymetry_text10,'visible','on');
                set(h.bathymetry_text11a,'visible','off');
                set(h.bathymetry_edit_mindepthshore,'visible','off');
                set(h.bathymetry_text11b,'visible','off');

                set(h.bathymetry_text11,'visible','off');
                set(h.bathymetry_edit_maxdepthshore,'visible','off');
                set(h.bathymetry_text12,'visible','off');
                set(h.bathymetry_text13,'visible','off');
                set(h.bathymetry_edit_slopeshore,'visible','off');
                set(h.bathymetry_text14,'visible','off');
                set(h.bathymetry_text15,'visible','off');
                set(h.bathymetry_edit_shoreposition,'visible','off');
                set(h.bathymetry_text16,'visible','off');
                set(h.bathymetry_text17,'visible','off');
                set(h.bathymetry_edit_filename,'visible','off');
                set(h.bathymetry_button_load,'visible','off');
                set(h.bathymetry_text18,'visible','on');
                set(h.bathymetry_popup_interpdepth,'visible','on');
               	set(h.bathymetry_text19,'visible','on'); % nunu
                set(h.bathymetry_edit_middepth,'visible','on'); % nunu
                set(h.bathymetry_text20,'visible','on'); % nunu

            case 'Slope (in y-axis)'
                set(h.bathymetry_edit_mindepth,'tooltipString','Specify depth in ymin');
                set(h.bathymetry_edit_maxdepth,'tooltipString','Specify depth in ymax');
                set(h.bathymetry_edit_startslope,'tooltipString','Specify the y-position of the foot of slope ')
		        set(h.bathymetry_edit_middepth,'tooltipString','Specify depth-reference for 3-depth interpolation'); % nunu

                set(h.bathymetry_text1, 'Visible', 'off');
                set(h.bathymetry_edit_depth, 'Visible', 'off');
                set(h.bathymetry_text2, 'Visible', 'off');
                set(h.bathymetry_text3,'visible','on','string','Depth (ymin):');
                set(h.bathymetry_edit_mindepth,'visible','on');
                set(h.bathymetry_text4,'visible','on');
                set(h.bathymetry_text5,'visible','on','string','(ymax):');
                set(h.bathymetry_edit_maxdepth,'visible','on');
                set(h.bathymetry_text6,'visible','on');
                set(h.bathymetry_text7,'visible','on');
                set(h.bathymetry_edit_slope,'visible','on');
                set(h.bathymetry_text8,'visible','off');
                set(h.bathymetry_text9,'visible','on');
                set(h.bathymetry_edit_startslope,'visible','on');
                set(h.bathymetry_text10,'visible','on');
                set(h.bathymetry_text11a,'visible','off');
                set(h.bathymetry_edit_mindepthshore,'visible','off');
                set(h.bathymetry_text11b,'visible','off');

                set(h.bathymetry_text11,'visible','off');
                set(h.bathymetry_edit_maxdepthshore,'visible','off');
                set(h.bathymetry_text12,'visible','off');
                set(h.bathymetry_text13,'visible','off');
                set(h.bathymetry_edit_slopeshore,'visible','off');
                set(h.bathymetry_text14,'visible','off');
                set(h.bathymetry_text15,'visible','off');
                set(h.bathymetry_edit_shoreposition,'visible','off');
                set(h.bathymetry_text16,'visible','off');
                set(h.bathymetry_text17,'visible','off');
                set(h.bathymetry_edit_filename,'visible','off');
                set(h.bathymetry_button_load,'visible','off');
                set(h.bathymetry_text18,'visible','on');
                set(h.bathymetry_popup_interpdepth,'visible','on');
                set(h.bathymetry_text19,'visible','on'); % nunu
                set(h.bathymetry_edit_middepth,'visible','on'); % nunu
                set(h.bathymetry_text20,'visible','on'); % nunu

            case 'Shore (in x-axis)'
                set(h.bathymetry_edit_slopeshore,'tooltipString','Specify the slope gradient ')
		        set(h.bathymetry_edit_shoreposition,'tooltipString','Specify a x shore position');
                set(h.bathymetry_edit_middepth,'tooltipString','Specify depth-reference for 3-depth interpolation'); % nunu
                set(h.bathymetry_edit_maxdepthshore,'tooltipString','Specify a maximum depth');
                set(h.bathymetry_edit_mindepthshore,'tooltipString','Specify a minimum depth (i.e 2% Hs)');
                
                set(h.bathymetry_text1, 'Visible', 'off');
                set(h.bathymetry_edit_depth, 'Visible', 'off');
                set(h.bathymetry_text2, 'Visible', 'off');
                set(h.bathymetry_text3,'visible','off');
                set(h.bathymetry_edit_mindepth,'visible','off');
                set(h.bathymetry_text4,'visible','off');
                set(h.bathymetry_text5,'visible','off');
                set(h.bathymetry_edit_maxdepth,'visible','off');
                set(h.bathymetry_text6,'visible','off');
                set(h.bathymetry_text7,'visible','off');
                set(h.bathymetry_edit_slope,'visible','off');
                set(h.bathymetry_text8,'visible','off');
                set(h.bathymetry_text9,'visible','off');
                set(h.bathymetry_edit_startslope,'visible','off');
                set(h.bathymetry_text10,'visible','off');
                set(h.bathymetry_text11,'visible','on');
                set(h.bathymetry_edit_maxdepthshore,'visible','on');
                set(h.bathymetry_text11a,'visible','on', 'Position', [.3 .63 .1 .1]);
                set(h.bathymetry_edit_mindepthshore,'visible','on', 'Position', [.4 .685 .1  .05]);
                set(h.bathymetry_text11b,'visible','on', 'Position', [.52 .63 .03 .1]);
                set(h.bathymetry_text12,'visible','on');
                set(h.bathymetry_text13,'visible','on');
                set(h.bathymetry_edit_slopeshore,'visible','on');
                set(h.bathymetry_text14,'visible','on');
                set(h.bathymetry_text15,'visible','off');
                set(h.bathymetry_edit_shoreposition,'visible','on');
                set(h.bathymetry_text16,'visible','on');
                set(h.bathymetry_text17,'visible','off');
                set(h.bathymetry_edit_filename,'visible','off');
                set(h.bathymetry_button_load,'visible','off');
                set(h.bathymetry_text18,'visible','on');
                set(h.bathymetry_popup_interpdepth,'visible','on');
                set(h.bathymetry_text19,'visible','on'); % nunu
                set(h.bathymetry_edit_middepth,'visible','on'); % nunu
                set(h.bathymetry_text20,'visible','on'); % nunu

                
            case 'Shore (in y-axis)'
                 set(h.bathymetry_edit_slopeshore,'tooltipString','Specify the slope gradient ')
		        set(h.bathymetry_edit_shoreposition,'tooltipString','Specify a x shore position');
                set(h.bathymetry_edit_middepth,'tooltipString','Specify depth-reference for 3-depth interpolation'); % nunu
                set(h.bathymetry_edit_maxdepthshore,'tooltipString','Specify a maximum depth');
                set(h.bathymetry_edit_mindepthshore,'tooltipString','Specify a minimum depth (i.e 2% Hs)');
                
                set(h.bathymetry_text1, 'Visible', 'off');
                set(h.bathymetry_edit_depth, 'Visible', 'off');
                set(h.bathymetry_text2, 'Visible', 'off');
                set(h.bathymetry_text3,'visible','off');
                set(h.bathymetry_edit_mindepth,'visible','on');
                set(h.bathymetry_text4,'visible','off');
                set(h.bathymetry_text5,'visible','off');
                set(h.bathymetry_edit_maxdepth,'visible','off');
                set(h.bathymetry_text6,'visible','off');
                set(h.bathymetry_text7,'visible','off');
                set(h.bathymetry_edit_slope,'visible','off');
                set(h.bathymetry_text8,'visible','off');
                set(h.bathymetry_text9,'visible','off');
                set(h.bathymetry_edit_startslope,'visible','off');
                set(h.bathymetry_text10,'visible','off');
                set(h.bathymetry_text11a,'visible','on', 'Position', [.3 .63 .1 .1]);
                set(h.bathymetry_edit_mindepthshore,'visible','on', 'Position', [.4 .685 .1  .05]);
                set(h.bathymetry_text11b,'visible','on', 'Position', [.52 .63 .03 .1]);
                set(h.bathymetry_text11,'visible','on');
                set(h.bathymetry_edit_maxdepthshore,'visible','on');
                set(h.bathymetry_text12,'visible','on');
                set(h.bathymetry_text13,'visible','on');
                set(h.bathymetry_edit_slopeshore,'visible','on');
                set(h.bathymetry_text14,'visible','off');
                set(h.bathymetry_text15,'visible','on');
                set(h.bathymetry_edit_shoreposition,'visible','on');
                set(h.bathymetry_text16,'visible','on');
                set(h.bathymetry_text17,'visible','off');
                set(h.bathymetry_edit_filename,'visible','off');
                set(h.bathymetry_button_load,'visible','off');
                set(h.bathymetry_text18,'visible','on');
                set(h.bathymetry_popup_interpdepth,'visible','on');
                set(h.bathymetry_text19,'visible','on'); % nunu
                set(h.bathymetry_edit_middepth,'visible','on'); % nunu
                set(h.bathymetry_text20,'visible','on'); % nunu
               
            case 'User-defined'
                 set(h.bathymetry_edit_mindepthshore,'tooltipString','Specify a minimum depth (i.e 2% Hs)');
                 set(h.bathymetry_edit_middepth,'tooltipString','Specify depth-reference for 3-depth interpolation'); % nunu

                set(h.bathymetry_text1, 'Visible', 'off');
                set(h.bathymetry_edit_depth, 'Visible', 'off');
                set(h.bathymetry_text2, 'Visible', 'off');
                set(h.bathymetry_text3,'visible','off');
                set(h.bathymetry_edit_mindepth,'visible','off');
                set(h.bathymetry_text4,'visible','off');
                set(h.bathymetry_text5,'visible','off');
                set(h.bathymetry_edit_maxdepth,'visible','off');
                set(h.bathymetry_text6,'visible','off');
                set(h.bathymetry_text7,'visible','off');
                set(h.bathymetry_edit_slope,'visible','off');
                set(h.bathymetry_text8,'visible','off');
                set(h.bathymetry_text9,'visible','off');
                set(h.bathymetry_edit_startslope,'visible','off');
                set(h.bathymetry_text10,'visible','off');
                
                set(h.bathymetry_text11,'visible','off');
                set(h.bathymetry_edit_maxdepthshore,'visible','off');
                set(h.bathymetry_text12,'visible','off');
                set(h.bathymetry_text13,'visible','off');
                set(h.bathymetry_edit_slopeshore,'visible','off');
                set(h.bathymetry_text14,'visible','off');
                set(h.bathymetry_text15,'visible','off');
                 set(h.bathymetry_edit_shoreposition,'visible','off');
                 set(h.bathymetry_text16,'visible','off');
                 set(h.bathymetry_text17,'visible','on');
                 set(h.bathymetry_edit_filename,'visible','on');
                 set(h.bathymetry_button_load,'visible','on');
                 set(h.bathymetry_text18,'visible','on');
                 set(h.bathymetry_popup_interpdepth,'visible','on');
                 set(h.bathymetry_text19,'visible','on'); % nunu
                 set(h.bathymetry_edit_middepth,'visible','on'); % nunu
                 set(h.bathymetry_text20,'visible','on'); % nunu
                
%                 set(h.bathymetry_text11a,'visible','off');
%                 set(h.bathymetry_edit_mindepthshore,'visible','off');
%                 set(h.bathymetry_text11b,'visible','off');

                
                dataUser=get(h.bathymetry_button_load,'userdata');
                if ~isempty(dataUser)
                    if any(dataUser.data(:,3)>0)
                        set(h.bathymetry_text11a,'visible','on','enable','on', 'Position', [.04 .55 .1 .1]);
                        set(h.bathymetry_edit_mindepthshore,'visible','on','enable','on', 'Position', [.15 .605 .1  .05]);
                        set(h.bathymetry_text11b,'visible','on', 'enable','on','Position', [.26 .55 .03 .1]);
                    else
                         set(h.bathymetry_text11a,'visible','on','enable','off', 'Position',  [.04 .55 .1 .1]);
                         set(h.bathymetry_edit_mindepthshore,'visible','on','enable','off', 'Position',[.15 .605 .1  .05]);
                         set(h.bathymetry_text11b,'visible','on','enable','off', 'Position', [.26 .55 .03 .1]);

                    end
                else
                   set(h.bathymetry_text11a,'visible','on','enable','off', 'Position', [.04 .55 .1 .1]);
                  set(h.bathymetry_edit_mindepthshore,'visible','on','enable','off', 'Position', [.15 .605 .1  .05]);
                   set(h.bathymetry_text11b,'visible','on','enable','off', 'Position', [.26 .55 .03 .1]);
               end

                              
        end
    end

    function callback_bathy_depth(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_depth,'Userdata',param);
    end

    function callback_bathy_mindepth(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_mindepth,'Userdata',param);
    end

   % nunu
    function callback_bathy_middepth(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_middepth,'Userdata',param);
    end

	% nunu
    function callback_bathy_popup_interpdepth(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=get(hObj,'Value');
		if param > 1
			set(h.bathymetry_edit_middepth,'enable','on')
		else
			set(h.bathymetry_edit_middepth,'enable','off')
		end
        set(h.bathymetry_popup_interpdepth,'Userdata',param);
    end

    function callback_bathy_maxdepth(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_maxdepth,'Userdata',param);
    end

    function callback_bathy_slope(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_slope,'Userdata',param);
    end

    function callback_bathy_startslope(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_startslope,'Userdata',param);
    end

    function callback_bathy_maxdepthshore(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_maxdepthshore,'Userdata',param);
    end

    function callback_bathy_mindepthshore(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_mindepthshore,'Userdata',param);
    end

    function callback_bathy_slopeshore(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_slopeshore,'Userdata',param);
    end

    function callback_bathy_shoreposition(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'String'));
        set(h.bathymetry_edit_shoreposition,'Userdata',param);
    end

    function callback_bathy_userdefined(hObj,eventdata,h)
        [file_name,directory]=uigetfile([h.projectdirectory,'\','*.txt; *.dat; *.mat; *.asc'],...
            'Load bathymetry');
        
        if directory~=0
            set(h.monitorbox,'String','>>loading data','foregroundcolor','k');
            temp=load([directory,file_name]);
            
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            
            if isstruct(my_data)==0
                bath.filename=file_name;
                bath.data=my_data;
                
                if length(my_data(1,:))~=3
                    set(h.monitorbox,'String','>>Wrong format of the bathymetry data',...
                        'foregroundcolor','r');
                else
                     if any(my_data(:,3)>0)
                        set(h.monitorbox,'String',['>> A land (positive value of bathymetry) is detected, This is runup case.'],...
                            'foregroundcolor','r');
                        set(h.bathymetry_edit_mindepthshore,'enable','on')
                     else
                         set(h.bathymetry_edit_mindepthshore,'enable','off')
                     end   
                     
                         set(h.bathymetry_button_load,'Userdata',bath); %record function name in handle
                         set(h.bathymetry_edit_filename,'string',file_name)
                         set(h.monitorbox,'String',['>> ', file_name,' has been loaded'],...
                             'foregroundcolor',[0.8 0.8 0.8]);
               
                end
                
            else
                set(h.monitorbox,'String',['>> Wrong data file'],...
                    'foregroundcolor','r');
            end
            
        else
            set(h.monitorbox,'String',['>> data is not loaded'],...
                'foregroundcolor','r');
        end
        
    end

    function callback_friction_check(hObj,eventdata,h)
        CheckBotFric=get(hObj,'Value');
        
        if CheckBotFric == 1
            set(h.bathymetry_table_friction, 'Visible', 'on');
            set(h.bathymetry_button_addrow, 'Visible', 'on');
            set(h.bathymetry_button_deleterow,'Visible','on');
        else
            set(h.bathymetry_table_friction, 'Visible', 'off');
            set(h.bathymetry_button_addrow, 'Visible', 'off');
            set(h.bathymetry_button_deleterow,'Visible','off');
        end
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_friction_addrow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        oldData = get(h.bathymetry_table_friction,'Data');
        newData = [oldData; {'Rectangle' 0 0 0 0 0}];
        set(h.bathymetry_table_friction,'Data',newData);
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_friction_deleterow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        oldData = get(h.bathymetry_table_friction,'Data');
        Ndata=length(oldData(:,1));
        if Ndata>1
            newData = oldData(1:Ndata-1,:);
            set(h.bathymetry_table_friction,'Data',newData);
        end
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_ivp_popup(hObj,eventdata,h)
        str = get(hObj, 'String');
        val = get(hObj,'Value');
        
        switch str{val}
            case 'Zero'
                set(h.waveinput_ivp_text1, 'Visible', 'off');
                set(h.waveinput_ivp_edit_A, 'Visible', 'off');
                set(h.waveinput_ivp_text2, 'Visible', 'off');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'off');
                set(h.waveinput_ivp_text3, 'Visible', 'off');
                set(h.waveinput_ivp_text3a,'visible','off');
                set(h.waveinput_ivp_text4, 'Visible', 'off');
                set(h.waveinput_ivp_text5, 'Visible', 'off');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'off');
                set(h.waveinput_ivp_text6, 'Visible', 'off');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'off');
                set(h.waveinput_ivp_button_load, 'Visible', 'off');
            case 'Gaussian (along x-axis)'
                set(h.waveinput_ivp_text1, 'Visible', 'on');
                set(h.waveinput_ivp_edit_A, 'Visible', 'on');
                set(h.waveinput_ivp_text2, 'Visible', 'on');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'on');
                set(h.waveinput_ivp_text3, 'Visible', 'on');
                set(h.waveinput_ivp_text3a,'visible','off');
                set(h.waveinput_ivp_text4, 'Visible', 'on');
                set(h.waveinput_ivp_text5, 'Visible', 'on');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'on');
                set(h.waveinput_ivp_text6, 'Visible', 'off');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'off');
                set(h.waveinput_ivp_button_load, 'Visible', 'off');
            case 'Gaussian (along y-axis)'
                set(h.waveinput_ivp_text1, 'Visible', 'on');
                set(h.waveinput_ivp_edit_A, 'Visible', 'on');
                set(h.waveinput_ivp_text2, 'Visible', 'on');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'on');
                set(h.waveinput_ivp_text3, 'Visible', 'off');
                set(h.waveinput_ivp_text3a,'visible','on');
                set(h.waveinput_ivp_text4, 'Visible', 'on');
                set(h.waveinput_ivp_text5, 'Visible', 'on');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'on');
                set(h.waveinput_ivp_text6, 'Visible', 'off');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'off');
                set(h.waveinput_ivp_button_load, 'Visible', 'off');
            case 'N-Wave (along x-axis)'
                set(h.waveinput_ivp_text1, 'Visible', 'on');
                set(h.waveinput_ivp_edit_A, 'Visible', 'on');
                set(h.waveinput_ivp_text2, 'Visible', 'on');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'on');
                set(h.waveinput_ivp_text3, 'Visible', 'on');
                set(h.waveinput_ivp_text3a,'visible','off');
                set(h.waveinput_ivp_text4, 'Visible', 'on');
                set(h.waveinput_ivp_text5, 'Visible', 'on');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'on');
                set(h.waveinput_ivp_text6, 'Visible', 'off');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'off');
                set(h.waveinput_ivp_button_load, 'Visible', 'off');
            case 'N-Wave (along y-axis)'
                set(h.waveinput_ivp_text1, 'Visible', 'on');
                set(h.waveinput_ivp_edit_A, 'Visible', 'on');
                set(h.waveinput_ivp_text2, 'Visible', 'on');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'on');
                set(h.waveinput_ivp_text3, 'Visible', 'off');
                set(h.waveinput_ivp_text3a,'visible','on');
                set(h.waveinput_ivp_text4, 'Visible', 'on');
                set(h.waveinput_ivp_text5, 'Visible', 'on');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'on');
                set(h.waveinput_ivp_text6, 'Visible', 'off');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'off');
                set(h.waveinput_ivp_button_load, 'Visible', 'off');
            case 'User-defined'
                set(h.waveinput_ivp_text1, 'Visible', 'off');
                set(h.waveinput_ivp_edit_A, 'Visible', 'off');
                set(h.waveinput_ivp_text2, 'Visible', 'off');
                set(h.waveinput_ivp_edit_centre_position, 'Visible', 'off');
                set(h.waveinput_ivp_text3, 'Visible', 'off');
                set(h.waveinput_ivp_text3a,'visible','off');
                set(h.waveinput_ivp_text4, 'Visible', 'off');
                set(h.waveinput_ivp_text5, 'Visible', 'off');
                set(h.waveinput_ivp_edit_sd, 'Visible', 'off');
                set(h.waveinput_ivp_text6, 'Visible', 'on');
                set(h.waveinput_ivp_edit_filename, 'Visible', 'on');
                set(h.waveinput_ivp_button_load, 'Visible', 'on');
        end
    end

    function callback_ivp_amplitude(hObj,eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.waveinput_ivp_edit_A,'Userdata',param);
    end

    function callback_ivp_centreposition(hObj,eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.waveinput_ivp_edit_centre_position,'Userdata',param);
    end


    function callback_ivp_stdev(hObj,eventdata,h)
        param=str2num(get(hObj,'String'));
        set(h.waveinput_ivp_edit_sd,'Userdata',param);
    end


    function callback_ivp_loaddata(hObj,eventdata,h)
        
        [file_name,directory]=uigetfile([h.projectdirectory,'\','*.txt; *.dat; *.mat; *.asc'],...
            'Load initial condition');
        
        if directory~=0
            set(h.monitorbox,'String','>>loading data','foregroundcolor','k');
            temp=load([directory,file_name]);
            
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
            
            if isstruct(my_data)==0
                if length(my_data(1,:))<3
                    set(h.monitorbox,'String','>>Wrong format of the initial condition data',...
                        'foregroundcolor','r');
                else
                    IVPdat.data=my_data;
                    IVPdat.name=file_name;
                    set(h.waveinput_ivp_button_load,'Userdata',IVPdat); %record function name in handle
                    set(h.waveinput_ivp_edit_filename,'string',file_name);
                    set(h.monitorbox,'String',[file_name,' has been loaded'],...
                        'foregroundcolor','k');
                end
            else
                set(h.monitorbox,'String','>>Wrong format of the initial condition data',...
                    'foregroundcolor','r');
            end
        else
            set(h.monitorbox,'String','>>data is not loaded','foregroundcolor','r');
        end
    end

    function callback_influx_popup(hObj,eventdata,h)
        str = get(hObj, 'String');
        val = get(hObj,'Value');
        
        
        switch str{val}
            case 'No'
                set(h.waveinput_influx_checkbox_ramp,'visible','off',...
                    'enable','off')
                set(h.waveinput_influx_checkbox_nonlinadj,'visible','off',...
                    'enable','off')
                set(h.waveinput_influx_edit_ramp,'visible','off')
                set(h.waveinput_influx_edit_nonlinadj,'visible','off')
                set(h.waveinput_influx_text_ramp,'visible','off');
                set(h.waveinput_influx_text_nonlinadj,'visible','off');
                set(h.waveinput_influxline_checkbox_ramp,'visible','off',...
                    'enable','off')
                set(h.waveinput_influxline_edit_ramp,'visible','off')
                set(h.waveinput_influxline_text_ramp,'visible','off');
                
                set(h.waveinput_influx_button_addrow,'visible','off')
                set(h.waveinput_influx_button_deleterow,'visible','off')
                
                set(h.waveinput_influx_proptable,'visible','off')
                set(h.waveinput_influx_methodtable,'visible','off');
                
            case 'Yes'
                set(h.waveinput_influx_checkbox_ramp,'visible','on',...
                    'enable','on')
                if get(h.waveinput_influx_checkbox_ramp,'value')==1
                    set(h.waveinput_influx_edit_ramp,'visible','on','enable','on')
                    set(h.waveinput_influx_text_ramp,'visible','on','enable','on');
                else
                    set(h.waveinput_influx_edit_ramp,'visible','on','enable','off')
                    set(h.waveinput_influx_text_ramp,'visible','on','enable','off');
                end
                
                IdHS=get(h.model_popup_dynamic,'value');
                if IdHS>1
                    set(h.waveinput_influx_checkbox_nonlinadj,'visible','on',...
                        'enable','on')
                    if get(h.waveinput_influx_checkbox_nonlinadj,'value')==1
                        set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','on')
                        set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','on');
                    else
                        set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','off')
                        set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','off');
                    end
                else
                    set(h.waveinput_influx_checkbox_nonlinadj,'visible','on',...
                        'enable','off')
                    set(h.waveinput_influx_edit_nonlinadj,'visible','on','enable','off')
                    set(h.waveinput_influx_text_nonlinadj,'visible','on','enable','off');
                end
                
                
                set(h.waveinput_influxline_checkbox_ramp,'visible','on',...
                    'enable','on')
                if get(h.waveinput_influxline_checkbox_ramp,'value')==1
                    set(h.waveinput_influxline_edit_ramp,'visible','on','enable','on')
                    set(h.waveinput_influxline_text_ramp,'visible','on','enable','on');
                else
                    set(h.waveinput_influxline_edit_ramp,'visible','on','enable','off')
                    set(h.waveinput_influxline_text_ramp,'visible','on','enable','off');
                end
                
                set(h.waveinput_influx_button_addrow,'visible','on')
                set(h.waveinput_influx_button_deleterow,'visible','on')
                
                set(h.waveinput_influx_proptable,'visible','on')
                set(h.waveinput_influx_methodtable,'visible','on');
                
                
        end
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a

    end



    function callback_influx_addrow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        oldData1 = get(h.waveinput_influx_proptable,'Data');
        newData1 = [oldData1; {'Harmonic', 0, '-', 0, '-', '-', 0}];
        set(h.waveinput_influx_proptable,'Data',newData1);
        
        oldData2 = get(h.waveinput_influx_methodtable,'Data');
        newData2 = [oldData2;{'Area', 'Straight', 0, 0, 0, 0, '-','-','-','-','-', 0, 0, 0}];
        set(h.waveinput_influx_methodtable,'Data',newData2);
        guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
    end

    function callback_influx_deleterow(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Ind_del=get(hObj,'userdata');
        if ~isempty(Ind_del)
            if Ind_del.Id==1
                oldData1 = get(h.waveinput_influx_proptable,'Data');
                oldData2 = get(h.waveinput_influx_methodtable,'Data');
                
                mask = (1:size(oldData1,1))';
                mask(Ind_del.row) = [];
                
                newData1 = oldData1(mask,:);
                set(h.waveinput_influx_proptable,'Data',newData1);
                newData2 = oldData2(mask,:);
                set(h.waveinput_influx_methodtable,'Data',newData2);
                
                guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            else
                set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
            end
        else
            set(h.monitorbox,'String','>>Please select a row in the table to be deleted','foregroundcolor','r');
        end
     end



    function callback_influx_check_ramp(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Id=get(hObj,'value');
        if Id==1
            set(h.waveinput_influx_edit_ramp,'enable','on');
            set(h.waveinput_influx_text_ramp,'enable','on');
            
        else
            set(h.waveinput_influx_edit_ramp,'enable','off');
            set(h.waveinput_influx_text_ramp,'enable','off');
        end
    end
    function callback_influx_ramp(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        if length(param)==1
            if param <=0
                set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
                uicontrol(h.waveinput_influx_edit_ramp);
            else
                set(h.waveinput_influx_edit_ramp,'userdata',param);
            end
        else
            set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
            uicontrol(h.waveinput_influx_edit_ramp);
        end
    end

    function callback_influxline_check_ramp(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Id=get(hObj,'value');
        if Id==1
            set(h.waveinput_influxline_edit_ramp,'enable','on');
            set(h.waveinput_influxline_text_ramp,'enable','on');
            
        else
            set(h.waveinput_influxline_edit_ramp,'enable','off');
            set(h.waveinput_influxline_text_ramp,'enable','off');
        end
    end
    function callback_influxline_ramp(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        if length(param)==1
            if param <=0
                set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
                uicontrol(h.waveinput_influxline_edit_ramp);
            else
                set(h.waveinput_influxline_edit_ramp,'userdata',param);
            end
        else
            set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
            uicontrol(h.waveinput_influxline_edit_ramp);
        end
    end

    function callback_influx_check_nonlinadj(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        Id=get(hObj,'value');
        if Id==1
            set(h.waveinput_influx_edit_nonlinadj,'enable','on');
            set(h.waveinput_influx_text_nonlinadj,'enable','on');
        else
            set(h.waveinput_influx_edit_nonlinadj,'enable','off');
            set(h.waveinput_influx_text_nonlinadj,'enable','off');
        end
    end

    function callback_influx_nonlinadj(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        if length(param)==1
            if param <=0
                set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
                uicontrol(h.waveinput_influx_edit_nonlinadj);
            else
                set(h.waveinput_influx_edit_nonlinadj,'userdata',param);
            end
        else
            set(h.monitorbox,'String','>> wrong input','foregroundcolor','r');
            uicontrol(h.waveinput_influx_edit_nonlinadj);
        end
    end

    function callback_edit_param(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param);
    end

    function callback_time_tinit(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.time_edit_interval_init,'userdata',param);
    end

    function callback_time_tend(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.time_edit_interval_end,'userdata',param);
    end

    function callback_time_dt(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.time_edit_step,'userdata',param);
    end

    function callback_option_inflow_check(hObj,eventdata,h)
    Id=get(hObj,'value');
    if Id==1
     set(h.options_edit_interval_init,'enable','on');
     set(h.options_edit_interval_end,'enable','on');
     set(h.options_edit_step,'enable','on');
    else
      set(h.options_edit_interval_init,'enable','off');
     set(h.options_edit_interval_end,'enable','off');
     set(h.options_edit_step,'enable','off');  
    end
    end

    function callback_option_inflow_tinit(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.options_edit_interval_init,'userdata',param);
    end

    function callback_option_inflow_tend(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.options_edit_interval_end,'userdata',param);
    end

    function callback_option_inflow_dtstep(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.options_edit_step,'userdata',param);
    end

    function callback_ode_partition(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.options_edit_partition,'enable','off')
            set(h.options_text6,'enable','off');
        else
            set(h.options_edit_partition,'enable','on')
            set(h.options_text6,'enable','on');
        end
    end
    function callback_option_partition(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.options_edit_partition,'userdata',param);
    end

    function callback_ode_tol(hObj,eventdata,h)
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        param=str2num(get(hObj,'string'));
        set(h.options_edit_ode_tol,'userdata',param);
    end

    function callback_mc_cb_check(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.options_montecarlo_cb_waveinput,'enable','on','value',0);
            set(h.options_montecarlo_cb_numset,'enable','on','value',0);
            
            inflData=get(h.waveinput_influx_proptable,'data');
            Ninflux=length(inflData(:,1));
            IDflag=0;
            for ii=1:Ninflux
                if strcmpi(inflData(ii,1),'JONSWAP')||...
                        strcmpi(inflData(ii,1),'User-defined (variance density spectrum)')
                    IDflag=1;
                end
            end
            if IDflag==1
                set(h.options_montecarlo_cb_Nsimul,'enable','on');
            else
                set(h.options_montecarlo_cb_Nsimul,'enable','off','value',0);
            end
        else
            set(h.options_montecarlo_cb_waveinput,'enable','off','value',0);
            set(h.options_montecarlo_cb_numset,'enable','off','value',0);
            set(h.options_montecarlo_cb_Nsimul,'enable','off','value',0);
            set(h.options_montecarlo_edit_Nsimul,'enable','off','string','');
        end
        set(h.options_montecarlo_pb_waveinput,'enable','off');
        set(h.options_montecarlo_popup_waveinput,'enable','off','string','--')
        set(h.options_montecarlo_cb_A ,'enable','off')
        set(h.options_montecarlo_edit_A,'enable','off')
        set(h.options_montecarlo_cb_Tp,'enable','off')
        set(h.options_montecarlo_edit_Tp,'enable','off')
        set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
        set(h.options_montecarlo_edit_gamma,'enable','off')
        set(h.options_montecarlo_cb_s,'enable','off','visible','off')
        set(h.options_montecarlo_edit_s,'enable','off','visible','off')
        set(h.options_montecarlo_text_px,'enable','off')
        set(h.options_montecarlo_edit_px,'enable','off')
        set(h.options_montecarlo_text_py,'enable','off')
        set(h.options_montecarlo_edit_py,'enable','off')
    end

    function callback_mc_N_run(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.options_montecarlo_edit_Nsimul,'enable','on')
            set(h.options_montecarlo_cb_numset,'value',0,'enable','off');
            set(h.options_montecarlo_cb_waveinput,'value',0,'enable','off')
        else
            set(h.options_montecarlo_edit_Nsimul,'enable','off')
            set(h.options_montecarlo_cb_numset,'value',0,'enable','on');
            set(h.options_montecarlo_cb_waveinput,'value',0,'enable','on')
        end
        callback_mc_cb_waveinput(h.options_montecarlo_cb_waveinput,[],h);
        callback_mc_cb_numset(h.options_montecarlo_cb_numset,[],h);
    end

    function callback_mc_cb_numset(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.options_montecarlo_text_px,'enable','on')
            set(h.options_montecarlo_edit_px,'enable','on')
            set(h.options_montecarlo_text_py,'enable','on')
            set(h.options_montecarlo_edit_py,'enable','on')
        else
            set(h.options_montecarlo_text_px,'enable','off')
            set(h.options_montecarlo_edit_px,'enable','off')
            set(h.options_montecarlo_text_py,'enable','off')
            set(h.options_montecarlo_edit_py,'enable','off')
        end
    end

    function callback_mc_cb_waveinput(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        Id=get(hObj,'value');
        inflData=get(h.waveinput_influx_proptable,'data');
        Ninflux=length(inflData(:,1));
        if Id==1
            if Ninflux>0
                for I=1:Ninflux
                    Str{I}=['#',num2str(I)];
                end
                
                set(h.options_montecarlo_popup_waveinput,'enable','on','string',Str);
                set(h.options_montecarlo_pb_waveinput,'enable','on');
                if strcmpi(inflData(1,1),'Harmonic')
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Amplitude')
                    set(h.options_montecarlo_edit_A,'enable','off')
                    set(h.options_montecarlo_cb_Tp,'enable','on')
                    set(h.options_montecarlo_edit_Tp,'enable','off')
                    set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
                    set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
                    set(h.options_montecarlo_cb_s,'enable','off','visible','off')
                    set(h.options_montecarlo_edit_s,'enable','off','visible','off')
                elseif strcmpi(inflData(1,1),'JONSWAP')
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Hs')
                    set(h.options_montecarlo_edit_A,'enable','off')
                    set(h.options_montecarlo_cb_Tp,'enable','on')
                    set(h.options_montecarlo_edit_Tp,'enable','off')
                    set(h.options_montecarlo_cb_gamma,'enable','on','visible','on')
                    set(h.options_montecarlo_edit_gamma,'enable','off','visible','on')
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on')
                    set(h.options_montecarlo_edit_s,'enable','off','visible','on')
                elseif strcmpi(inflData(1,1),'User-defined (variance density spectrum)')
                    set(h.options_montecarlo_cb_A ,'enable','off','string','Hs')
                    set(h.options_montecarlo_edit_A,'enable','off')
                    set(h.options_montecarlo_cb_Tp,'enable','off')
                    set(h.options_montecarlo_edit_Tp,'enable','off')
                    set(h.options_montecarlo_cb_gamma,'enable','off','visible','on')
                    set(h.options_montecarlo_edit_gamma,'enable','off','visible','on')
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on')
                    set(h.options_montecarlo_edit_s,'enable','off','visible','on')
                else
                    set(h.options_montecarlo_cb_A ,'enable','off','string','Amplitude')
                    set(h.options_montecarlo_edit_A,'enable','off')
                    set(h.options_montecarlo_cb_Tp,'enable','off')
                    set(h.options_montecarlo_edit_Tp,'enable','off')
                    set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
                    set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
                    set(h.options_montecarlo_cb_s,'enable','off','visible','off')
                    set(h.options_montecarlo_edit_s,'enable','off','visible','off')
                end
            else
                set(h.options_montecarlo_pb_waveinput,'enable','off');
                set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave influx data!')
                set(callback_mc_cb_waveinput,'value',0)
                set(h.options_montecarlo_popup_waaveinput,'enable','off','string','--')
                set(h.options_montecarlo_cb_A ,'enable','off','string','Amplitude')
                set(h.options_montecarlo_edit_A,'enable','off')
                set(h.options_montecarlo_cb_Tp,'enable','off')
                set(h.options_montecarlo_edit_Tp,'enable','off')
                set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
                set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
                set(h.options_montecarlo_cb_s,'enable','off','visible','off')
                set(h.options_montecarlo_edit_s,'enable','off','visible','off')
            end
        else
            set(h.options_montecarlo_pb_waveinput,'enable','off');
            set(h.options_montecarlo_cb_A ,'enable','off','string','Amplitude')
            set(h.options_montecarlo_cb_A ,'enable','off','value',0)
            set(h.options_montecarlo_edit_A,'enable','off')
            set(h.options_montecarlo_cb_Tp,'enable','off','value',0)
            set(h.options_montecarlo_edit_Tp,'enable','off')
            set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_cb_s,'enable','off','visible','off')
            set(h.options_montecarlo_edit_s,'enable','off','visible','off')
            set(h.options_montecarlo_popup_waveinput,'enable','off','value',1,'string','--','visible','on')
            
        end
        
    end

    function callback_mc_pb_storedata(hObj,eventdata,h)
        Id=get(h.options_montecarlo_popup_waveinput,'value');
        inflData=get(h.waveinput_influx_proptable,'data');
        data=get(h.options_montecarlo_pb_waveinput,'userdata');
        
        set(h.monitorbox,'foregroundcolor','k','string','>>');
        
        if strcmpi(inflData(Id,1),'Harmonic')
            data(Id).type='Harmonic';
            data(Id).A_check=get(h.options_montecarlo_cb_A,'value');
            if data(Id).A_check==1
                data(Id).A_param=get(h.options_montecarlo_edit_A,'userdata');
                if isempty(data(Id).A_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify amplitude(s)']);
                    return;
                end
            end
            data(Id).Tp_check=get(h.options_montecarlo_cb_Tp,'value');
            if data(Id).Tp_check==1
                data(Id).Tp_param=get(h.options_montecarlo_edit_Tp,'userdata');
                if isempty(data(Id).Tp_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify period(s)']);
                    return;
                end
            end
        elseif strcmpi(inflData(Id,1),'JONSWAP')
            data(Id).Hs_check=get(h.options_montecarlo_cb_A,'value');
            if data(Id).Hs_check==1
                data(Id).Hs_param=get(h.options_montecarlo_edit_A,'userdata');
                if isempty(data(Id).Hs_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify Hs']);
                    return;
                end
            end
            data(Id).Tp_check=get(h.options_montecarlo_cb_Tp,'value');
            if data(Id).Tp_check==1
                data(Id).Tp_param=get(h.options_montecarlo_edit_Tp,'userdata');
                if isempty(data(Id).Tp_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify period(s)']);
                    return;
                end
            end
            data(Id).gamma_check=get(h.options_montecarlo_cb_gamma,'value');
            if data(Id).gamma_check==1
                data(Id).gamma_param=get(h.options_montecarlo_edit_gamma,'userdata');
                if isempty(data(Id).gamma_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify gamma']);
                    return;
                end
            end
            
            data(Id).s_check=get(h.options_montecarlo_cb_s,'value');
            if data(Id).s_check==1
                data(Id).s_param=get(h.options_montecarlo_edit_s,'userdata');
                if isempty(data(Id).s_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify s']);
                    return;
                end
            end
        elseif strcmpi(inflData(Id,1),'User-defined (variance density spectrum)')
            data(Id).s_check=get(h.options_montecarlo_cb_s,'value');
            if data(Id).s_check==1
                data(Id).s_param=get(h.options_montecarlo_edit_s,'userdata');
                if isempty(data(Id).s_param)
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Specify s']);
                    return;
                end
            end
        end
        set(hObj,'userdata',data);
    end

    function callback_mc_popup_waveinput(hObj,eventdata,h)
        Id=get(hObj,'value');
        inflData=get(h.waveinput_influx_proptable,'data');
        data=get(h.options_montecarlo_pb_waveinput,'userdata');
        
        set(h.monitorbox,'foregroundcolor','k','string','>>');
        if ~isempty(data)
            if Id>1
                if isempty(data(Id-1))
                    set(hObj,'value',Id-1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>> Please store data for waveinput #',num2str(Id-1)]);
                    return;
                end
            end
        else
            set(hObj,'value',1);
            set(h.monitorbox,'foregroundcolor','r','string',['>> Please store data for waveinput #1']);
            return;
        end
        
        set(h.monitorbox,'foregroundcolor','k','string','>>');
        
        if strcmpi(inflData(Id,1),'Harmonic')
            if length(data)>=Id
                if isfield(data(Id),'A_param')
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Amplitude','value',data(Id).A_check)
                    if data(Id).A_check==1
                        Ndat=length(data(Id).A_param(:));
                        if Ndat>1
                            for i=1:Ndat-1
                                strA{i}=[num2str(data(Id).A_param(i)),';',];
                            end
                            strA{Ndat}=num2str(data(Id).A_param(Ndat));
                        else
                            strA{Ndat}=num2str(data(Id).A_param(1));
                        end
                        set(h.options_montecarlo_edit_A,'enable','on','userdata',data(Id).A_param,...
                            'string',[strA{1:Ndat}])
                    else
                        set(h.options_montecarlo_edit_A,'enable','off','string',...
                            '','userdata',[],'visible','on')
                    end
                else
                    set(h.options_montecarlo_edit_A,'enable','off','string','','userdata',[])
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Amplitude','value',0)
                end
            else
                set(h.options_montecarlo_edit_A,'enable','off','string','','userdata',[])
                set(h.options_montecarlo_cb_A ,'enable','on','string','Amplitude','value',0)
            end
            
            if length(data)>=Id
                if isfield(data(Id),'Tp_param')
                    set(h.options_montecarlo_cb_Tp ,'enable','on','value',data(Id).Tp_check)
                    if data(Id).Tp_check==1
                        Ndat=length(data(Id).Tp_param(:));
                        if Ndat>1
                            for i=1:Ndat-1
                                strTp{i}=[num2str(data(Id).Tp_param(i)),';',];
                            end
                            strTp{Ndat}=num2str(data(Id).Tp_param(Ndat));
                        else
                            strTp{Ndat}=num2str(data(Id).Tp_param(1));
                        end
                        set(h.options_montecarlo_edit_Tp,'enable','on','string',...
                            [strTp{1:Ndat}],'userdata',data(Id).Tp_param)
                    else
                        set(h.options_montecarlo_edit_Tp,'enable','off','string',...
                            '','userdata',[],'visible','on')
                    end
                else
                    set(h.options_montecarlo_edit_Tp,'enable','off','string','','userdata',[])
                    set(h.options_montecarlo_cb_Tp,'enable','on','value',0);
                end
            else
                set(h.options_montecarlo_edit_Tp,'enable','off','string','','userdata',[])
                set(h.options_montecarlo_cb_Tp,'enable','on','value',0);
            end
            set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_cb_s,'enable','off','visible','off')
            set(h.options_montecarlo_edit_s,'enable','off','visible','off')
        elseif strcmpi(inflData(Id,1),'JONSWAP')
            if length(data)>=Id
                if isfield(data(Id),'Hs_param')
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Hs','value',data(Id).Hs_check)
                    if data(Id).Hs_check==1
                        Ndat=length(data(Id).Hs_param(:));
                        if Ndat>1
                            for i=1:Ndat-1
                                strA{i}=[num2str(data(Id).Hs_param(i)),';',];
                            end
                            strA{Ndat}=num2str(data(Id).Hs_param(Ndat));
                        else
                            strA{Ndat}=num2str(data(Id).Hs_param(1));
                        end
                        
                        set(h.options_montecarlo_edit_A,'enable','on',...
                            'string',[strA{1:Ndat}],'userdata',data(Id).Hs_param)
                    else
                        set(h.options_montecarlo_edit_A,'enable','off',...
                            'string','','userdata',[])
                    end
                else
                    set(h.options_montecarlo_edit_A,'enable','off',...
                        'string','','userdata',[])
                    set(h.options_montecarlo_cb_A ,'enable','on','string','Hs','value',0)
                end
            else
                set(h.options_montecarlo_edit_A,'enable','off',...
                    'string','','userdata',[])
                set(h.options_montecarlo_cb_A ,'enable','on','string','Hs','value',0)
            end
            if length(data)>=Id
                if isfield(data(Id),'Tp_param')
                    set(h.options_montecarlo_cb_Tp,'enable','on','value',data(Id).Tp_check);
                    if data(Id).Tp_check==1
                        Ndat=length(data(Id).Tp_param(:));
                        if Ndat>1
                            for i=1:Ndat-1
                                strTp{i}=[num2str(data(Id).Tp_param(i)),';',];
                            end
                            strTp{Ndat}=num2str(data(Id).Tp_param(Ndat));
                        else
                            strTp{Ndat}=num2str(data(Id).Tp_param(1));
                        end
                        
                        set(h.options_montecarlo_edit_Tp,'enable','on','string',...
                            [strTp{1:Ndat}],'userdata',data(Id).Tp_param)
                    else
                        set(h.options_montecarlo_edit_Tp,'enable','off','string',...
                            '','userdata',[])
                    end
                else
                    set(h.options_montecarlo_edit_Tp,'enable','off','string','','userdata',[])
                    set(h.options_montecarlo_cb_Tp,'enable','on','value',0);
                end
            else
                set(h.options_montecarlo_edit_Tp,'enable','off','string','','userdata',[])
                set(h.options_montecarlo_cb_Tp,'enable','on','value',0);
            end
            
            if length(data)>=Id
                if isfield(data(Id),'gamma_param')
                    set(h.options_montecarlo_cb_gamma,'enable','on','visible','on','value',data(Id).gamma_check)
                    if data(Id).gamma_check==1
                        Ndat=length(data(Id).gamma_param(:));
                        if Ndat>1
                            for i=1:Ndat-1
                                strG{i}=[num2str(data(Id).gamma_param(i)),';',];
                            end
                            strG{Ndat}=num2str(data(Id).gamma_param(Ndat));
                        else
                            strG{Ndat}=num2str(data(Id).gamma_param(1));
                        end
                        set(h.options_montecarlo_edit_gamma,'enable','on','visible','on',...
                            'string',[strG{1:Ndat}],'userdata',data(Id).gamma_param)
                    else
                        set(h.options_montecarlo_edit_gamma,'enable','off','visible','on',...
                            'string','','userdata',[])
                    end
                else
                    set(h.options_montecarlo_edit_gamma,'enable','off','visible','on',...
                        'string','','userdata',[])
                    set(h.options_montecarlo_cb_gamma,'enable','on','visible','on','value',0)
                end
            else
                set(h.options_montecarlo_edit_gamma,'enable','off','visible','on',...
                    'string','','userdata',[])
                set(h.options_montecarlo_cb_gamma,'enable','on','visible','on','value',0)
            end
            
            if length(data)>=Id
                if isfield(data(Id),'s_param')
                    Ndat=length(data(Id).s_param(:));
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',data(Id).s_check);
                    if data(Id).s_check==1
                        if Ndat>1
                            for i=1:Ndat-1
                                strS{i}=[num2str(data(Id).s_param(i)),';',];
                            end
                            strS{Ndat}=num2str(data(Id).s_param(Ndat));
                        else
                            strS{Ndat}=num2str(data(Id).s_param(1));
                        end
                        set(h.options_montecarlo_edit_s,'enable','on','visible','on',...
                            'string',[strS{1:Ndat}],'userdata',data(Id).s_param);
                    else
                        set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                            'string','','userdata',[]);
                    end
                else
                    set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                        'string','','userdata',[])
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',0);
                end
            else
                set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                    'string','','userdata',[])
                set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',0);
            end
        elseif strcmpi(inflData(1,1),'User-defined (variance density spectrum)')
            if length(data)>=Id
                if isfield(data(Id),'s_param')
                    Ndat=length(data(Id).s_param(:));
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',data(Id).s_check);
                    if data(Id).s_check==1
                        if Ndat>1
                            for i=1:Ndat-1
                                strS{i}=[num2str(data(Id).s_param(i)),';',];
                            end
                            strS{Ndat}=num2str(data(Id).s_param(Ndat));
                        else
                            strS{Ndat}=num2str(data(Id).s_param(1));
                        end
                        set(h.options_montecarlo_edit_s,'enable','on','visible','on',...
                            'string',[strS{1:Ndat}],'userdata',data(Id).s_param);
                    else
                        set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                            'string','','userdata',[]);
                    end
                else
                    set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                        'string','','userdata',[])
                    set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',0);
                end
            else
                set(h.options_montecarlo_edit_s,'enable','off','visible','on',...
                    'string','','userdata',[])
                set(h.options_montecarlo_cb_s,'enable','on','visible','on','value',0);
            end
        else
            set(h.options_montecarlo_cb_A ,'enable','off','string','Amplitude')
            set(h.options_montecarlo_edit_A,'enable','off')
            set(h.options_montecarlo_cb_Tp,'enable','off')
            set(h.options_montecarlo_edit_Tp,'enable','off')
            set(h.options_montecarlo_cb_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_edit_gamma,'enable','off','visible','off')
            set(h.options_montecarlo_cb_s,'enable','off','visible','off')
            set(h.options_montecarlo_edit_s,'enable','off','visible','off')
        end
        set(hObj,'visible','on')
    end

    function callback_mc_cb_param(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on','visible','on');
        else
            set(href,'enable','off','visible','on');
        end
    end

    function callback_mc_edit_param(hObj,eventdata)
        param=get(hObj,'string');
        set(hObj,'userdata',str2num(param));
    end

    function callback_preproc(hObj,eventdata,h)
        [statusbarObj]=funGui_JavaFrame_handling();
        statusbarObj.setText('');
        set(h.monitorbox,'String','>>preprocessing...','foregroundcolor','k');
        pause(0.001);
        %%%% passing parameter %%%%%%%%%%%%%
        GUIinput=fun_passing_handles_preproc(h);
        GUIinput.modelphiForm=ID_model_phiform;
        global postproc_flag
        postproc_flag=0;
        %%%%%%%Control handles%%%%%%%%%%%%%%%%%%%%%%%%
        set(h.monitorbox,'String','>>','foregroundcolor','k');
        
        if isempty(GUIinput.proj.name)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.project_edit_name);
            return;
        end
        
        if isempty(GUIinput.proj.projdir) || strcmpi(GUIinput.proj.projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.project_popup_projdir);
            return;
        end
        
            
        if strcmpi(GUIinput.modelbreak.check,'yes')
            if strcmp(GUIinput.modeldyn,'HS1')
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Use only a nonlinear model for breaking wave simulation!')
                uicontrol(h.model_popup_dynamic);
                return;
            end
            
            if isempty(GUIinput.modelbreak.initiation)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify an initiation breaking coef. !')
                uicontrol(h.model_break_edit_initiation);
                return;
            end
            
            if GUIinput.modelbreak.initiation<0.4
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify an initiation breaking coef. in [0.6;1] !')
                uicontrol(h.model_break_edit_initiation);
                return;
            end
            
            
            if isempty(GUIinput.modelbreak.termination)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a termination breaking coef. !')
                uicontrol(h.model_break_edit_termination);
                return;
            end
            
            if GUIinput.modelbreak.termination<0 || GUIinput.modelbreak.termination>1
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> The termination breaking coef. must be in interval [0;1] !')
                uicontrol(h.model_break_edit_termination);
                return;
            end
            
            if isempty(GUIinput.modelbreak.Tstar)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a T_star breaking coef. !')
                uicontrol(h.model_break_edit_Tstar);
                return;
            end
            
            if GUIinput.modelbreak.Tstar<0 || GUIinput.modelbreak.Tstar>1
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> The  T_star breaking coef. must be in interval [0;1] !')
                uicontrol(h.model_break_edit_Tstar);
                return;
            end
            
        end
        
        
        if GUIinput.modelcurrent.check==1
            if isempty(GUIinput.modelcurrent.ux)||length(GUIinput.modelcurrent.ux)~=1
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a velocity in x-direction')
                uicontrol(h.model_current_ux_edit);
                return;
            end
            
            if isempty(GUIinput.modelcurrent.uy)||length(GUIinput.modelcurrent.uy)~=1
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a velocity in y-direction')
                uicontrol(h.model_current_uy_edit);
                return;
            end
        end
        
        if isempty(GUIinput.spatialinterv.xmin)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a x minimum in spatial-interval');
            uicontrol(h.spatial_edit_xmin);
            return;
        end
        
        
        if isempty(GUIinput.spatialinterv.xmax)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a x maximum in spatial-interval');
            uicontrol(h.spatial_edit_xmax);
            return;
        end
        
        
        if GUIinput.spatialinterv.xmin>GUIinput.spatialinterv.xmax
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> x min> x max. Specify a correct input for x min and x max in spatial-interval');
            uicontrol(h.spatial_edit_xmin);
            return;
        end
        
        if isempty(GUIinput.spatialinterv.ymin)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a y minimum in spatial-interval');
            uicontrol(h.spatial_edit_ymin);
            return;
        end
        
        if isempty(GUIinput.spatialinterv.ymax)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a y maximum in spatial-interval');
            uicontrol(h.spatial_edit_ymax);
            return;
        end
        
        if GUIinput.spatialinterv.ymin>GUIinput.spatialinterv.ymax
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>  y min> y max. Specify a correct input for y min and y max in spatial-interval');
            uicontrol(h.spatial_edit_ymin);
            return;
        end
        
        if isempty(GUIinput.fourier_cutfrac.k)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a cut-fraction of wave number');
            uicontrol(h.spatial_edit_cutfrac_kx);
            return;
        end
        
        if isempty(GUIinput.spatialgrid.dx)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a dx grid size');
            uicontrol(h.spatial_edit_dx);
            return;
        end
        
        if isempty(GUIinput.spatialgrid.dy)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a dy grid size');
            uicontrol(h.spatial_edit_dy);
            return;
        end
        
        Nx=round((GUIinput.spatialinterv.xmax-GUIinput.spatialinterv.xmin)*10000/...
            (GUIinput.spatialgrid.dx*10000))+1;
        Ny=round((GUIinput.spatialinterv.ymax-GUIinput.spatialinterv.ymin)*10000/...
            (GUIinput.spatialgrid.dy*10000))+1;

        if any(factor(Nx)>13)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Number of discritization points in x axis has a largest prime factor>13. That may slow down the computation, adjust the dx to avoid that!');
            uicontrol(h.spatial_edit_dx);
            return;
        end
        
        if any(factor(Ny)>13)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Number of discritization points in y axis has a prime factor> 13. That may slow down the computation, adjust the dy to avoid that!');
            uicontrol(h.spatial_edit_dy);
            return;
        end
        
        if strcmpi(GUIinput.spatialwall.check,'yes')
            
            walldata=GUIinput.spatialwall.data;
            Nwall=length(walldata(:,1));
            
            if isempty(walldata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify wall parameter in the table');
                return;
            end
            
            for i=1:Nwall
                if strcmp(walldata(i,1),'Circle')
                    flagwarn=0;
                    if isempty(cell2mat(walldata(i,2)))==1,flagwarn=1;end
                    if isempty(cell2mat(walldata(i,3)))==1,flagwarn=1;end
                    if isempty(cell2mat(walldata(i,4)))==1,flagwarn=1;end
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify wall parameters in the table ');
                        return;
                    end
                    if cell2mat(walldata(i,2))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(walldata(i,2))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(walldata(i,3))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(walldata(i,3))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    if cell2mat(walldata(i,4))<=0
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a radius of the wall > 0 ');
                        return;
                    end
                    
                elseif strcmp(walldata(i,1),'Rectangle')
                    flagwarn=0;
                    if isempty(cell2mat(walldata(i,5)))==1,flagwarn=1;end
                    if isempty(cell2mat(walldata(i,6)))==1,flagwarn=1;end
                    if isempty(cell2mat(walldata(i,7)))==1,flagwarn=1;end
                    if isempty(cell2mat(walldata(i,8)))==1,flagwarn=1;end
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify wall parameters in the table ');
                        return;
                    end
                    
                    if cell2mat(walldata(i,5))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(walldata(i,5))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(walldata(i,6))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(walldata(i,6))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(walldata(i,5))>cell2mat(walldata(i,6))
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Wall position x min> x max. Specify a correct input for x min and x max.');
                        return;
                    end
                    
                    if cell2mat(walldata(i,7))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(walldata(i,7))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    if cell2mat(walldata(i,8))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(walldata(i,8))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(walldata(i,7))>cell2mat(walldata(i,8))
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Wall position y min> y max. Specify a correct input for y min and y max.');
                        return;
                    end
                    
                    
                else
                    walluser=GUIinput.spatialwall.userdata;               
                    if isempty(walluser(i).bdry)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify an user-defined wall boundary ');
                        return;
                    end
%                      assignin('base','walluser',walluser)
                    if min(walluser(i).bdry(:,1))<GUIinput.spatialinterv.xmin || ...
                            min(walluser(i).bdry(:,1))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if max(walluser(i).bdry(:,1))<GUIinput.spatialinterv.xmin || ...
                            max(walluser(i).bdry(:,1))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if min(walluser(i).bdry(:,2))<GUIinput.spatialinterv.ymin || ...
                            min(walluser(i).bdry(:,2))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                    
                    if max(walluser(i).bdry(:,2))<GUIinput.spatialinterv.ymin || ...
                            max(walluser(i).bdry(:,2))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> A wall is at outside the simulation domain ');
                        return;
                    end
                end
                
                if strcmpi(walldata(i,10),'Uniform')
                    reflcoef=cell2mat(walldata(i,11));
                    if isempty(reflcoef)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a reflection coef. in the table ');
                        return;
                    else
                        if (reflcoef>1) || reflcoef<=0
                            set(h.monitorbox,'foregroundcolor','r','string',...
                                '>> Specify a reflection coef. in (0,1] ');
                            return;
                        end
                    end
                else
                    reflEq=cell2mat(walldata(i,12));
                    if isempty(reflEq)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a reflection coef. equation as a function of radial frequency (w) [rad/s] ');
                        return;
                    end
                    
                end
            end
            set(h.monitorbox,'foregroundcolor','k','string',...
                '>> ');
        end
        
        if strcmpi(GUIinput.spatialdamp.check,'yes')
            dampdata=GUIinput.spatialdamp.data;
            Ndamp=length(dampdata(:,1));
            
            for i=1:Ndamp
                if strcmp(dampdata(i,1),'Fourier Bdry.')
                    flagwarn=0;
                    if isempty(cell2mat(dampdata(i,2)))==1,flagwarn=1;end
                    if isempty(cell2mat(dampdata(i,3)))==1,flagwarn=1;end
                    if isempty(cell2mat(dampdata(i,4)))==1,flagwarn=1;end
                    if isempty(cell2mat(dampdata(i,5)))==1,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify damping parameters in the table ');
                        return;
                    end
                else
                    dampuser=GUIinput.spatialdamp.userdata;
                    
                    if isempty(isempty(dampuser(i).bdry))
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify an user-defined damping boundary ');
                        return;
                    end
                    
                    if isempty(cell2mat(dampdata(i,6)))==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a smooth fact. for the userdefined damping ');
                        return;
                    end
                    
                    if  cell2mat(dampdata(i,6))<1 ||mod(cell2mat(dampdata(i,6)),2)~=1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> The smooth fact. of the userdefined damping must be an odd positive integer.');
                        return;
                    end
                    
                end
            end
            set(h.monitorbox,'foregroundcolor','k','string',...
                '>> ');
        end
        
        if strcmpi(GUIinput.spatialbathy.type,'Flat')
            if isempty(GUIinput.spatialbathy.depth)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a water depth ');
                uicontrol(h.bathymetry_edit_depth);
                return;
            end
            
            if GUIinput.spatialbathy.depth<=0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a water depth >0 ');
                uicontrol(h.bathymetry_edit_depth);
                return;
            end
        elseif strcmpi(GUIinput.spatialbathy.type,'Slope (in x-axis)') || ...
                strcmpi(GUIinput.spatialbathy.type,'Slope (in y-axis)')
            
            if isempty(GUIinput.spatialbathy.xymindepth)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a minimum water depth ');
                uicontrol(h.bathymetry_edit_mindepth);
                return;
            end
            
            if GUIinput.spatialbathy.xymindepth<0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a minimum water depth >0 ');
                uicontrol(h.bathymetry_edit_mindepth);
                return;
            end
            
            if isempty(GUIinput.spatialbathy.xymiddepth)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify the depth-reference for 3-depth interpolation');
                uicontrol(h.bathymetry_edit_middepth);
                return;
            end
            if GUIinput.spatialbathy.interp > 2
                if GUIinput.spatialbathy.xymiddepth < min(GUIinput.spatialbathy.xymindepth,GUIinput.spatialbathy.xymaxdepth) || ...
                        GUIinput.spatialbathy.xymiddepth > max(GUIinput.spatialbathy.xymindepth,GUIinput.spatialbathy.xymaxdepth)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> the minimum depth < THE MID-DEPTH REFERENCE < the maximum depth ');
                    uicontrol(h.bathymetry_edit_middepth);
                    return;
                end
            end
            
            if isempty(GUIinput.spatialbathy.xymaxdepth)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a maximum water depth ');
                uicontrol(h.bathymetry_edit_maxdepth);
                return;
            end
            
            if GUIinput.spatialbathy.xymaxdepth<0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a maximum water depth >0 ');
                uicontrol(h.bathymetry_edit_mindepth);
                return;
            end
            
            if isempty(GUIinput.spatialbathy.slope)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a slope of bathymetry');
                uicontrol(h.bathymetry_edit_slope);
                return;
            end
            
            if isempty(GUIinput.spatialbathy.startslope)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a start position of the slope');
                uicontrol(h.bathymetry_edit_startslope);
                return;
            end
            if strcmpi(GUIinput.spatialbathy.type,'Slope (in x-axis)')
                if GUIinput.spatialbathy.startslope<GUIinput.spatialinterv.xmin||...
                        GUIinput.spatialbathy.startslope>GUIinput.spatialinterv.xmax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> A slope is at the ouside simulation domain');
                    uicontrol(h.bathymetry_edit_startslope);
                    return;
                end
            else
                if GUIinput.spatialbathy.startslope<GUIinput.spatialinterv.ymin||...
                        GUIinput.spatialbathy.startslope>GUIinput.spatialinterv.ymax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> A slope is at the ouside simulation domain');
                    uicontrol(h.bathymetry_edit_startslope);
                    return;
                end
                
            end
            
        elseif strcmpi(GUIinput.spatialbathy.type,'Shore (in x-axis)') || ...
                strcmpi(GUIinput.spatialbathy.type,'Shore (in y-axis)')
            if isempty(GUIinput.spatialbathy.maxdepthshore)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a maximum water depth');
                uicontrol(h.bathymetry_edit_maxdepthshore);
                return;
            end
            
            if GUIinput.spatialbathy.maxdepthshore<0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a maximum water depth >0 ');
                uicontrol(h.bathymetry_edit_maxdepthshore);
                return;
            end
            
             if isempty(GUIinput.spatialbathy.mindepthshore)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a minimum water depth (i.e 2% Hs)');
                uicontrol(h.bathymetry_edit_maxdepthshore);
                return;
             end
            
            if GUIinput.spatialbathy.mindepthshore<0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a minimum water depth >0 ');
                uicontrol(h.bathymetry_edit_mindepthshore);
                return;
            end
            
            if isempty(GUIinput.spatialbathy.slopeshore)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a slope of bathymetry');
                uicontrol(h.bathymetry_edit_slopeshore);
                return;
            end
            
            
            if isempty(GUIinput.spatialbathy.shoreposition)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a shore position');
                uicontrol(h.bathymetry_edit_shoreposition);
                return;
            end
            
            if isempty(GUIinput.spatialbathy.xymiddepth)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify the depth-reference for 3-depth interpolation');
                uicontrol(h.bathymetry_edit_middepth);
                return;
            end
            
            
            if strcmpi(GUIinput.spatialbathy.type,'Shore (in x-axis)' )
                if GUIinput.spatialbathy.shoreposition<GUIinput.spatialinterv.xmin||...
                        GUIinput.spatialbathy.shoreposition>GUIinput.spatialinterv.xmax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Shore position is at the ouside simulation domain');
                    uicontrol(h.bathymetry_edit_startslope);
                    return;
                end
            else
                if GUIinput.spatialbathy.shoreposition<GUIinput.spatialinterv.ymin||...
                        GUIinput.spatialbathy.shoreposition>GUIinput.spatialinterv.ymax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>>Shore position is at the ouside simulation domain');
                    uicontrol(h.bathymetry_edit_startslope);
                    return;
                end
                
            end
            
        elseif strcmpi(GUIinput.spatialbathy.type,'User-defined')
             if GUIinput.spatialbathy.interp == 3
                if isempty(GUIinput.spatialbathy.xymiddepth)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify the depth-reference for 3-depth interpolation');
                    uicontrol(h.bathymetry_edit_middepth);
                    return;
                end
             end
            
            if isempty(GUIinput.spatialbathy.userdata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify bathymetry data');
                uicontrol(h.bathymetry_button_load);
                return;
            end
            
            if any(GUIinput.spatialbathy.userdata.data(:,3)>0)
                if isempty(GUIinput.spatialbathy.mindepthshore)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify a minimum water depth (i.e 2% Hs)');
                    uicontrol(h.bathymetry_edit_mindepthshore);
                    return;
                end
                
                if GUIinput.spatialbathy.mindepthshore<0
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify a minimum water depth >0 ');
                    uicontrol(h.bathymetry_edit_mindepthshore);
                    return;
                end
                
            end
%             maxX=max(GUIinput.spatialbathy.userdata.data(:,1));
%             minX=min(GUIinput.spatialbathy.userdata.data(:,1));
%             maxY=max(GUIinput.spatialbathy.userdata.data(:,2));
%             minY=min(GUIinput.spatialbathy.userdata.data(:,2));
            
        end
       
        
        if GUIinput.spatialfriction.check==1
            fricdata=GUIinput.spatialfriction.data;
            Nfrict=length(fricdata(:,1));
            
            if isempty(fricdata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify bottom friction parameter in the table');
                return;
            end
            
            for i=1:Nfrict
                if strcmp(fricdata(i,1),'Rectangle')
                    flagwarn=0;
                    if isempty(cell2mat(fricdata(i,2)))==1,flagwarn=1;end
                    if isempty(cell2mat(fricdata(i,3)))==1,flagwarn=1;end
                    if isempty(cell2mat(fricdata(i,4)))==1,flagwarn=1;end
                    if isempty(cell2mat(fricdata(i,5)))==1,flagwarn=1;end
                    if isempty(cell2mat(fricdata(i,6)))==1,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify bottom friction parameters in the table ');
                        return;
                    end
                    if cell2mat(fricdata(i,2))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(fricdata(i,2))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(fricdata(i,3))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(fricdata(i,3))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary is at outside the simulation domain ');
                        return;
                    end
                    if cell2mat(fricdata(i,2)) > cell2mat(fricdata(i,3))
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary,  specify a correct input for x_min and x_max (x_min<x_max) ');
                        return; 
                    end
                    
                    if cell2mat(fricdata(i,4))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(fricdata(i,4))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(fricdata(i,5))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(fricdata(i,5))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary is at outside the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(fricdata(i,4)) > cell2mat(fricdata(i,5))
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Bottom friction boundary,  specify a correct input for y_min and y_max (y_min<y_max) ');
                        return; 
                    end
                    
                    if cell2mat(fricdata(i,6))<=0
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a coef. of bottom friction > 0 ');
                        return;
                    end
                                
                else
                    frictuser=GUIinput.spatialfriction.userdata;
                   
                    if isempty(frictuser(i).bdry)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify an user-defined bottom friction boundary ');
                        return;
                    end
                    
                    if isempty(cell2mat(fricdata(i,6)))==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a coef. of bottom friction');
                        return;
                    end
                    if cell2mat(fricdata(i,6))<=0
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify a coef. of bottom friction > 0 ');
                        return;
                    end
                end
            end
            set(h.monitorbox,'foregroundcolor','k','string',...
                '>> ');
        end
        
        
        if strcmpi(GUIinput.wave_ivp.type,'Gaussian (along x-axis)')||...
                strcmpi(GUIinput.wave_ivp.type,'Gaussian (along y-axis)')||...
                strcmpi(GUIinput.wave_ivp.type,'N-wave (along x-axis)')||...
                strcmpi(GUIinput.wave_ivp.type,'N-wave (along y-axis)')
            if isempty(GUIinput.wave_ivp.A)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify wave amplitude for the initial condition');
                uicontrol(h.h.waveinput_ivp_edit_A);
                return;
            end
            
            if isempty(GUIinput.wave_ivp.centerposition)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify center position for the initial condition');
                uicontrol(h.waveinput_ivp_edit_centre_position);
                return;
            end
            
            if isempty(GUIinput.wave_ivp.stdev)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify standard deviation for the initial condition');
                uicontrol(h.waveinput_ivp_edit_sd);
                return;
            end
        elseif strcmpi(GUIinput.wave_ivp.type,'User-defined')
            if isempty(GUIinput.wave_ivp.userdata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify user data for the initial condition');
                uicontrol(h.waveinput_ivp_button_load);
                return;
            end
        end
        
        if strcmpi(GUIinput.wave_influx.type,'Yes')
            wavepropdata=GUIinput.wave_influx.propertiesdata;
            influxmethoddata=GUIinput.wave_influx.methoddata;
            
            if isempty(wavepropdata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify wave influx parameter in the table');
                return;
            end
            
            
            Ninflux=length(wavepropdata(:,1));
            
            for i=1:Ninflux
                
                %%%% wave properties
                if strcmp(wavepropdata(i,1),'Harmonic')
                    if iscellstr(wavepropdata(i,2))
                     wavepropdata(i,2)={str2num(cell2mat(wavepropdata(i,2)))};
                    end
                    if iscellstr(wavepropdata(i,4))
                     wavepropdata(i,4)={str2num(cell2mat(wavepropdata(i,4)))};
                    end
                    if iscellstr(wavepropdata(i,5))
                     wavepropdata(i,5)={str2num(cell2mat(wavepropdata(i,5)))};
                    end
                    if iscellstr(wavepropdata(i,6))
                     wavepropdata(i,6)={str2num(cell2mat(wavepropdata(i,6)))};
                    end
                    if iscellstr(wavepropdata(i,7))
                     wavepropdata(i,7)={str2num(cell2mat(wavepropdata(i,7)))};
                    end
                    
                    
                    flagwarn=0;
                    if isempty(cell2mat(wavepropdata(i,2)))==1 ||cell2mat(wavepropdata(i,2))==0
                        flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,4)))==1 || cell2mat(wavepropdata(i,4))==0
                        flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,7)))==1,flagwarn=1;end
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify wave properties in the table ');
                        return;
                    end
                    
                elseif  strcmp(wavepropdata(i,1),'JONSWAP')
                    flagwarn=0;
                    if iscellstr(wavepropdata(i,3))
                     wavepropdata(i,3)={str2num(cell2mat(wavepropdata(i,3)))};
                    end
                    if iscellstr(wavepropdata(i,4))
                     wavepropdata(i,4)={str2num(cell2mat(wavepropdata(i,4)))};
                    end
                    if iscellstr(wavepropdata(i,5))
                     wavepropdata(i,5)={str2num(cell2mat(wavepropdata(i,5)))};
                    end
                    if iscellstr(wavepropdata(i,6))
                     wavepropdata(i,6)={str2num(cell2mat(wavepropdata(i,6)))};
                    end
                    if iscellstr(wavepropdata(i,7))
                     wavepropdata(i,7)={str2num(cell2mat(wavepropdata(i,7)))};
                    end
                    
                    if isempty(cell2mat(wavepropdata(i,3)))==1 || cell2mat(wavepropdata(i,3))==0
                        flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,4)))==1 || cell2mat(wavepropdata(i,4))==0
                        flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,5)))==1,flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,6)))==1,flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,7)))==1,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify wave properties in the table ');
                        return;
                    end
                elseif strcmp(wavepropdata(i,1),'User-defined (signal)')
                    influxuser=GUIinput.wave_influx.userdata;
                    
                    if isempty(influxuser(i).eta)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify an user-defined influx signal(s) ');
                        return;
                    end
                else
                    influxuser=GUIinput.wave_influx.userdata;
                    if isempty(influxuser(i).varianceDensity)
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify an user-defined variance density spectrum');
                        return;
                    end
                    flagwarn=0;
                    
                    
                    if iscellstr(wavepropdata(i,6))
                     wavepropdata(i,6)={str2num(cell2mat(wavepropdata(i,6)))};
                    end
                    if iscellstr(wavepropdata(i,7))
                     wavepropdata(i,7)={str2num(cell2mat(wavepropdata(i,7)))};
                    end
                    
                    if isempty(cell2mat(wavepropdata(i,6)))==1,flagwarn=1;end
                    if isempty(cell2mat(wavepropdata(i,7)))==1,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify wave properties in the table ');
                        return;
                    end
                    
                end
                
                %%%% Influx method and position of influx line
                if strcmpi(influxmethoddata(i,2),'Straight')
                    flagwarn=0;
                    
                  
                    if iscellstr(influxmethoddata(i,3))
                     influxmethoddata(i,3)={str2num(cell2mat(influxmethoddata(i,3)))};
                    end
                    if iscellstr(influxmethoddata(i,4))
                     influxmethoddata(i,4)={str2num(cell2mat(influxmethoddata(i,4)))};
                    end
                    if iscellstr(influxmethoddata(i,5))
                     influxmethoddata(i,5)={str2num(cell2mat(influxmethoddata(i,5)))};
                    end
                    if iscellstr(influxmethoddata(i,6))
                     influxmethoddata(i,6)={str2num(cell2mat(influxmethoddata(i,6)))};
                    end
                    
                    if isempty(cell2mat(influxmethoddata(i,3)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,4)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,5)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,6)))==1,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify the influx line position ');
                        return;
                    end
                    
                     if cell2mat(influxmethoddata(i,3))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(influxmethoddata(i,3))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                     end
                    
                    if cell2mat(influxmethoddata(i,4))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(influxmethoddata(i,4))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(influxmethoddata(i,5))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(influxmethoddata(i,5))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                    end
                    
                      if cell2mat(influxmethoddata(i,6))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(influxmethoddata(i,6))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                    end
 
                    flagwarn=0;
                    
                    if iscellstr(influxmethoddata(i,12))
                     influxmethoddata(i,12)={str2num(cell2mat(influxmethoddata(i,12)))};
                    end
                    if iscellstr(influxmethoddata(i,13))
                     influxmethoddata(i,13)={str2num(cell2mat(influxmethoddata(i,13)))};
                    end
                    if iscellstr(influxmethoddata(i,14))
                     influxmethoddata(i,14)={str2num(cell2mat(influxmethoddata(i,14)))};
                    end
                    
                    if isempty(cell2mat(influxmethoddata(i,12)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,13)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,14)))==1,flagwarn=1;end
                    if cell2mat(influxmethoddata(i,12))>=cell2mat(influxmethoddata(i,13))
                        flagwarn=1;end
                    if cell2mat(influxmethoddata(i,14))==0,flagwarn=1;end
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify the time interval for the coresponding signal ');
                        return;
                    end
                else
                    flagwarn=0;
                    
                    if iscellstr(influxmethoddata(i,7))
                     influxmethoddata(i,7)={str2num(cell2mat(influxmethoddata(i,7)))};
                    end
                    if iscellstr(influxmethoddata(i,8))
                     influxmethoddata(i,8)={str2num(cell2mat(influxmethoddata(i,8)))};
                    end
                    if iscellstr(influxmethoddata(i,9))
                     influxmethoddata(i,9)={str2num(cell2mat(influxmethoddata(i,9)))};
                    end
                    if iscellstr(influxmethoddata(i,10))
                     influxmethoddata(i,10)={str2num(cell2mat(influxmethoddata(i,10)))};
                    end
                    if iscellstr(influxmethoddata(i,11))
                     influxmethoddata(i,11)={str2num(cell2mat(influxmethoddata(i,11)))};
                    end
                    
                    if isempty(cell2mat(influxmethoddata(i,7)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,8)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,9)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,10)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,11)))==1,flagwarn=1;end
                    
                    
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify the influx line position ');
                        return;
                    end
                    
                    if cell2mat(influxmethoddata(i,7))<GUIinput.spatialinterv.xmin || ...
                            cell2mat(influxmethoddata(i,7))>GUIinput.spatialinterv.xmax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                    end
                    
                    if cell2mat(influxmethoddata(i,8))<GUIinput.spatialinterv.ymin || ...
                            cell2mat(influxmethoddata(i,8))>GUIinput.spatialinterv.ymax
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Influx line is in the outside of the simulation domain ');
                        return;
                    end
                    
                   if cell2mat(influxmethoddata(i,10))>cell2mat(influxmethoddata(i,11))
                      set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Theta 1 must be smaller than theta 2 ');
                        return; 
                   end
                    theta1=cell2mat(influxmethoddata(i,10));
                    theta2=cell2mat(influxmethoddata(i,11));
                    
                    Kwad1=findG_kwadran_angle(theta1);
                    Kwad2=findG_kwadran_angle(theta2);
                    flagwarn=0;
                    if Kwad1<=3
                        if Kwad2-Kwad1>1, flagwarn=1; end
                    elseif Kwad1==4
                        if Kwad2>1, flagwarn=1; end
                    end   
                    if theta1==theta2
                      set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Theta 1 cannot be same as theta 2.');
                        return;  
                        
                    end
                    if flagwarn==1
                      set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Theta 2 must be in the same kuadran or one kuadran higher of theta 1.');
                        return; 
                   end
                        
                    flagwarn=0;
                   
                    if iscellstr(influxmethoddata(i,12))
                     influxmethoddata(i,12)={str2num(cell2mat(influxmethoddata(i,12)))};
                    end
                    if iscellstr(influxmethoddata(i,13))
                     influxmethoddata(i,13)={str2num(cell2mat(influxmethoddata(i,13)))};
                    end
                    if iscellstr(influxmethoddata(i,14))
                     influxmethoddata(i,14)={str2num(cell2mat(influxmethoddata(i,14)))};
                    end
                    
                    if isempty(cell2mat(influxmethoddata(i,12)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,13)))==1,flagwarn=1;end
                    if isempty(cell2mat(influxmethoddata(i,14)))==1,flagwarn=1;end
                    if cell2mat(influxmethoddata(i,12))>=cell2mat(influxmethoddata(i,13))
                        flagwarn=1;end
                    if cell2mat(influxmethoddata(i,14))==0,flagwarn=1;end
                             
                    if flagwarn==1
                        set(h.monitorbox,'foregroundcolor','r','string',...
                            '>> Specify the time interval for the coresponding signal ');
                        return;
                    end
                    
                end
            end
            set(h.monitorbox,'foregroundcolor','k','string','>> ');
            
            if GUIinput.wave_influx.rampcheck==1
                if isempty(GUIinput.wave_influx.rampfactor)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify the ramp factor for the time signal');
                    uicontrol(h.waveinput_influx_edit_ramp);
                    return;
                end
                
                if GUIinput.wave_influx.rampfactor<0
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> a ramp factor for the time signal must be larger than 0');
                    uicontrol(h.waveinput_influx_edit_ramp);
                    return;
                end
            end
            
            if GUIinput.wave_influx.ramplinecheck==1
                if isempty(GUIinput.wave_influx.ramplinefactor)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify the ramp factor for the influxing line');
                    uicontrol(h.waveinput_influxline_edit_ramp);
                    return;
                end
                if GUIinput.wave_influx.ramplinefactor<0
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> a ramp factor for the influxing line must be larger than 0');
                    uicontrol(h.waveinput_influxline_edit_ramp);
                    return;
                end
                
                if GUIinput.wave_influx.ramplinefactor>=0.5
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> a ramp factor for the influxing line must be smaller than 0.5');
                    uicontrol(h.waveinput_influxline_edit_ramp);
                    return;
                end
            end
            
            if GUIinput.wave_influx.nonlinadjcheck==1
                if isempty(GUIinput.wave_influx.nonlinadjfactor)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify the influxing adjustment factor for a nonlinear simulation');
                    uicontrol(h.waveinput_influx_edit_nonlinadj);
                    return;
                end
                if GUIinput.wave_influx.nonlinadjfactor<0
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> a nonlinear adjustment factor must be larger than 0');
                    uicontrol(h.waveinput_influxline_edit_ramp);
                    return;
                end
            end
        end
        
        
        if GUIinput.wave_bdy_assim.checkVal==2
            if GUIinput.wave_bdy_assim.shapeCheck==1
                if isempty(GUIinput.wave_bdy_assim.R1)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify an outer radius for bdy. assimilation');
                    uicontrol(h.waveinput_bdy_assim_edit_R1);
                    return;
                end
                
                if isempty(GUIinput.wave_bdy_assim.xc)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify x center for bdy. assimilation');
                    uicontrol(h.waveinput_bdy_assim_edit_xc);
                    return
                end
                if isempty(GUIinput.wave_bdy_assim.yc)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify y center for bdy. assimilation');
                    uicontrol(h.waveinput_bdy_assim_edit_yc);
                    return;
                end
                
                
                if GUIinput.wave_bdy_assim.yc<GUIinput.spatialinterv.ymin || ...
                        GUIinput.wave_bdy_assim.yc>GUIinput.spatialinterv.ymax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> y center for bdy. assimilation is  in outside of simulation domain');
                    uicontrol(h.waveinput_bdy_assim_edit_yc);
                    return;
                end
                
                if GUIinput.wave_bdy_assim.xc<GUIinput.spatialinterv.xmin || ...
                        GUIinput.wave_bdy_assim.xc>GUIinput.spatialinterv.xmax
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> x center for bdy. assimilation is  in outside of simulation domain');
                    uicontrol(h.waveinput_bdy_assim_edit_yc);
                    return;
                end
                
            else
                if isempty(GUIinput.wave_bdy_assim.shapeUserdata)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>>Specify boundary data for assimilation');
                    uicontrol(h.waveinput_bdy_assim_shape_popup);
                    return;
                end
            end
            
            if isempty(GUIinput.wave_bdy_assim.smoothfact)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>>Specify a standard deviation as smoothing factor of boundary assimilation');
                uicontrol(h.waveinput_bdy_assim_edit_smooth);
                return;
            end
            if GUIinput.wave_bdy_assim.smoothfact<0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>>Specify a non-negative standard deviation as smoothing factior of boundary assimilation');
                uicontrol(h.waveinput_bdy_assim_edit_smooth);
                return;
            end
            
            if isempty(GUIinput.wave_bdy_assim.assimdata)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>>Specify data of boundary assimilation');
                uicontrol(h.waveinput_bdy_assim_button_load);
                return;
            end
           if ID_model_phiform==1
            if GUIinput.wave_bdy_assim.cb_phi==1
                if isempty(GUIinput.wave_bdy_assim.assimdata_phi)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>>Specify potential data of boundary assimilation');
                    uicontrol(h.waveinput_bdy_assim_button_load);
                    return;
                end
            end
           else
               if GUIinput.wave_bdy_assim.cb_vel==1
                   if isempty(GUIinput.wave_bdy_assim.assimdata_u)
                       set(h.monitorbox,'foregroundcolor','r','string',...
                           '>>Specify potential data of boundary assimilation');
                       uicontrol(h.waveinput_bdy_assim_button_load_u);
                       return;
                   end
                   if isempty(GUIinput.wave_bdy_assim.assimdata_v)
                       set(h.monitorbox,'foregroundcolor','r','string',...
                           '>>Specify potential data of boundary assimilation');
                       uicontrol(h.waveinput_bdy_assim_button_load_v);
                       return;
                   end
               end
           end
            
            if isempty(GUIinput.wave_bdy_assim.tinit)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify time interval for bdy. assimilation');
                uicontrol(h.waveinput_bdy_assim_time_edit_interval_init);
                return;
            end
            
            if isempty(GUIinput.wave_bdy_assim.tend)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify time interval for bdy. assimilation');
                uicontrol(h.waveinput_bdy_assim_time_edit_interval_init);
                return;
            end
            
             if GUIinput.wave_bdy_assim.tend< GUIinput.wave_bdy_assim.tinit
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a correct time interval for bdy. assimilation ');
                uicontrol(h.waveinput_bdy_assim_time_edit_interval_init);
                return;
            end
            
            if isempty(GUIinput.wave_bdy_assim.dt)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a time step for bdy. assimilation');
                uicontrol(h.waveinput_bdy_assim_time_edit_step);
                return;
            end
            
            set(h.monitorbox,'foregroundcolor','k','string',...
                '>> ');
        end
        
        
        if strcmpi(GUIinput.wave_influx.type,'No') && ...
                strcmpi(GUIinput.wave_ivp.type,'Zero') && ...
                GUIinput.wave_bdy_assim.checkVal==0
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> There is no wave input! Specify the input by an initial condition and/or influxing.');
            return;
        end
        
        if isempty(GUIinput.time.t_start)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify the simulation time interval');
            uicontrol(h.time_edit_interval_init);
            return;
        end
        
        if isempty(GUIinput.time.t_end)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify the simulation time interval');
            uicontrol(h.time_edit_interval_end);
            return;
        end
        
        
        if GUIinput.time.t_start>=GUIinput.time.t_end
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify the simulation time interval');
            uicontrol(h.time_edit_interval_start);
            return;
        end
        
        if isempty(GUIinput.time.dt)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a time step for the simulation');
            uicontrol(h.time_edit_step);
            return;
        end
        
        if GUIinput.time.dt<=0
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Specify a time step for the simulation');
            uicontrol(h.time_edit_step);
            return;
        end
        
        
        if GUIinput.wave_bdy_assim.checkVal==2
            if GUIinput.time.t_start>GUIinput.wave_bdy_assim.tinit
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> The start time of simulation must be less than or same as the start time of assimilation');
                uicontrol(h.time_edit_interval_init);
                return;
            end
            if GUIinput.time.t_end<GUIinput.wave_bdy_assim.tend
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> The end time of simulation must be larger than or same as the end time of assimilation');
                uicontrol(h.time_edit_interval_init);
                return;
            end
            if GUIinput.wave_bdy_assim.dt<GUIinput.time.dt
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> time step of simulation must be smaller than time step of assimilation');
                uicontrol(h.time_edit_interval_init);
                return;
            end
            if mod(GUIinput.wave_bdy_assim.dt,GUIinput.time.dt)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> modulo of time step of assimilation and simulation must be 0');
                uicontrol(h.time_edit_interval_init);
                return;
            end
        end
        
        if GUIinput.option_intflow.check==1
            if isempty(GUIinput.option_intflow.tinit)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify the simulation time interval for the interior flow calculation');
                uicontrol(h.options_edit_interval_init);
                return;
            end
            
            if isempty(GUIinput.option_intflow.tend)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify the simulation time interval for the interior flow calculation');
                uicontrol(h.options_edit_interval_end);
                return;
            end
            
            if isempty(GUIinput.option_intflow.dt_fact)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify a time step  factor for the interior flow calculation');
                uicontrol(h.options_edit_step);
                return;
            end
        end
        
        
        if GUIinput.option_partition.check_def==0
            if isempty(GUIinput.option_partition.total)
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify the number partition of ODE calculation');
                uicontrol(h.options_edit_partition);
                return;
            end
        end
        
        
        if GUIinput.option_mc.check==1
            if  GUIinput.option_mc.waveinput_check+ GUIinput.option_mc.numset_check ...
                    +GUIinput.option_mc.numrun_check==0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>> Specify parameters for monte carlo simulation');
                return
            end
            
            if GUIinput.option_mc.waveinput_check==1
                if isempty(GUIinput.option_mc.waveinput_userdata)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify wave input parameters for monte carlo simulation');
                    return;
                end
            end
            
            if GUIinput.option_mc.numset_check==1
                 if isempty(GUIinput.option_mc.numset_px)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify grid size for monte carlo simulation');
                    uicontrol(h.options_montecarlo_edit_px);
                    return;
                end
                
                
                if isempty(GUIinput.option_mc.numset_py)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify grid size for monte carlo simulation');
                    uicontrol(h.options_montecarlo_edit_py);
                    return;
                end
                
                if length(GUIinput.option_mc.numset_px)~= length(GUIinput.option_mc.numset_py)
                   set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Number of grid size in x must be same as in y');
                    uicontrol(h.options_montecarlo_edit_px);
                    return; 
                end
                
                Nmc=length(GUIinput.option_mc.numset_px);
                
                for jj=1:Nmc
                dxx=GUIinput.option_mc.numset_px(jj);    
                Nx=round((GUIinput.spatialinterv.xmax-GUIinput.spatialinterv.xmin)*10000/...
                    (dxx*10000))+1;
                dyy=GUIinput.option_mc.numset_py(jj);
                Ny=round((GUIinput.spatialinterv.ymax-GUIinput.spatialinterv.ymin)*10000/...
                    (dyy*10000))+1;
                
                if any(factor(Nx)>13)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Number of discritization points in x axis for monte carlo simulation has a largest prime factor>13. That may slow down the computation, adjust the dx to avoid that!');
                    uicontrol(h.spatial_edit_dx);
                    return;
                end
                
                if any(factor(Ny)>13)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Number of discritization points in y axis for monte carlo simulation has a prime factor> 13. That may slow down the computation, adjust the dy to avoid that!');
                    uicontrol(h.spatial_edit_dy);
                    return;
                end
                end
                
            end
            
            if GUIinput.option_mc.numrun_check==1
                if isempty(GUIinput.option_mc_numrun_edit)
                    set(h.monitorbox,'foregroundcolor','r','string',...
                        '>> Specify number of discritization for monte carlo simulation');
                    uicontrol(h.options_montecarlo_edit_Nsimul);
                    return;
                end
            end
        end
        
        
        
        
        set(h.monitorbox,'foregroundcolor','k','string','>>');
        if exist(GUIinput.proj.workdir,'dir')
            if flag_warn_workdir==0
                set(h.monitorbox,'foregroundcolor','r','string',...
                    '>>Warning: Project exists already, it will be overwritten');
                uicontrol(h.project_edit_name);
                flag_warn_workdir=1;
                return
            else
                flag_warn_workdir=0;
                set(h.monitorbox,'foregroundcolor','k','string','>>');
            end
        else
            flag_warn_workdir=0;
            set(h.monitorbox,'foregroundcolor','k','string','>>');
        end
        
        func_save_guistate(GUIinput);
        
        [jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
        set(jProgressBar,'Maximum',100, 'Value',0);
        jProgressBar.setStringPainted( true );
        statusbarObj.setText('pre-processing...');
        
        
        Preproc=MainPreparation(GUIinput,jProgressBar,statusbarObj);
        
        
        input=Preproc(1).input;
        if strcmpi(input.wave.option,'Yes')
            if ~strcmpi(Preproc(1).ivp.typename,'Zero')
                Ninp=input.wave.N+1;
            else
                Ninp=input.wave.N;
            end
            if input.wave.N>1
                Ninp=Ninp+1;
                flagCombineSpect=1;
                if ~strcmpi(Preproc(1).ivp.typename,'Zero')
                    Nend=Ninp-1;
                else Nend=Ninp;
                end
            else
                flagCombineSpect=0;
            end
        else
            Ninp=1;
            flagCombineSpect=0;
        end
        strW=cell(1,Ninp);
        jj=1;
        for ii=1:Ninp
            if flagCombineSpect==1
                
                if ii<Nend
                    strC{jj}=num2str(ii);
                end
                if ii<Nend-1
                    strC{jj+1}='+';
                    jj=jj+2;
                end
                
                if ii==Nend
                    strW(ii)={num2str([strC{1:end}])};
                else
                    strW(ii)={num2str(ii)};
                    if ~strcmpi(Preproc(1).ivp.typename,'Zero') && ii==Ninp
                        strW(ii)={num2str(ii-1)};
                    end
                end
            else
                strW(ii)={num2str(ii)};
            end
        end
        set(h.preview.wave_popup_var,'string',strW,'value',1,'enable','on','userdata',Preproc)
        set(h.preview.wave_popup_spect,'value',1,'enable','on')
        
        statusbarObj.setText('Saving preview images...');
        PreProcView(h,Preproc);
        set(h.preview.spatial_popup_var,'userdata',Preproc);
        
        set(h.menu_run_preproc,'userdata',Preproc);
        set(lmenu.tabgroup ,'selectedtab',lmenu.tabpreview);
        
        set(jProgressBar,'Maximum',100, 'Value',100);
        jProgressBar.setStringPainted( true );
        statusbarObj.setText('Pre-processing done.');
        postproc_flag=1;
    end

    function GUIinput=fun_passing_handles_preproc(h)
        %%%% passing parameter %%%%%%%%%%%%%
        GUIinput.proj.path=h.pathnow;
        GUIinput.proj.name=get(h.project_edit_name,'string');
        GUIinput.proj.note=get(h.project_edit_note,'string');
        GUIinput.proj.projdir=get(h.project_popup_projdir,'string');
        GUIinput.proj.module =h.module;
        GUIinput.proj.workdir=[GUIinput.proj.projdir,'\',GUIinput.proj.name,'\'];
        GUIinput.proj.projhist=h.projechist;
        
        dyntype=cellstr(get(h.model_popup_dynamic,'String'));
        GUIinput.modeldynVal=get(h.model_popup_dynamic,'value');
        GUIinput.modeldyn=dyntype{get(h.model_popup_dynamic,'value')};
        disptype=cellstr(get(h.model_popup_dispersion,'String'));
        GUIinput.modeldispVal=get(h.model_popup_dispersion,'value');
        GUIinput.modeldisp=disptype{get(h.model_popup_dispersion,'value')};
            
        breakcheck=cellstr(get(h.model_popup_breaking,'String'));
        GUIinput.modelbreak.checkVal=get(h.model_popup_breaking,'value');
        GUIinput.modelbreak.check=breakcheck{get(h.model_popup_breaking,'value')};
        GUIinput.modelbreak.initiation=get(h.model_break_edit_initiation,'userdata');
        GUIinput.modelbreak.termination=get(h.model_break_edit_termination,'userdata');
        GUIinput.modelbreak.Tstar=get(h.model_break_edit_Tstar,'userdata');
        
        GUIinput.modelcurrent.check=get(h.model_current_cb,'value');
        GUIinput.modelcurrent.ux=get(h.model_current_ux_edit,'userdata');
        GUIinput.modelcurrent.uy=get(h.model_current_uy_edit,'userdata');
        
        
        GUIinput.spatialinterv.xmin=get(h.spatial_edit_xmin,'userdata');
        GUIinput.spatialinterv.xmax=get(h.spatial_edit_xmax,'userdata');
        GUIinput.spatialinterv.ymin=get(h.spatial_edit_ymin,'userdata');
        GUIinput.spatialinterv.ymax=get(h.spatial_edit_ymax,'userdata');
        
        GUIinput.spatialgrid.dx=get(h.spatial_edit_dx,'userdata');
        GUIinput.spatialgrid.dy=get(h.spatial_edit_dy,'userdata');
        
        GUIinput.fourier_cutfrac.k=get(h.cutfrac_k,'userdata');
        
        wallcheck=cellstr(get(h.wall_popup,'string'));
        GUIinput.spatialwall.checkVal=get(h.wall_popup,'value');
        GUIinput.spatialwall.check=wallcheck{get(h.wall_popup,'value')};
        GUIinput.spatialwall.data=get(h.wall_table,'data');
        GUIinput.spatialwall.userdata=get(h.wall_table,'userdata');
        dampcheck=cellstr(get(h.damping_popup,'string'));
        GUIinput.spatialdamp.checkVal=get(h.damping_popup,'value');
        GUIinput.spatialdamp.check=dampcheck{get(h.damping_popup,'value')};
        GUIinput.spatialdamp.data=get(h.damping_table,'data');
        GUIinput.spatialdamp.userdata=get(h.damping_table,'userdata');
        
        bathtype=cellstr(get(h.bathymetry_popup_type,'string'));
        GUIinput.spatialbathy.typeVal=get(h.bathymetry_popup_type,'value');
        GUIinput.spatialbathy.type=bathtype{get(h.bathymetry_popup_type,'value')};
        GUIinput.spatialbathy.depth=get(h.bathymetry_edit_depth,'userdata');
        GUIinput.spatialbathy.xymindepth=get(h.bathymetry_edit_mindepth,'userdata');
        GUIinput.spatialbathy.xymaxdepth=get(h.bathymetry_edit_maxdepth,'userdata');
        GUIinput.spatialbathy.slope=get(h.bathymetry_edit_slope,'userdata');
        GUIinput.spatialbathy.startslope=get(h.bathymetry_edit_startslope,'userdata');
        GUIinput.spatialbathy.mindepthshore=get(h.bathymetry_edit_mindepthshore,'userdata');
        GUIinput.spatialbathy.maxdepthshore=get(h.bathymetry_edit_maxdepthshore,'userdata');
        GUIinput.spatialbathy.slopeshore=get(h.bathymetry_edit_slopeshore,'userdata');
        GUIinput.spatialbathy.shoreposition=get(h.bathymetry_edit_shoreposition,'userdata');
        GUIinput.spatialbathy.userdata=get(h.bathymetry_button_load ,'userdata');
        GUIinput.spatialbathy.interp  =1+get(h.bathymetry_popup_interpdepth,'value');
        GUIinput.spatialbathy.xymiddepth=get(h.bathymetry_edit_middepth,'userdata'); % nunu
        GUIinput.spatialfriction.check=get(h.bathymetry_checkbox_friction,'value');
        GUIinput.spatialfriction.data=get(h.bathymetry_table_friction,'data');
        GUIinput.spatialfriction.userdata=get(h.bathymetry_table_friction,'userdata');
        
        
        ivptype=cellstr(get(h.waveinput_ivp_popup_type,'string'));
        GUIinput.wave_ivp.typeVal=get(h.waveinput_ivp_popup_type,'value');
        GUIinput.wave_ivp.type=ivptype{get(h.waveinput_ivp_popup_type,'value')};
        GUIinput.wave_ivp.A    =get(h.waveinput_ivp_edit_A,'userdata');
        GUIinput.wave_ivp.centerposition=get(h.waveinput_ivp_edit_centre_position,'userdata');
        GUIinput.wave_ivp.stdev=get(h.waveinput_ivp_edit_sd,'userdata');
        GUIinput.wave_ivp.userdata=get(h.waveinput_ivp_button_load,'userdata');
        
        influxtype=cellstr(get(h.waveinput_influx_popup,'string'));
        GUIinput.wave_influx.typeVal=get(h.waveinput_influx_popup,'value');
        GUIinput.wave_influx.type=influxtype{get(h.waveinput_influx_popup,'value')};
        GUIinput.wave_influx.propertiesdata=get(h.waveinput_influx_proptable,'data');
        GUIinput.wave_influx.userdata=get(h.waveinput_influx_proptable,'userdata');
        GUIinput.wave_influx.methoddata=get(h.waveinput_influx_methodtable,'data');
        GUIinput.wave_influx.rampcheck=get(h.waveinput_influx_checkbox_ramp,'value');
        GUIinput.wave_influx.rampfactor=get(h.waveinput_influx_edit_ramp,'userdata');
        GUIinput.wave_influx.ramplinecheck=get(h.waveinput_influxline_checkbox_ramp,'value');
        GUIinput.wave_influx.ramplinefactor=get(h.waveinput_influxline_edit_ramp,'userdata');
        GUIinput.wave_influx.nonlinadjcheck=get(h.waveinput_influx_checkbox_nonlinadj,'value');
        GUIinput.wave_influx.nonlinadjfactor=get(h.waveinput_influx_edit_nonlinadj,'userdata');
        
        GUIinput.wave_bdy_assim.checkVal=get(h.waveinput_bdy_assim_popup,'value');
        GUIinput.wave_bdy_assim.shapeCheck=get(h.waveinput_bdy_assim_shape_popup,'value');
        GUIinput.wave_bdy_assim.shapeUserdata=get(h.waveinput_bdy_assim_shape_popup,'userdata');
        GUIinput.wave_bdy_assim.R1=get(h.waveinput_bdy_assim_edit_R1,'userdata');
        GUIinput.wave_bdy_assim.xc=get(h.waveinput_bdy_assim_edit_xc,'userdata');
        GUIinput.wave_bdy_assim.yc=get(h.waveinput_bdy_assim_edit_yc,'userdata');
        GUIinput.wave_bdy_assim.smoothfact=get(h.waveinput_bdy_assim_edit_smooth,'userdata');
        GUIinput.wave_bdy_assim.assimdata=get(h.waveinput_bdy_assim_button_load,'userdata');
        if ID_model_phiform==1
        GUIinput.wave_bdy_assim.assimdata_phi=get(h.waveinput_bdy_assim_button_load_phi,'userdata');
        GUIinput.wave_bdy_assim.cb_phi=get(h.waveinput_bdy_assim_cb_phi ,'value');
        else
        GUIinput.wave_bdy_assim.assimdata_u=get(h.waveinput_bdy_assim_button_load_u,'userdata');
        GUIinput.wave_bdy_assim.assimdata_v=get(h.waveinput_bdy_assim_button_load_v,'userdata');
        GUIinput.wave_bdy_assim.cb_vel=get(h.waveinput_bdy_assim_cb_vel ,'value');    
        end
        GUIinput.wave_bdy_assim.propdir=get(h.waveinput_bdy_assim_propdir_popup,'value');
        GUIinput.wave_bdy_assim.cb_nonlinAdj=get(h.waveinput_bdy_assim_cb_nonlinAdj,'value');
        GUIinput.wave_bdy_assim.nonlinAdj_distance=get(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'userdata');
        GUIinput.wave_bdy_assim.nonlinAdj_smooth=get(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'userdata');
        
        
        
        GUIinput.wave_bdy_assim.tinit=get(h.waveinput_bdy_assim_time_edit_interval_init,'userdata');
        GUIinput.wave_bdy_assim.tend=get(h.waveinput_bdy_assim_time_edit_interval_end,'userdata');
        GUIinput.wave_bdy_assim.dt=get(h.waveinput_bdy_assim_time_edit_step,'userdata');
        
        GUIinput.time.t_start=get(h.time_edit_interval_init,'userdata');
        GUIinput.time.t_end=get(h.time_edit_interval_end,'userdata');
        GUIinput.time.dt=get(h.time_edit_step,'userdata');
        
        GUIinput.option_outputvar=get(h.options_outputvar_popup,'value');
        GUIinput.option_intflow.check=get(h.options_checkbox_intflow,'value');
        GUIinput.option_intflow.tinit=get(h.options_edit_interval_init,'userdata');
        GUIinput.option_intflow.tend=get(h.options_edit_interval_end,'userdata');
        GUIinput.option_intflow.dt_fact=get(h.options_edit_step,'userdata');
        
        GUIinput.option_partition.check_def=get(h.options_checkbox_default,'value');
        GUIinput.option_partition.total=get(h.options_edit_partition,'userdata');
        GUIinput.option_partition.check_combine=get(h.options_checkbox_combine,'value');
        GUIinput.option_ode_tol=get(h.options_edit_ode_tol,'userdata');
        GUIinput.option_ode_sol=get(h.options_odesolv_popup,'value');
        
        GUIinput.option_mc.check=get(h.options_montecarlo_cb_check,'value');
        GUIinput.option_mc.waveinput_check=get(h.options_montecarlo_cb_waveinput,'value');
        GUIinput.option_mc.waveinput_userdata=get(h.options_montecarlo_pb_waveinput,'userdata');
        GUIinput.option_mc.numset_check=get(h.options_montecarlo_cb_numset,'value');
        GUIinput.option_mc.numset_px=get(h.options_montecarlo_edit_px,'userdata');
        GUIinput.option_mc.numset_py=get(h.options_montecarlo_edit_py,'userdata');
        GUIinput.option_mc.numrun_check=get(h.options_montecarlo_cb_Nsimul,'value');
        GUIinput.option_mc_numrun_edit=str2num(get(h.options_montecarlo_edit_Nsimul,'string'));
        
    end

    function callback_run_odesolver(hObj,eventdata,h)
        global postproc_flag postprocIF_flag
        set(h.monitorbox,'foregroundcolor','k','string','>>');
        if postproc_flag==1
            preproc=get(h.menu_run_preproc,'userdata');
            set(h.menu_run_preproc,'userdata',[]);%clear memory
            set(h.waveinput_bdy_assim_button_load,'userdata',[]);%clear memory
            if ID_model_phiform==1
            set(h.waveinput_bdy_assim_button_load_phi,'userdata',[]);%clear memory
            else
            set(h.waveinput_bdy_assim_button_load_u,'userdata',[]);%clear memory
            set(h.waveinput_bdy_assim_button_load_v,'userdata',[]);%clear memory
            end
            %[output,log]=ODEsolver_uv(h,preproc); %%for runup UV by Nida
            [output,outkinematic,log]=ODEsolver(h,preproc);
            set(h.preview.log,'string',{log.string{:}})
            set(h.menu_run_simul,'userdata',output);
            set(h.waveinput_bdy_assim_edit_data,'string','');
            my_data.Proj=preproc(end).Proj;
            my_data.bath=preproc(end).bath;
            my_data.dom=preproc(end).dom;
            my_data.influx=preproc(end).influx;
            preproc(end).input.bdyassim.assimdata=[];%clear memory;
            preproc(end).input.bdyassim.assimdata_phi=[];%clear memory;
            my_data.input=preproc(end).input;
            my_data.model=preproc(end).model;
            my_data.par=preproc(end).par;
            my_data.ivp=preproc(end).ivp;
            my_data.bdyassim=preproc(end).bdyassim;
            my_data.timeSimul=preproc(end).timeSimul;
            options=preproc(end).options;
            
            if options.interior.check==1
                reset_handles_postprocIF(h);
                set(h.IF_project_popup_projdir,'string',preproc(end).Proj.projdir);
                
                if ~isempty(outkinematic)
                    set(h.IF_project_edit_name,'string',preproc(end).Proj.savename);
                    set(h.IF_project_edit_note,'string',my_data.Proj.usernote);
                     if options.mc.check==0
                        set(h.IF_project_edit_data,'string',[preproc(end).Proj.savename,'_simul_Interior2D','.mat loaded']);
                    else
                        set(h.IF_project_edit_data,'string',[preproc(end).Proj.savename,'_simul_Interior2D_mc_',num2str(length(preproc)),'.mat loaded']);
                    end
                end
               
                my_data_temp=my_data;
                my_data_temp.outkinematic=outkinematic;
                set(h.IF_project_load_data,'userdata',my_data_temp);
                if ~isempty(outkinematic)
                input_handles_postprocIF(h,my_data_temp);
                end
                postprocIF_flag=0;
               % set(lmenu.tabgroup ,'selectedtab',lmenu.tabpostproc);
                
            end
            
             
            my_data.output=output;
            %%%%Passing variables to Postproc panel%%%%%%%%%%%%%%%%%%%%%%%%%%%
            reset_handles_postproc(h);
            set(h.pp_project_popup_projdir,'string',preproc(end).Proj.projdir);
            set(h.pp_project_edit_name,'string',preproc(end).Proj.savename);
            
            if options.mc.check==0
                set(h.pp_project_edit_data,'string',[preproc(end).Proj.savename,'_simul','.mat loaded']);
            else
                set(h.pp_project_edit_data,'string',[preproc(end).Proj.savename,'_simul_mc_',num2str(length(preproc)),'.mat loaded']);
            end
            set(h.pp_project_edit_note,'string',my_data.Proj.usernote);
            
            set(h.pp_project_load_data,'userdata',my_data);
            input_handles_postproc(h,my_data);
            postproc_flag=0;
            set(lmenu.tabgroup ,'selectedtab',lmenu.tabpostproc);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if my_data.model.phiForm==1
             set(h.pp_density_signal_popup_var,'String',{'elevation','potential'});   
             set(h.pp_line_buoy_popup_var,'String',{'elevation','potential'}); 
             set(h.pp_line_MTAA_popup_var,'String',{'elevation','potential'}); 
             set(h.pp_anim_line_popup_var,'String',{'elevation','potential'});
             set(h.pp_validation_buoy_popup_var,'String',{'elevation','potential'})
            else
             set(h.pp_density_signal_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})    
             set(h.pp_line_buoy_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})    
             set(h.pp_line_MTAA_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})    
             set(h.pp_anim_line_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
             set(h.pp_validation_buoy_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
            end
            
        else
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>> Please run pre-processing!');
            return;
        end
    end

    function callback_pp_projectname(hObj,eventdata,h)
        projname=get(hObj,'String');
        projdir=get(h.pp_project_popup_projdir,'string');
        
        if isempty(projdir)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Specify a working directory');
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if exist([projdir,'/',projname],'dir')
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Warning: Project exists already, it will be overwritten');
            uicontrol(h.pp_project_edit_name);
            flag_warn_workdir=1;
        else
            flag_warn_workdir=0;
            set(h.monitorbox,'foregroundcolor','k','string','>>');
        end
    end

    function callback_pp_project_button_projdir(hObj,eventdata,h)
        pathnow=h.pathnow;
        workingdir = uigetdir(pathnow,'browse a directory');
        if workingdir ~= 0
            set(h.pp_project_popup_projdir,'String',workingdir);
        else
            set(h.monitorbox,'String',['>>Warning: No directory is loaded'],'foregroundcolor','k');
            uicontrol(h.pp_project_button_projdir)
            return;
        end
        
    end



    function callback_pp_loaddata(hObj,eventdata,h)
        projdir=get(h.pp_project_popup_projdir,'string');
        
        if isempty(projdir)||strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify a project directory for post-processing');
            return;
        end
        
        
        
        [file_name,directory]=uigetfile([projdir,'\','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
        
        if directory~=0
            set(h.monitorbox,'foregroundcolor','k','string','>>loading...');
            pause(0.001)
            my_data=load([directory,file_name]);
            my_data.Load=1;
            if isempty(my_data)
                set(h.pp_project_edit_data,'string','');
                set(h.monitorbox,'foregroundcolor','r','string','>>No data loaded');
            elseif ~isfield(my_data,'dom')||~isfield(my_data,'output')
                set(h.pp_project_edit_data,'string','')
                set(h.monitorbox,'foregroundcolor','r','string','>>Wrong input file');
            else
                reset_handles_postproc(h);
                set(h.pp_project_popup_projdir,'string',projdir);
                set(h.pp_project_edit_name,'string',my_data.Proj.savename);
                set(h.pp_project_edit_data,'string',file_name);
                set(h.pp_project_edit_note,'string',my_data.Proj.usernote);
                set(h.monitorbox,'foregroundcolor','k','string','>>Data loaded');
                 if isfield(my_data.model,'phiForm')
                    if my_data.model.phiForm==1
                        set(h.pp_density_signal_popup_var,'String',{'elevation','potential'});
                        set(h.pp_line_buoy_popup_var,'String',{'elevation','potential'});
                        set(h.pp_line_MTAA_popup_var,'String',{'elevation','potential'});
                        set(h.pp_anim_line_popup_var,'String',{'elevation','potential'});
                        set(h.pp_validation_buoy_popup_var,'String',{'elevation','potential'})
                    else
                        set(h.pp_density_signal_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
                        set(h.pp_line_buoy_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
                        set(h.pp_line_MTAA_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
                        set(h.pp_anim_line_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
                        set(h.pp_validation_buoy_popup_var,'String',{'elevation','velocity (in x dir)','velocity (in y dir)'})
                    end
                else
                   my_data.model.phiForm=1; 
                end
                set(h.pp_project_load_data,'userdata',my_data);
                input_handles_postproc(h,my_data);
            end
            
        end
    end

    function input_handles_postproc(h,my_data)
        fbl=my_data.dom.fbl;
        
        set(h.pp_density_profile_setting_cb_xlim,'value',1);
        set(h.pp_density_profile_setting_cb_ylim,'value',1);
        set(h.pp_density_profile_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_density_profile_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_density_profile_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_density_profile_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        
        set(h.pp_density_profile_cb_level_eta,'enable','on')
        set(h.pp_density_profile_cb_level_phi,'enable','on')
        set(h.pp_density_profile_cb_level_extremeCrest,'enable','on')
        set(h.pp_density_profile_cb_level_quiver,'enable','off','value',0)
        
        set(h.pp_density_signal_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_density_signal_setting_cb_tlim,'value',1);
        
        
        set(h.pp_density_statistic_setting_cb_xlim,'value',1);
        set(h.pp_density_statistic_setting_cb_ylim,'value',1);
        set(h.pp_density_statistic_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_density_statistic_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_density_statistic_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_density_statistic_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        set(h.pp_density_statistic_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_density_statistic_setting_cb_tlim,'value',1);
        
        set(h.pp_line_buoy_setting_cb_horznlim,'value',1);
        set(h.pp_line_buoy_setting_edit_horznlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_line_statistic_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_line_statistic_setting_cb_tlim,'value',1);
        
        
        set(h.pp_line_ham_mom_setting_cb_xlim,'value',1);
        set(h.pp_line_ham_mom_setting_cb_ylim,'value',1);
        set(h.pp_line_ham_mom_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_line_ham_mom_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_line_ham_mom_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_line_ham_mom_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        set(h.pp_line_ham_mom_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_line_ham_mom_setting_cb_tlim,'value',1);
        
        set(h.pp_line_MTAA_setting_cb_xlim,'value',1);
        set(h.pp_line_MTAA_setting_cb_ylim,'value',1);
        set(h.pp_line_MTAA_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_line_MTAA_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_line_MTAA_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_line_MTAA_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        set(h.pp_line_MTAA_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_line_MTAA_setting_cb_tlim,'value',1);
        
        
        set(h.pp_line_Extreme_setting_cb_xlim,'value',1);
        set(h.pp_line_Extreme_setting_cb_ylim,'value',1);
        set(h.pp_line_Extreme_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_line_Extreme_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_line_Extreme_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_line_Extreme_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        set(h.pp_line_Extreme_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_line_Extreme_setting_cb_tlim,'value',1);
        
        
        set(h.pp_line_breaking_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_line_breaking_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_line_breaking_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_line_breaking_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        if ~isempty(my_data.output.break_nodes)
            set(h.pp_line_breaking_setting_edit_tlim,'userdata',...
                [my_data.output.break_nodes(1,1);my_data.output.break_nodes(end,1)],'enable','on',...
                'string',[num2str(my_data.output.break_nodes(1,1)),';',num2str(my_data.output.break_nodes(end,1))]);
        end
        set(h.pp_line_breaking_setting_cb_tlim,'value',1);
        callback_line_breaking_popup_var(h.pp_line_breaking_popup_var,[],h);
        set(h.pp_line_breaking_setting_cb_xlim,'value',1);
        set(h.pp_line_breaking_setting_cb_ylim,'value',1);
        set(h.pp_line_breaking_setting_edit_xlim,'enable','on');
        set(h.pp_line_breaking_setting_edit_ylim,'enable','on');
        set(h.pp_line_breaking_setting_cb_coarse,'value',0,'enable','off');
        set(h.pp_line_breaking_setting_edit_coarse,'enable','off')
        
        
        set(h.pp_quant_buoy_tinterv_edit,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        
        
        set(h.pp_anim_density_setting_cb_xlim,'value',1);
        set(h.pp_anim_density_setting_cb_ylim,'value',1);
        set(h.pp_anim_density_setting_edit_xlim,'string',...
            [num2str(my_data.dom.X(1)+fbl.l),';',num2str(my_data.dom.X(end)-fbl.r)]);
        set(h.pp_anim_density_setting_edit_xlim,'userdata',...
            [my_data.dom.X(1)+fbl.l;my_data.dom.X(end)-fbl.r],'enable','on');
        
        set(h.pp_anim_density_setting_edit_ylim,'string',...
            [num2str(my_data.dom.Y(1)+fbl.b),';',num2str(my_data.dom.Y(end)-fbl.t)]);
        set(h.pp_anim_density_setting_edit_ylim,'userdata',...
            [my_data.dom.Y(1)+fbl.b;my_data.dom.Y(end)-fbl.t],'enable','on');
        
        
        
        set(h.pp_anim_density_cb_level_eta,'enable','on')
        set(h.pp_anim_density_cb_level_phi,'enable','on')
        set(h.pp_anim_density_cb_level_quiver,'enable','off','value',0)
        
        set(h.pp_anim_density_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_anim_density_setting_cb_tlim,'value',1);
        
        set(h.pp_anim_line_setting_edit_tlim,'userdata',...
            [my_data.output.time(1);my_data.output.time(end)],'enable','on',...
            'string',[num2str(my_data.output.time(1)),';',num2str(my_data.output.time(end))]);
        set(h.pp_anim_line_setting_cb_tlim,'value',1);
        
        
        
    end

    function callback_pp_setting_cb_view(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
            set(href,'userdata',2,'string','2');
        end
    end

    function callback_pp_setting_edit_view(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)>2||length(param)<1
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_clim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end


    function callback_pp_setting_edit_clim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if isempty(param)||length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify colorbar axes limit [zmin;zmax]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_zlim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_zlim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if isempty(param)||length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify vertical axes limit [zmin;zmax]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_xlim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_xlim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_ylim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_ylim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
            uicontrol(hObj);
            return;
        end
    end


    function callback_pp_setting_cb_spatlim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_spatlim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify spatial axes limit [xmin;xmax] or [ymin;ymax]'])
            uicontrol(hObj);
            return;
        end
    end



    function callback_pp_setting_cb_tlim(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_horzlim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>']) 
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify horizontal axes limit [min;max]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_edit_tlim(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_saveanim(hObj,eventdata,href,href1)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
            set(href1,'enable','on');
            set(href,'string','0.001;inf','userdata',str2num('0.001;inf'));
        else
            set(href,'enable','off');
            set(href1,'enable','off');
        end
    end

    function callback_pp_setting_edit_gifset(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=2
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify gif parameters:  [delay factor; number of loop]'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_setting_cb_coarse(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_coarse(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
            uicontrol(hObj);
            return; 
        end
    end

    function callback_pp_setting_cb_ampliref(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_ampliref(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
            
        if length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude ref.'])
            uicontrol(hObj);
            return;
        elseif param<=0
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude ref. that is larger than zero'])
            uicontrol(hObj);
            return;      
        end
    end


    function callback_pp_setting_cb_level(hObj,eventdata,h,href,href1,href2,href3,href4,href5)
        Id=get(hObj,'value');
        if Id==1
            Id1=get(href1,'value');
            Id2=get(href2,'value');
            Id3=get(href3,'value');
            if Id1+Id2+Id3>0
                set(href,'enable','on');
            else
                set(href,'enable','off');
                set(hObj,'value',0);
            end
            Id4=get(href4,'value');
            Id5=get(href5,'value');
            if Id4==1 || Id5==1
                if Id1+Id2+Id3==0
                    set(href,'enable','off');
                    set(hObj,'value',0);
                end
            end
        else
            set(href,'enable','off');
        end
    end

    function callback_pp_setting_edit_level(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
    end

    function callback_pp_setting_cb_savefig(hObj,eventdata,href)
        Id=get(hObj,'value');
        if Id==1
            set(href,'enable','on');
        else
            set(href,'enable','off');
        end
    end



    function callback_density_profile_at_time(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        if length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify a time')
            uicontrol(hObj);
        end
    end

    function callback_density_profile_level_eta(hObj,eventdata,h,href1,href2,href3,href4)
        Id=get(hObj,'value');
        if Id==1
            set(href1,'value',0);set(href2,'value',0);
        end
        
        set(href3,'value',0);
        set(href4,'enable','off');
    end

    function callback_density_profile_level_phi(hObj,eventdata,h,href1,href2,href3,href4)
        Id=get(hObj,'value');
        simuldata=get(h.pp_project_load_data,'userdata');
        if Id==1
            set(href1,'value',0);set(href2,'value',0);
            if isfield(simuldata,'output')
                if ID_model_phiform==1
                    if ~isfield(simuldata.output,'phi')
                        set(hObj,'value',0);
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, only the elevation level-line is available'])
                        set(href1,'value',1);set(href2,'value',0);set(href3,'value',0);
                    end
                    if length(simuldata.output.phi(1,:,:))==1
                        set(hObj,'value',0);
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, only the elevation level-line is available'])
                        set(href1,'value',1);set(href2,'value',0);set(href3,'value',0);
                    end
                else
                    if ~isfield(simuldata.output,'u')
                        set(hObj,'value',0);
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, only the elevation level-line is available'])
                        set(href1,'value',1);set(href2,'value',0);set(href3,'value',0);
                    end
                    if length(simuldata.output.u(1,:,:))==1
                        set(hObj,'value',0);
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, only the elevation level-line is available'])
                        set(href1,'value',1);set(href2,'value',0);set(href3,'value',0);
                    end
                    
                end
            end
        end
        
        set(href3,'value',0);
        set(href4,'enable','off');
    end

    function callback_density_profile_quiver(hObj,eventdata,href1,href2,href3,href4)
        Id=get(hObj,'value');
        if Id==1
            set(href1,'value',0);set(href2,'value',0);
        end
        set(href3,'value',0);
        set(href4,'enable','off');
    end

    function callback_density_profile_level_extremecrest(hObj,eventdata,href1,href2,href3,href4)
        Id=get(hObj,'value');
        set(href1,'value',0);
        set(href2,'enable','off');
        if Id==1
            set(href3,'enable','on');
            Idd=get(href3,'value');
            if Idd==1
                set(href4,'enable','on');
            else
                set(href4,'enable','off');
            end
        else
            set(href3,'enable','off');
            set(href4,'enable','off');
        end
    end

    function callback_density_profile_popup_var(hObj,eventdata,h,href1,href2,href3,href4,href5)
        Id=get(hObj,'value');
        simuldata=get(h.pp_project_load_data,'userdata');
        
        if isempty(simuldata)
            set(hObj,'value',1);
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a simulation data'])
            return;
        else
            if ID_model_phiform==1
                 if ~isfield(simuldata.output,'phi') &&Id>2
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, please specify in the option for output variables and then re-run the simulation.'])
                    return;
                end
                if length(simuldata.output.phi(:,1,1))==1 &&Id>2
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, please specify in the option for output variables and then re-run the simulation.'])
                    return;
                end
            else
                if ~isfield(simuldata.output,'u') &&Id>2
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option for output variables and then re-run the simulation.'])
                    return;
                end
                if length(simuldata.output.u(:,1,1))==1 &&Id>2
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option for output variables and then re-run the simulation.'])
                    return;
                end
            end
        end
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        if Id==1
            set(href1,'enable','on')
            set(href2,'enable','on')
            set(href3,'enable','off','value',0)
            set(href4,'enable','on')
            set(href5,'string','level opt.')
        elseif Id==2
            set(href1,'enable','on')
            set(href2,'enable','on')
            set(href3,'enable','off','value',0)
            set(href4,'enable','on')
            set(href5,'string','level opt.')
        elseif Id==3
            set(href1,'enable','on')
            set(href2,'enable','on')
            set(href3,'enable','off','value',0)
            set(href4,'enable','on')
            set(href5,'string','level opt.')
        else
            set(href1,'enable','off','value',0)
            set(href2,'enable','off','value',0)
            set(href3,'enable','on')
            set(href4,'enable','on')
            set(href5,'string','quiver opt.')
        end
    end

    function callback_pp_plot_density_profile(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_density_profile_popup_var,'value');
        setting.view.check=get(h.pp_density_profile_setting_cb_view,'value');
        setting.view.param=get(h.pp_density_profile_setting_edit_view,'userdata');
        setting.clim.check=get(h.pp_density_profile_setting_cb_clim,'value');
        setting.clim.param=get(h.pp_density_profile_setting_edit_clim,'userdata');
        setting.xlim.check=get(h.pp_density_profile_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_density_profile_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_density_profile_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_density_profile_setting_edit_ylim,'userdata');
        setting.coarse.check=get(h.pp_density_profile_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_density_profile_setting_edit_coarse,'userdata');
        setting.ampliref.check=get(h.pp_density_profile_setting_cb_ampliref,'value');
        setting.ampliref.param=get(h.pp_density_profile_setting_edit_ampliref,'userdata');
        setting.levelquiver.check=get(h.pp_density_profile_setting_cb_level,'value');
        setting.levelquiver.param=get(h.pp_density_profile_setting_edit_level,'userdata');
        setting.savefig.check=get(h.pp_density_profile_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_density_profile_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_density_profile_setting_popup_savefig,'value'));
        setting.level.eta_check=get(h.pp_density_profile_cb_level_eta,'value');
        setting.level.phi_check=get(h.pp_density_profile_cb_level_phi,'value');
        setting.level.quiver_check=get(h.pp_density_profile_cb_level_quiver,'value');
        setting.level.extremecrest=get(h.pp_density_profile_cb_level_extremeCrest,'value');
        setting.level.extremetrough=get(h.pp_density_profile_cb_level_extremeTrough,'value');
        setting.bathy.cb=get(h.pp_density_profile_cb_bathy,'value');
        setting.bathy.scale=get(h.pp_density_profile_edit_bathyScale,'userdata');

        
        colorfig=cellstr(get(h.pp_density_profile_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.pp_density_profile_setting_popup_colormap,'value'));
        axesfig=h.pp_density_profile_axes;
        tsnap=get(h.pp_density_profile_time_edit,'userdata');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if setting.var>2 || setting.level.phi_check==1
            if ID_model_phiform==1
                if ~isfield(simuldata.output,'phi')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_density_profile_popup_var);
                    return;
                end
                if length(simuldata.output.phi(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_density_profile_popup_var);
                    return;
                end
            else
                 if ~isfield(simuldata.output,'u')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_density_profile_popup_var);
                    return;
                end
                if length(simuldata.output.u(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_density_profile_popup_var);
                    return;
                end
            end
        end
        
       
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if isempty(tsnap)||length(tsnap)~=1
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a time')
            uicontrol(h.pp_density_profile_time_edit);
            return;
        end
        
        if  tsnap<T(1) || tsnap>T(end)
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in the interval ['...
                ,num2str(T(1)),';',num2str(T(end)),'] (s)'])
            uicontrol(h.pp_density_profile_time_edit);
            return;
        end
        
        if setting.view.check==1
            if isempty(setting.view.param)||length(setting.view.param)>2||length(setting.view.param)<1
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_profile_setting_edit_view);
                return;
            end
            if length(setting.view.param)==1 && (setting.view.param>3||setting.view.param<2)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_profile_setting_edit_view);
                return;
            end
        end
        if setting.clim.check==1
            if isempty(setting.clim.param)||length(setting.clim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify colorbar axes limit [zmin;zmax]'])
                uicontrol(h.pp_density_profile_setting_edit_clim);
                return;
            end
        end
        if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_density_profile_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_density_profile_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_density_profile_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_density_profile_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_density_profile_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_density_profile_setting_edit_ylim);
                return; 
            end
        end
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_density_profile_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_density_profile_setting_edit_coarse);
                return;
            end
        end
        
        if setting.ampliref.check==1
            if length(setting.ampliref.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference'])
                uicontrol(h.pp_density_profile_setting_edit_ampliref);
                return;
            end
            if  setting.ampliref.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference that is > 0'])
                uicontrol(h.pp_density_profile_setting_edit_ampliref);
                return;
            end
        end
        
        if setting.levelquiver.check==1
            if isempty(setting.levelquiver.param)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify level/quiver option'])
                uicontrol(h.pp_density_profile_setting_edit_level);
                return;
            end
        end
      
         if setting.bathy.cb==1 && (isempty(setting.bathy.scale)||length(setting.bathy.scale)>1||setting.bathy.scale<0)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a scale for the bathymetry plot')
            uicontrol(h.pp_density_profile_edit_bathyScale);
            return;   
         end
        
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_density_plot_profile(tsnap,simuldata,setting,axesfig,h);
    end

    function callback_pp_cb_x(hObj,eventdata,h,href1,href2)
        Id=get(hObj,'value');
        if Id==1
            set(href1,'enable','on');
            set(href2,'enable','off');
        else
            set(href1,'enable','off');
            set(href2,'enable','on');
        end
    end

    function callback_pp_edit_x(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])

        if length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position at x axis'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_cb_y(hObj,eventdata,h,href1,href2)
        Id=get(hObj,'value');
        if Id==1
            set(href1,'enable','on','value',1);
            set(href2,'enable','off','value',0);
        else
            set(href1,'enable','off');
            set(href2,'enable','on');
            
        end
    end

    function callback_pp_edit_y(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])

        if length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position at y axis'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_density_signal_popup_var(hObj,eventdata,h)
        Id=get(hObj,'value');
        simuldata=get(h.pp_project_load_data,'userdata');
        
        if isempty(simuldata)
            set(hObj,'value',1);
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a simulation data'])
            return;
        else
            if ID_model_phiform==1
                if ~isfield(simuldata.output,'phi') &&Id>1
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, please specify in the option output variables and then re-run the simulation.'])
                    return;
                 end
                 if length(simuldata.output.phi(:,1,1))==1 &&Id>1
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, please specify in the option output variables and then re-run the simulation.'])
                    return;
                end
            else
                 if ~isfield(simuldata.output,'u') &&Id>1
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option output variables and then re-run the simulation.'])
                    return;
                 end
                 if length(simuldata.output.u(:,1,1))==1 &&Id>1
                    set(hObj,'value',1);
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option output variables and then re-run the simulation.'])
                    return;
                end
            end
        end
        set(h.monitorbox,'foregroundcolor','k','string','>>')
    end

    function callback_pp_plot_density_signal(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_density_signal_popup_var,'value');
        setting.view.check=get(h.pp_density_signal_setting_cb_view,'value');
        setting.view.param=get(h.pp_density_signal_setting_edit_view,'userdata');
        setting.clim.check=get(h.pp_density_signal_setting_cb_clim,'value');
        setting.clim.param=get(h.pp_density_signal_setting_edit_clim,'userdata');
        setting.spatlim.check=get(h.pp_density_signal_setting_cb_spatlim,'value');
        setting.spatlim.param=get(h.pp_density_signal_setting_edit_spatlim,'userdata');
        setting.tlim.check=get(h.pp_density_signal_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_density_signal_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_density_signal_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_density_signal_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_density_signal_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_density_signal_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_density_signal_setting_popup_savefig,'value'));
        colorfig=cellstr(get(h.pp_density_signal_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.pp_density_signal_setting_popup_colormap,'value'));
        axesfig=h.pp_density_signal_axes;
        spatsnap.x_check=get(h.pp_density_signal_x_cb,'value');
        spatsnap.x_snap=get(h.pp_density_signal_x_edit,'userdata');
        spatsnap.y_check=get(h.pp_density_signal_y_cb,'value');
        spatsnap.y_snap=get(h.pp_density_signal_y_edit,'userdata');
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if setting.var==2
            if ~isfield(simuldata.output,'u')
                set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                uicontrol(h.pp_density_signal_popup_var);
                return;
            end
            if length(simuldata.output.u(:,1,1))==1
                set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                uicontrol(h.pp_density_signal_popup_var);
                return;
            end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if spatsnap.x_check==0 && spatsnap.y_check==0
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x or y axis')
            uicontrol(h.pp_density_signal_x_edit);
            return;
        end
        
        if spatsnap.x_check==1
            if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x axis')
                uicontrol(h.pp_density_signal_x_edit);
                return;
            end
            
             if spatsnap.x_snap<X(1) || spatsnap.x_snap>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_density_signal_x_edit);
                return; 
            end
            
        end
        
        if spatsnap.y_check==1
            if isempty(spatsnap.y_snap)||length(spatsnap.y_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in y axis')
                uicontrol(h.pp_density_signal_y_edit);
                return;
            end
            
            if spatsnap.y_snap<Y(1) || spatsnap.y_snap>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_density_signal_y_edit);
                return; 
            end
        end
        if setting.view.check==1
            if isempty(setting.view.param)||length(setting.view.param)>2||length(setting.view.param)<1
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_signal_setting_edit_view);
                return;
            end
            if length(setting.view.param)==1 && (setting.view.param>3||setting.view.param<2)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_signal_setting_edit_view);
                return;
            end
        end
        if setting.clim.check==1
            if isempty(setting.clim.param)||length(setting.clim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify colorbar axes limit [zmin;zmax]'])
                uicontrol(h.pp_density_signal_setting_edit_clim);
                return;
            end
        end
        
        if setting.spatlim.check==1
            
            if spatsnap.y_check==1
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify spatial axes limit [xmin;xmax]'])
                uicontrol(h.pp_density_signal_setting_edit_spatlim);
                return;
                end
            
                if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1 || setting.spatlim.param(1)>setting.spatlim.param(2)
                    set(h.monitorbox,'foregroundcolor','r','string','>>Specify x axis limit [xmin;xmax]')
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<X(1) || setting.spatlim.param(1)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<X(1) || setting.spatlim.param(2)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
            else
                
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify spatial axes limit [ymin;ymax]'])
                uicontrol(h.pp_density_signal_setting_edit_spatlim);
                return;
                end
                
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2 || setting.spatlim.param(1)>setting.spatlim.param(2)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<Y(1) || setting.spatlim.param(1)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<Y(1) || setting.spatlim.param(2)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_density_signal_setting_edit_spatlim);
                    return;
                end
            end
            
        end
        
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_density_signal_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_density_signal_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_density_signal_setting_edit_tlim);
                return;
            end
            
        end
        if setting.coarse.check==1
            
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_density_signal_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_density_signal_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_density_plot_signal(spatsnap,simuldata,setting,axesfig,h.fig);
    end

%     function callback_cb_hs(hObj,evetdata,h)
%         Id=get(hObj,'value');
%         if Id==1
%             set(h.pp_density_statistic_hs_cb,'enable','on','value',1)
%             set(h.pp_density_statistic_mtc_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','off','value',0)
%         else
%             set(h.pp_density_statistic_hs_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','on','value',0)
%         end
%     end



%     function callback_cb_mtc(hObj,eventdata,h)
%         Id=get(hObj,'value');
%         if Id==1
%             set(h.pp_density_statistic_hs_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','on','value',1)
%             set(h.pp_density_statistic_mtt_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','off','value',0)
%         else
%             set(h.pp_density_statistic_hs_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','on','value',0)
%         end
%     end


%     function callback_cb_mtt(hObj,eventdata,h)
%         Id=get(hObj,'value');
%         if Id==1
%             set(h.pp_density_statistic_hs_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','on','value',1)
%             set(h.pp_density_statistic_mwl_cb,'enable','off','value',0)
%         else
%             set(h.pp_density_statistic_hs_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','on','value',0)
%         end
%     end
% 
%     function callback_cb_mwl(hObj,eventdata,h)
%         Id=get(hObj,'value');
%         if Id==1
%             set(h.pp_density_statistic_hs_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','off','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','on','value',1)
%         else
%             set(h.pp_density_statistic_hs_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtc_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mtt_cb,'enable','on','value',0)
%             set(h.pp_density_statistic_mwl_cb,'enable','on','value',0)
%         end
%     end



    function callback_pp_plot_density_statistic(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>> calculating')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.view.check=get(h.pp_density_statistic_setting_cb_view,'value');
        setting.view.param=get(h.pp_density_statistic_setting_edit_view,'userdata');
        setting.clim.check=get(h.pp_density_statistic_setting_cb_clim,'value');
        setting.clim.param=get(h.pp_density_statistic_setting_edit_clim,'userdata');
        setting.xlim.check=get(h.pp_density_statistic_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_density_statistic_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_density_statistic_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_density_statistic_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_density_statistic_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_density_statistic_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_density_statistic_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_density_statistic_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_density_statistic_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_density_statistic_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_density_statistic_setting_popup_savefig,'value'));
        colorfig=cellstr(get(h.pp_density_statistic_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.pp_density_statistic_setting_popup_colormap,'value'));
        axesfig=h.pp_density_statistic_axes;
        setting.savedata=get(h.pp_density_statistic_savedata,'value');
        
        cbvar=get(h.pp_density_statistic_popup_var,'value');
%         if get(h.pp_density_statistic_hs_cb,'value')==1,cbvar=1;
%         elseif get(h.pp_density_statistic_mtc_cb,'value')==1,cbvar=2;
%         elseif get(h.pp_density_statistic_mtt_cb,'value')==1,cbvar=3;
%         elseif get(h.pp_density_statistic_mwl_cb,'value')==1,cbvar=4;
%         end
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        
        if cbvar==0
            set(h.monitorbox,'foregroundcolor','r','string','>> Choose a statistic variable')
            uicontrol(h.pp_density_statistic_hs_cb);
            return;
        end
        
        if setting.view.check==1
            if isempty(setting.view.param)||length(setting.view.param)>2||length(setting.view.param)<1
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_statistic_setting_edit_view);
                return;
            end
            if length(setting.view.param)==1 && (setting.view.param>3||setting.view.param<2)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_density_statistic_setting_edit_view);
                return;
            end
        end
        if setting.clim.check==1
            if isempty(setting.clim.param)||length(setting.clim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a colorbar axes limit [zmin;zmax]'])
                uicontrol(h.pp_density_statistic_setting_edit_clim);
                return;
            end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;

         if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_density_statistic_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_density_statistic_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_density_statistic_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_density_statistic_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_density_statistic_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_density_statistic_setting_edit_ylim);
                return; 
            end
        end
        
  
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_density_statistic_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_density_statistic_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_density_statistic_setting_edit_tlim);
                return;
            end
            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_density_statistic_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_density_statistic_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_density_plot_statistic(cbvar,simuldata,setting,axesfig,h);
        
        
    end

    function callback_line_profile_at_time(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])

        if isempty(param)
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify a time or multiple time: t1;t2;t3')
            uicontrol(hObj);
        end
    end

    function callback_pp_cb_signal(hObj,~,hscb,hspopup,hssmooth1,hssmooth2)
        Id=get(hObj,'value');
        if Id==1
            set(hscb,'enable','off','value',0);
            set(hspopup,'enable','off','value',1,'visible','on');
            set(hssmooth1,'enable','off','value',0);set(hssmooth2,'enable','off');
        else
            set(hscb,'enable','on','value',1);
            set(hspopup,'enable','on','value',1,'visible','on');
            set(hObj,'enable','off','value',0)
            set(hssmooth1,'enable','on','value',0);set(hssmooth2,'enable','off');
        end
    end

    function callback_pp_cb_spectrum(hObj,~,hspopup,hsignal,hssmooth1,hssmooth2)
        Id=get(hObj,'value');
        if Id==1
            set(hsignal,'enable','off','value',0);
            set(hspopup,'enable','on','value',1,'visible','on');
            set(hssmooth1,'enable','on','value',0);set(hssmooth2,'enable','off');
        else
            set(hsignal,'enable','on','value',1);
            set(hspopup,'enable','off','value',1,'visible','on');
            set(hObj,'enable','off','value',0)
            set(hssmooth1,'enable','off','value',0);set(hssmooth2,'enable','off');
        end
    end

    function callback_pp_MTA_check(hObj,eventdata,h)
        Id1=get(h.pp_line_profile_MTC_cb,'value');
        Id2=get(h.pp_line_profile_MTT_cb,'value');
        if Id1+Id2>0
         set(h.pp_line_profile_setting_cb_MTAtlim,'enable','on')
        else
         set(h.pp_line_profile_setting_cb_MTAtlim,'enable','off','value',0) 
         set(h.pp_line_profile_setting_edit_MTAtlim,'enable','off');
        end
    end

    function callback_pp_MTA_tlim_check(hObj,eventdata,h)
       Id=get(hObj,'value');
       if Id==1
          set(h.pp_line_profile_setting_edit_MTAtlim,'enable','on') 
       else
          set(h.pp_line_profile_setting_edit_MTAtlim,'enable','off') 
       end
    end

    function callback_pp_plot_line_profile(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_profile_popup_var,'value');
        
        setting.zlim.check=get(h.pp_line_profile_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_profile_setting_edit_zlim,'userdata');
        setting.spatlim.check=get(h.pp_line_profile_setting_cb_spatlim,'value');
        setting.spatlim.param=get(h.pp_line_profile_setting_edit_spatlim,'userdata');
        setting.MTAtlim.check=get(h.pp_line_profile_setting_cb_MTAtlim,'value');
        setting.MTAtlim.param=get(h.pp_line_profile_setting_edit_MTAtlim,'userdata');
        setting.coarse.check=get(h.pp_line_profile_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_profile_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_line_profile_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_profile_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_profile_setting_popup_savefig,'value'));
        axesfig=h.pp_line_profile_axes;
        tsnap=get(h.pp_line_profile_time_edit,'userdata');
        spatsnap.x_check=get(h.pp_line_profile_x_cb,'value');
        spatsnap.x_snap=get(h.pp_line_profile_x_edit,'userdata');
        spatsnap.y_check=get(h.pp_line_profile_y_cb,'value');
        spatsnap.y_snap=get(h.pp_line_profile_y_edit,'userdata');
        setting.MTCcheck=get(h.pp_line_profile_MTC_cb,'value');
        setting.MTTcheck=get(h.pp_line_profile_MTT_cb,'value');
        setting.bathy.cb=get(h.pp_line_profile_cb_bathy,'value');
        setting.bathy.scale=get(h.pp_line_profile_edit_bathyScale,'userdata');

        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if setting.var==2
           if ~isfield(simuldata.output,'u')
                   set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                   uicontrol(h.pp_line_profile_popup_var);
                   return;
           end
           if length(simuldata.output.u(:,1,1))==1
               set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
               uicontrol(h.pp_line_profile_popup_var);
               return;
           end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        if isempty(tsnap)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a single time or multiple time: t1;t2t3')
            uicontrol(h.pp_line_profile_time_edit);
            return;
        end
        if  any(tsnap<T(1)) || any(tsnap>T(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>> The specified time must be in the interval ['...
                ,num2str(T(1)),';',num2str(T(end)),'] (s)'])
            uicontrol(h.pp_line_profile_time_edit);
            return;
        end
        
        if spatsnap.x_check==0 && spatsnap.y_check==0
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x or y axis')
            uicontrol(h.pp_line_profile_x_edit);
            return;
        end
        
        if spatsnap.x_check==1
            if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x axis')
                uicontrol(h.pp_line_profile_x_edit);
                return;
            end
            
             if spatsnap.x_snap<X(1) || spatsnap.x_snap>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_profile_x_edit);
                return; 
            end
            
        end
        if spatsnap.y_check==1
            if isempty(spatsnap.y_snap)||length(spatsnap.y_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in y axis')
                uicontrol(h.pp_line_profile_y_edit);
                return;
            end
            
            if spatsnap.y_snap<Y(1) || spatsnap.y_snap>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_profile_y_edit);
                return; 
            end
        end
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_profile_setting_edit_clim);
                return;
            end
        end
        
        if setting.spatlim.check==1
            if spatsnap.y_check==1
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify spatial axes limit [xmin;xmax]'])
                uicontrol(h.pp_line_profile_setting_edit_spatlim);
                return;
                end
                
                %if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1 || setting.spatlim.param(1)>setting.spatlim.param(2)
                %Nida
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2 || setting.spatlim.param(1)>setting.spatlim.param(2)                       
                    set(h.monitorbox,'foregroundcolor','r','string','>>Specify x axis limit [xmin;xmax]')
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<X(1) || setting.spatlim.param(1)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<X(1) || setting.spatlim.param(2)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
            else
                
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify spatial axes limit [ymin;ymax]'])
                uicontrol(h.pp_line_profile_setting_edit_spatlim);
                return;
                end
                
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2 || setting.spatlim.param(1)>setting.spatlim.param(2)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<Y(1) || setting.spatlim.param(1)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<Y(1) || setting.spatlim.param(2)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_line_profile_setting_edit_spatlim);
                    return;
                end
            end
        end
        
        if setting.MTAtlim.check==1
             if isempty(setting.MTAtlim.param)||length(setting.MTAtlim.param)~=2 || setting.MTAtlim.param(1)>setting.MTAtlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify MTA time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_profile_setting_edit_MTAtlim);
                return;
            end
            
            if setting.MTAtlim.param(1)<T(1) || setting.MTAtlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify MTA time limit in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_profile_setting_edit_MTAtlim);
                return;
            end
            
            if setting.MTAtlim.param(2)<T(1) || setting.MTAtlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify MTA time limit in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_profile_setting_edit_MTAtlim);
                return;
            end
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_profile_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_profile_setting_edit_coarse);
                return;
            end
        end
        
        if setting.bathy.cb==1 && (isempty(setting.bathy.scale)||length(setting.bathy.scale)>1||setting.bathy.scale<0)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a scale for the bathymetry plot')
            uicontrol(h.pp_line_profile_edit_bathyScale);
            return;   
         end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_profile(tsnap,spatsnap,simuldata,setting,axesfig,h);
    end

    function callback_pp_line_buoy_edit_x(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])

        if isempty(param)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position at x axis'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_pp_line_buoy_edit_y(hObj,eventdata,h)
        param=str2num(get(hObj,'string'));
        set(hObj,'userdata',param)
        set(h.monitorbox,'foregroundcolor','k','string',['>>'])

        if isempty(param)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position at y axis'])
            uicontrol(hObj);
            return;
        end
    end

    function callback_cb_combine_buoys(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.pp_line_buoy_combine_tinterv_text,'enable','on');
            set(h.pp_line_buoy_combine_tinterv_edit,'enable','on');
        else
            set(h.pp_line_buoy_combine_tinterv_text,'enable','off');
            set(h.pp_line_buoy_combine_tinterv_edit,'enable','off');
        end
    end


    function callback_pp_plot_line_buoy(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_buoy_popup_var,'value');
        setting.cb.signal=get(h.pp_line_buoy_cb_signal,'value');
        setting.cb.spectrum=get(h.pp_line_buoy_cb_spectrum,'value');
        setting.cb.spectrum_var=get(h.pp_line_buoy_popup_spectrum_var,'value');
        setting.cb.combine_buoys=get(h.pp_line_buoy_cb_combine_buoys,'value');
        setting.cb.combine_tinterv=get(h.pp_line_buoy_combine_tinterv_edit,'userdata');
        
        setting.zlim.check=get(h.pp_line_buoy_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_buoy_setting_edit_zlim,'userdata');
        setting.horznlim.check=get(h.pp_line_buoy_setting_cb_horznlim,'value');
        setting.horznlim.param=get(h.pp_line_buoy_setting_edit_horznlim,'userdata');
        setting.tlim.check=get(h.pp_line_buoy_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_buoy_setting_edit_tlim,'userdata');
        
        setting.coarse.check=get(h.pp_line_buoy_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_buoy_setting_edit_coarse,'userdata');
        setting.spsmooth.check=get(h.pp_line_buoy_setting_cb_spsmooth,'value');
        setting.spsmooth.param=get(h.pp_line_buoy_setting_edit_spsmooth,'userdata');
        
        setting.savefig.check=get(h.pp_line_buoy_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_buoy_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_buoy_setting_popup_savefig,'value'));
        axesfig=h.pp_line_buoy_axes;
        spatsnap.x_snap=get(h.pp_line_buoy_x_edit,'userdata');
        spatsnap.y_snap=get(h.pp_line_buoy_y_edit,'userdata');
        
        setting.savedat=get(h.pp_line_buoy_savedata,'value');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
       
        if setting.var==2
           if ~isfield(simuldata.output,'u')
                   set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                   uicontrol(h.pp_line_buoy_popup_var);
                   return;
           end
           if length(simuldata.output.u(:,1,1))==1
               set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
               uicontrol(h.pp_line_buoy_popup_var);
               return;
           end
        end
        
        
        if isempty(spatsnap.x_snap)
            set(h.monitorbox,'foregroundcolor','r','string','Specify a position or multiple position in x axis, x:[x1;x2;x2]')
            uicontrol(h.pp_line_buoy_x_edit);
            return;
        end
        
        if any(spatsnap.x_snap<X(1)) || any(spatsnap.x_snap>X(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
            uicontrol(h.pp_line_profile_x_edit);
            return;
        end
            
        
        if isempty(spatsnap.y_snap)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position or multiple position in y axis, y:[y1;y2;y3]')
            uicontrol(h.pp_line_buoy_y_edit);
            return;
        end
        
        
        if any(spatsnap.y_snap<Y(1)) || any(spatsnap.y_snap>Y(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
            uicontrol(h.pp_line_profile_y_edit);
            return;
        end
        
        if length(spatsnap.x_snap)~=length(spatsnap.y_snap)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Number of position in x and y axes must be the same']);
            uicontrol(h.pp_line_profile_y_edit);
            return;
        end
        
        if setting.cb.combine_buoys==1
            if length(spatsnap.x_snap)<=1 || length(spatsnap.y_snap)<=1
                set(h.pp_line_buoy_cb_signal,'value',0)
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify another buoy position for combining signals')
                uicontrol(h.pp_line_buoy_x_edit);
                return;
            end
            if length(setting.cb.combine_tinterv)~=2
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a time interval for combining signals')
                uicontrol(h.pp_line_buoy_combine_tinterv_edit);
                return;
            end
            
            if setting.cb.combine_tinterv(2)<setting.cb.combine_tinterv(1)
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a time interval with format (t_min;t_max)')
                uicontrol(h.pp_line_buoy_combine_tinterv_edit);
                return;
            end
            
            if any(setting.cb.combine_tinterv<T(1)) || any(setting.cb.combine_tinterv>T(end))
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time interval in [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_buoy_setting_edit_tlim);
                return;
            end
            
            if setting.cb.spectrum_var>2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Combine buoys is not avialable for directional spect']);
                uicontrol(h.pp_line_buoy_popup_var);
                return;
            end
        end
        
        
        if setting.cb.spectrum>0 && setting.var==2
            set(h.monitorbox,'foregroundcolor','r','string','>> Choose only an elevation input variable for plotting spectrum');
            uicontrol(h.pp_line_buoy_popup_var);
            return;
        end
        
      
        
        if setting.cb.spectrum_var>2 && length(spatsnap.x_snap)<3
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify at least 3 buoy positions');
            uicontrol(h.pp_line_profile_x_edit);
            return;
        end
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_buoy_setting_edit_zlim);
                return;
            end
        end
        
        if setting.horznlim.check==1
            if isempty(setting.horznlim.param)||length(setting.horznlim.param)~=2 || setting.horznlim.param(1)>=setting.horznlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a horizontal axes limit, [min;max]'])
                uicontrol(h.pp_line_buoy_setting_edit_horznlim);
                return;
            end
        end
        
        
       if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_buoy_setting_edit_tlim);
                return;
            end
           
            if setting.cb.combine_buoys==0
                if any(setting.tlim.param<T(1)-0.0001) || any(setting.tlim.param>T(end)+0.0001)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                    uicontrol(h.pp_line_buoy_setting_edit_tlim);
                    return;
                end
            end
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
        end
        
        
        if setting.spsmooth.check==1
            if length(setting.spsmooth.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a smooth factor for spectrum'])
                uicontrol(h.pp_line_buoy_setting_cb_spsmooth);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_buoy(spatsnap,simuldata,setting,axesfig,h.fig);
        
    end


    function callback_pp_statistic_cb_along_axes(hObj,eventdata,h)
        Id=get(hObj,'value');
        IdVar=get(h.pp_line_statistic_popup_var,'value');
        if IdVar>10
            Id=0;
            set(hObj,'value',0);
        end
        if Id==1
            set(h.pp_line_statistic_x_cb,'enable','on');
            set(h.pp_line_statistic_x_edit,'enable','off');
            set(h.pp_line_statistic_y_cb,'enable','on');
            set(h.pp_line_statistic_y_edit,'enable','off');
            set(h.pp_line_statistic_buoys_cb,'value',0);
            set(h.pp_line_statistic_buoy_x_text,'enable','off')
            set(h.pp_line_statistic_buoy_x_edit,'enable','off')
            set(h.pp_line_statistic_buoy_y_text,'enable','off')
            set(h.pp_line_statistic_buoy_y_edit,'enable','off')
        else
            set(h.pp_line_statistic_x_cb,'enable','off');
            set(h.pp_line_statistic_x_edit,'enable','off');
            set(h.pp_line_statistic_y_cb,'enable','off');
            set(h.pp_line_statistic_y_edit,'enable','off');
            set(h.pp_line_statistic_buoys_cb,'value',1);
            set(h.pp_line_statistic_buoy_x_text,'enable','on')
            set(h.pp_line_statistic_buoy_x_edit,'enable','on')
            set(h.pp_line_statistic_buoy_y_text,'enable','on')
            set(h.pp_line_statistic_buoy_y_edit,'enable','on')
        end
    end

    function callback_pp_statistic_cb_at_buoys(hObj,eventdata,h)
        Id=get(hObj,'value');
        IdVar=get(h.pp_line_statistic_popup_var,'value');
        if IdVar==11
            Id=1;
            set(hObj,'value',1);
        end
        if Id==1
            set(h.pp_line_statistic_x_cb,'enable','off');
            set(h.pp_line_statistic_x_edit,'enable','off');
            set(h.pp_line_statistic_y_cb,'enable','off');
            set(h.pp_line_statistic_y_edit,'enable','off');
            set(h.pp_line_statistic_along_axes_cb,'value',0);
            set(h.pp_line_statistic_buoy_x_text,'enable','on')
            set(h.pp_line_statistic_buoy_x_edit,'enable','on')
            set(h.pp_line_statistic_buoy_y_text,'enable','on')
            set(h.pp_line_statistic_buoy_y_edit,'enable','on')
        else
            if IdVar<11
                set(h.pp_line_statistic_x_cb,'enable','on');
                set(h.pp_line_statistic_x_edit,'enable','off');
                set(h.pp_line_statistic_y_cb,'enable','on');
                set(h.pp_line_statistic_y_edit,'enable','off');
                set(h.pp_line_statistic_along_axes_cb,'value',1);
                
                set(h.pp_line_statistic_buoy_x_text,'enable','off')
                set(h.pp_line_statistic_buoy_x_edit,'enable','off')
                set(h.pp_line_statistic_buoy_y_text,'enable','off')
                set(h.pp_line_statistic_buoy_y_edit,'enable','off')
            else
                set(h.pp_line_statistic_x_cb,'enable','off');
                set(h.pp_line_statistic_x_edit,'enable','off');
                set(h.pp_line_statistic_y_cb,'enable','off');
                set(h.pp_line_statistic_y_edit,'enable','off');
                set(hObj,'value',1)
                set(h.pp_line_statistic_along_axes_cb,'value',0);
                set(h.pp_line_statistic_buoy_x_text,'enable','on')
                set(h.pp_line_statistic_buoy_x_edit,'enable','on')
                set(h.pp_line_statistic_buoy_y_text,'enable','on')
                set(h.pp_line_statistic_buoy_y_edit,'enable','on')
            end
            
        end
    end

    function callback_pp_plot_line_statistic_var(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id>10
            set(h.pp_line_statistic_x_cb,'enable','off');
            set(h.pp_line_statistic_x_edit,'enable','off');
            set(h.pp_line_statistic_y_cb,'enable','off');
            set(h.pp_line_statistic_y_edit,'enable','off');
            set(h.pp_line_statistic_along_axes_cb,'value',0);
            set(h.pp_line_statistic_buoys_cb,'value',1)
            set(h.pp_line_statistic_buoy_x_text,'enable','on')
            set(h.pp_line_statistic_buoy_x_edit,'enable','on')
            set(h.pp_line_statistic_buoy_y_text,'enable','on')
            set(h.pp_line_statistic_buoy_y_edit,'enable','on')
            if Id>11
                set(h.pp_line_statistic_buoys_cb,'string','Area interval:','value',1)
                set(h.pp_line_statistic_buoy_x_text,'enable','on')
                set(h.pp_line_statistic_buoy_x_edit,'enable','on')
                set(h.pp_line_statistic_buoy_y_text,'enable','on')
                set(h.pp_line_statistic_buoy_y_edit,'enable','on')
            else
             set(h.pp_line_statistic_buoys_cb,'string','at buoy(s):')   
            end
            set(h.pp_line_statistic_setting_cb_Hs_ref,'enable','on','value',0);
            set(h.pp_line_statistic_setting_edit_Hs_ref,'enable','off');
        else
            set(h.pp_line_statistic_buoys_cb,'string','at buoy(s):')
            set(h.pp_line_statistic_setting_cb_Hs_ref,'enable','off','value',0);
            set(h.pp_line_statistic_setting_edit_Hs_ref,'enable','off');
        end
    end

    function callback_pp_plot_line_statistic(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_statistic_popup_var,'value');
        
        setting.zlim.check=get(h.pp_line_statistic_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_statistic_setting_edit_zlim,'userdata');
        setting.spatlim.check=get(h.pp_line_statistic_setting_cb_spatlim,'value');
        setting.spatlim.param=get(h.pp_line_statistic_setting_edit_spatlim,'userdata');
        setting.tlim.check=get(h.pp_line_statistic_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_statistic_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_line_statistic_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_statistic_setting_edit_coarse,'userdata');
        setting.Hsref.check=get(h.pp_line_statistic_setting_cb_Hs_ref,'value');
        setting.Hsref.param=get(h.pp_line_statistic_setting_edit_Hs_ref,'userdata');
        setting.savefig.check=get(h.pp_line_statistic_setting_cb_savefig,'value');
        setting.savedat=get(h.pp_line_statistic_savedata,'value');
        
        format=cellstr(get(h.pp_line_statistic_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_statistic_setting_popup_savefig,'value'));
        axesfig=h.pp_line_statistic_axes;
        spatsnap.along_axes_check=get(h.pp_line_statistic_along_axes_cb,'value');
        spatsnap.x_check=get(h.pp_line_statistic_x_cb,'value');
        spatsnap.x_snap=get(h.pp_line_statistic_x_edit,'userdata');
        spatsnap.y_check=get(h.pp_line_statistic_y_cb,'value');
        spatsnap.y_snap=get(h.pp_line_statistic_y_edit,'userdata');
       
        if setting.var<12
            spatsnap.buoy_check=get(h.pp_line_statistic_buoys_cb,'value');
            spatsnap.buoy_x_snap=get(h.pp_line_statistic_buoy_x_edit,'userdata');
            spatsnap.buoy_y_snap=get(h.pp_line_statistic_buoy_y_edit,'userdata');
            
        else
            spatsnap.area_check=get(h.pp_line_statistic_buoys_cb,'value');
            spatsnap.area_x_interv=get(h.pp_line_statistic_buoy_x_edit,'userdata');
            spatsnap.area_y_interv=get(h.pp_line_statistic_buoy_y_edit,'userdata');
        end
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if spatsnap.along_axes_check==1
            if spatsnap.x_check==0 && spatsnap.y_check==0
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x or y axis')
                uicontrol(h.pp_line_statistic_x_edit);
                return;
            end
            
            if spatsnap.x_check==1
                if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1
                    set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x axis')
                    uicontrol(h.pp_line_statistic_x_edit);
                    return;
                end
                
                if spatsnap.x_snap<X(1) || spatsnap.x_snap>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_line_statistic_x_edit);
                    return;
                end
                
            end
            
            if spatsnap.y_check==1
                if isempty(spatsnap.y_snap)||length(spatsnap.y_snap)~=1
                    set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in y axis')
                    uicontrol(h.pp_line_statistic_y_edit);
                    return;
                end
                
                if spatsnap.y_snap<Y(1) || spatsnap.y_snap>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_line_statistic_y_edit);
                    return;
                end
            end
            
        else
            if setting.var<=11 
                if isempty(spatsnap.buoy_x_snap)
                    set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x axis')
                    uicontrol(h.pp_line_statistic_buoy_x_edit);
                    return;
                end
                
                if any(spatsnap.buoy_x_snap<X(1)) || any(spatsnap.buoy_x_snap>X(end))
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_line_statistic_buoy_x_edit);
                    return;
                end
                
                if isempty(spatsnap.buoy_y_snap)
                    set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in y axis')
                    uicontrol(h.pp_line_statistic_buoy_y_edit);
                    return;
                end
                
                if any(spatsnap.buoy_y_snap<Y(1)) || any(spatsnap.buoy_y_snap>Y(end))
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_line_statistic_buoy_y_edit);
                    return;
                end
                
                if length(spatsnap.buoy_x_snap)~=length(spatsnap.buoy_y_snap)
                    set(h.monitorbox,'foregroundcolor','r','string','>> Specify positions in (x,y)')
                    uicontrol(h.pp_line_statistic_buoy_y_edit);
                    return;
                end
            else
                
                if spatsnap.area_check==1
                    if isempty(spatsnap.area_x_interv)
                        set(h.monitorbox,'foregroundcolor','r','string','>> Specify x interval(xstart;xend)')
                        uicontrol(h.pp_line_statistic_buoy_x_edit);
                        return;
                    end
                    if isempty(spatsnap.area_y_interv)
                        set(h.monitorbox,'foregroundcolor','r','string','>> Specify y interval')
                        uicontrol(h.pp_line_statistic_buoy_y_edit);
                        return;
                    end
                    if length(spatsnap.area_x_interv)~=2
                        set(h.monitorbox,'foregroundcolor','r','string','>> Specify x interval (xstart;xend)')
                        uicontrol(h.pp_line_statistic_buoy_x_edit);
                        return;
                    end
                    
                    if length(spatsnap.area_y_interv)~=2
                        set(h.monitorbox,'foregroundcolor','r','string','>> Specify y interval (ystart;yend)')
                        uicontrol(h.pp_line_statistic_buoy_y_edit);
                        return;
                    end
                    
                     if any(spatsnap.area_x_interv<X(1)) || any(spatsnap.area_x_interv>X(end))
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x interval in [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                        uicontrol(h.pp_line_statistic_buoy_x_edit);
                        return;
                     end
                    
                     if any(spatsnap.area_y_interv<Y(1)) || any(spatsnap.area_y_interv>Y(end))
                        set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y interval in [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                        uicontrol(h.pp_line_statistic_buoy_y_edit);
                        return;
                    end
                    
                end
                
            end
            
        end
        
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_statistic_setting_edit_zlim);
                return;
            end
        end
        
        if setting.spatlim.check==1
            if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a horizontal axes limit [xmin;xmax] or [ymin;ymax]'])
                uicontrol(h.pp_line_statistic_setting_edit_xlim);
                return;
            end
        end
        
        if setting.tlim.check==1
             
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_statistic_profile_setting_edit_tlim);
                return;
            end
            
            if any(setting.tlim.param<T(1)-0.0001) || any(setting.tlim.param>T(end)+0.0001)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_statistic_profile_setting_edit_tlim);
                return;
            end
            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
            
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        if spatsnap.along_axes_check==1
            set(axesfig,'visible','on');
            set(findall(h.pp_line_statistic_toolbar_panel),'visible','on');
            set(h.pp_line_statistic_show_value_edit,'visible','off')
            funP_line_plot_statistic(spatsnap,simuldata,setting,axesfig,h);
        else
            cla(axesfig);
            if setting.var<=10
                set(axesfig,'visible','off');
                set(h.pp_line_statistic_show_value_edit,'visible','on')
                set(findall(h.pp_line_statistic_toolbar_panel),'visible','off');
            else
                set(axesfig,'visible','on');
                set(findall(h.pp_line_statistic_toolbar_panel),'visible','on');
                set(h.pp_line_statistic_show_value_edit,'visible','off')
            end
            if setting.var<=11 
                funP_line_plot_statistic_buoy(spatsnap,simuldata,setting,axesfig,h);
            elseif setting.var==12 
                %funP_line_plot_statistic_area_exceedance(simuldata,spatsnap,setting,axesfig,h);
                funP_line_plot_statistic_area_exceedance_of_elevation(simuldata,spatsnap,setting,axesfig,h)
            elseif setting.var==13 %% not used in GUI
                funP_line_plot_statistic_time_exceedance(simuldata,spatsnap,setting,axesfig,h)
            end
        end
        
    end

    function callback_pp_plot_line_ham_mom(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_ham_mom_popup_var,'value');
        
        setting.zlim.check=get(h.pp_line_ham_mom_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_ham_mom_setting_edit_zlim,'userdata');
        setting.xlim.check=get(h.pp_line_ham_mom_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_line_ham_mom_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_line_ham_mom_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_line_ham_mom_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_line_ham_mom_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_ham_mom_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_line_ham_mom_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_ham_mom_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_line_ham_mom_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_ham_mom_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_ham_mom_setting_popup_savefig,'value'));
        axesfig=h.pp_line_ham_mom_axes;
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if ID_model_phiform==1
           if ~isfield(simuldata.output,'phi')
                set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option for output variables and then re-run the simulation.'])
                return;
            end
            
            if length(squeeze(squeeze(simuldata.output.phi(:,1,1))))==1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Wave potential data is not found, please specify in the option for output variables and then re-run the simulation.'])
                return;
            end 
        else
            if ~isfield(simuldata.output,'u')
                set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option for output variables and then re-run the simulation.'])
                return;
            end
            
            if length(squeeze(squeeze(simuldata.output.u(:,1,1))))==1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Wave velocity data is not found, please specify in the option for output variables and then re-run the simulation.'])
                return;
            end
        end
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_ham_mom_setting_edit_zlim);
                return;
            end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;

         if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_line_ham_mom_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_line_ham_mom_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_ylim);
                return; 
            end
        end
        
  
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_ham_mom_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_ham_mom_setting_edit_tlim);
                return;
            end
            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_ham_mom_setting_edit_coarse);
                return;
            end
            
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_ham_mom_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_ham_mom(simuldata,setting,axesfig,h);
    end

    function callback_line_breaking_popup_var(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.pp_line_breaking_setting_cb_zlim,'enable','off','value',0);
            set(h.pp_line_breaking_setting_edit_zlim,'enable','off')
            set(h.pp_line_breaking_setting_cb_xlim,'enable','on','value',0);
            set(h.pp_line_breaking_setting_edit_xlim,'enable','off')
            set(h.pp_line_breaking_setting_cb_ylim,'enable','on','value',0);
            set(h.pp_line_breaking_setting_edit_ylim,'enable','off')
            set(h.pp_line_breaking_setting_cb_coarse,'value',0,'enable','off');
            set(h.pp_line_breaking_setting_edit_coarse,'enable','off')
        else
            set(h.pp_line_breaking_setting_cb_zlim,'enable','on','value',0);
            set(h.pp_line_breaking_setting_edit_zlim,'enable','off')
            set(h.pp_line_breaking_setting_cb_xlim,'enable','off','value',0);
            set(h.pp_line_breaking_setting_edit_xlim,'enable','off')
            set(h.pp_line_breaking_setting_cb_ylim,'enable','off','value',0);
            set(h.pp_line_breaking_setting_edit_ylim,'enable','off')
            set(h.pp_line_breaking_setting_cb_coarse,'value',0,'enable','on');
            set(h.pp_line_breaking_setting_edit_coarse,'enable','off')
        end
    end

    function callback_pp_plot_line_breaking(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_breaking_popup_var,'value');
        
        setting.zlim.check=get(h.pp_line_breaking_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_breaking_setting_edit_zlim,'userdata');
        setting.xlim.check=get(h.pp_line_breaking_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_line_breaking_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_line_breaking_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_line_breaking_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_line_breaking_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_breaking_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_line_breaking_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_breaking_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_line_breaking_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_breaking_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_breaking_setting_popup_savefig,'value'));
        axesfig=h.pp_line_breaking_axes;
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
    
        if ~isfield(simuldata.output,'break_nodes')
            set(h.monitorbox,'foregroundcolor','r','string','>> No breaking data available!')
            return;
        end
        
        if isempty(simuldata.output.break_nodes)
            set(h.monitorbox,'foregroundcolor','r','string','>> No breaking data available!')
            return;
        end
        
        if setting.var==1
            if ~isfield(simuldata.output,'break_crest')
                set(h.monitorbox,'foregroundcolor','r','string',['>>No breaking data available!'])
                return;
            end
            
            if ~isfield(simuldata.output.break_crest,'X')
                set(h.monitorbox,'foregroundcolor','r','string',['>>Old format data is detected, please re-run the simulation'])
                return;
            end
            
            if isempty(simuldata.output.break_crest)
                set(h.monitorbox,'foregroundcolor','r','string',['>>No breaking data available!'])
                return;
            end
            
        end
        
        if setting.var==2
            if ~isfield(simuldata.output,'break_speed')
                set(h.monitorbox,'foregroundcolor','r','string',['>>No breaking data available!'])
                return;
            end
            
            if isempty(simuldata.output.break_speed)
                set(h.monitorbox,'foregroundcolor','r','string',['>>No breaking data available!'])
                return;
            end
        end
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_breaking_setting_edit_zlim);
                return;
            end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;

         if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_line_breaking_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_breaking_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_breaking_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_line_breaking_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_breaking_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_breaking_setting_edit_ylim);
                return; 
            end
        end
        
  
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_breaking_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_breaking_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_breaking_setting_edit_tlim);
                return;
            end
            
        end
        
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_breaking_setting_edit_coarse);
                return;
            end
            
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_breaking_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_breaking(simuldata,setting,axesfig,h);
    end

    function callback_pp_plot_line_MTAA(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_line_MTAA_popup_var,'value');
        setting.plotopt=get(h.pp_line_MTAA_popup_plotopt,'value');
        setting.zlim.check=get(h.pp_line_MTAA_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_MTAA_setting_edit_zlim,'userdata');
        setting.xlim.check=get(h.pp_line_MTAA_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_line_MTAA_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_line_MTAA_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_line_MTAA_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_line_MTAA_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_MTAA_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_line_MTAA_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_MTAA_setting_edit_coarse,'userdata');
        setting.savefig.check=get(h.pp_line_MTAA_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_MTAA_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_MTAA_setting_popup_savefig,'value'));
        axesfig=h.pp_line_MTAA_axes;
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if setting.var==2
            if ID_model_phiform==1
                 if ~isfield(simuldata.output,'phi')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_line_MTAA_popup_var);
                    return;
                end
                if length(simuldata.output.phi(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_line_MTAA_popup_var);
                    return;
                end
            else
                if ~isfield(simuldata.output,'u')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_line_MTAA_popup_var);
                    return;
                end
                if length(simuldata.output.u(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_line_MTAA_popup_var);
                    return;
                end
            end
        end
        
        
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_MTAA_setting_edit_zlim);
                return;
            end
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;

         if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_line_MTAA_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_MTAA_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_MTAA_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_line_MTAA_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_MTAA_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_MTAA_setting_edit_ylim);
                return; 
            end
        end
        
  
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_MTAA_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_MTAA_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_MTAA_setting_edit_tlim);
                return;
            end
            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_MTAA_setting_edit_coarse);
                return;
            end
            
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_MTAA_setting_edit_coarse);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_MTAA(simuldata,setting,axesfig,h);
    end

    function callback_pp_plot_line_Extreme(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.plotopt=get(h.pp_line_Extreme_popup_plotopt,'value');
        setting.zlim.check=get(h.pp_line_Extreme_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_line_Extreme_setting_edit_zlim,'userdata');
        setting.xlim.check=get(h.pp_line_Extreme_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_line_Extreme_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_line_Extreme_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_line_Extreme_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_line_Extreme_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_line_Extreme_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_line_Extreme_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_line_Extreme_setting_edit_coarse,'userdata');
        setting.ampliref.check=get(h.pp_line_Extreme_setting_cb_ampliref,'value');
        setting.ampliref.param=get(h.pp_line_Extreme_setting_edit_ampliref,'userdata');
        setting.savefig.check=get(h.pp_line_Extreme_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_line_Extreme_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_line_Extreme_setting_popup_savefig,'value'));
        axesfig=h.pp_line_Extreme_axes;
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_line_Extreme_setting_edit_zlim);
                return;
            end
        end
        
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;

         if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_line_Extreme_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_Extreme_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_line_Extreme_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_line_Extreme_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_Extreme_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_line_Extreme_setting_edit_ylim);
                return; 
            end
        end
        
  
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_line_Extreme_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(1)<T(1) || setting.tlim.param(1)>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_Extreme_setting_edit_tlim);
                return;
            end
            
            if setting.tlim.param(2)<T(1) || setting.tlim.param(2)>T(end)
               set(h.monitorbox,'foregroundcolor','r','string',['>>Specify t lim in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_line_Extreme_setting_edit_tlim);
                return;
            end
            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_Extreme_setting_edit_coarse);
                return;
            end
            
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_Extreme_setting_edit_coarse);
                return;
            end
        end
        
        if setting.ampliref.check==1
            if length(setting.ampliref.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference'])
                uicontrol(h.pp_line_Extreme_setting_edit_ampliref);
                return;
            end
            
            if  setting.ampliref.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference that is > 0'])
                uicontrol(h.pp_line_Extreme_setting_edit_ampliref);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_line_plot_Extreme(simuldata,setting,axesfig,h);
    end

    function callback_quant_buoy_calculate(hObj,eventdata,h)
        buoy.x=get(h.pp_quant_buoy_x_edit,'userdata');
        buoy.y=get(h.pp_quant_buoy_y_edit,'userdata');
        buoy.tinterv=get(h.pp_quant_buoy_tinterv_edit,'userdata');
        buoy.combine=get(h.pp_quant_buoy_cb_combine_buoys,'value');
        buoy.savedata=get(h.pp_quant_buoy_cb_savedata,'value');
        simuldata=get(h.pp_project_load_data,'userdata');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if isempty(buoy.x)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify position of buoy(s) in x axis'])
            uicontrol(h.pp_quant_buoy_x_edit);
            return;
        end
        
         if any(buoy.x<X(1)) || any(buoy.x>X(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
            uicontrol(h.pp_quant_buoy_x_edit);
            return;
        end
        
        if isempty(buoy.y)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify position of buoy(s) in y axis'])
            uicontrol(h.pp_quant_buoy_y_edit);
            return;
        end
        
        if any(buoy.y<Y(1)) || any(buoy.y>Y(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
            uicontrol(h.pp_quant_buoy_y_edit);
            return;
        end
        
        if length(buoy.x)~=length(buoy.y)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Number of input in x and y has to be the same'])
            uicontrol(h.pp_quant_buoy_x_edit);
            return;
        end
        
        if length(buoy.tinterv)~=2 || buoy.tinterv(1)>=buoy.tinterv(2)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a time interval (t start;  t end)'])
            uicontrol(h.pp_quant_buoy_tinterv_edit);
            return;
        end
        
        if any(buoy.tinterv<T(1)) || any(buoy.tinterv>T(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time interval in [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
            uicontrol(h.pp_quant_buoy_tinterv_edit);
            return;
        end
        
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_quant_buoy(simuldata,buoy,h);
        
        
    end

    
    function callback_anim_edit_bathy_scale(hObj,eventdata)
       param=str2num(get(hObj,'string'));
       set(hObj,'userdata',param)
    end

    
    function callback_anim_cb_bathy(hObj,eventdata,hobs)
       Id=get(hObj,'value');
       if Id==1
       set(hobs,'enable','on')
       else
       set(hobs,'enable','off')    
       end
    end

    function callback_pp_plot_anim_density(hObj,eventdata,h)
        if get(h.pp_anim_pause,'userdata')==1
            set(h.monitorbox,'foregroundcolor','r','string','>> Please press stop button.')
            return;
        end
        
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_anim_density_popup_var,'value');
        setting.view.check=get(h.pp_anim_density_setting_cb_view,'value');
        setting.view.param=get(h.pp_anim_density_setting_edit_view,'userdata');
        setting.clim.check=get(h.pp_anim_density_setting_cb_clim,'value');
        setting.clim.param=get(h.pp_anim_density_setting_edit_clim,'userdata');
        setting.xlim.check=get(h.pp_anim_density_setting_cb_xlim,'value');
        setting.xlim.param=get(h.pp_anim_density_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.pp_anim_density_setting_cb_ylim,'value');
        setting.ylim.param=get(h.pp_anim_density_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.pp_anim_density_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_anim_density_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_anim_density_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_anim_density_setting_edit_coarse,'userdata');
        setting.ampliref.check=get(h.pp_anim_density_setting_cb_ampliref,'value');
        setting.ampliref.param=get(h.pp_anim_density_setting_edit_ampliref,'userdata');
        setting.levelquiver.check=get(h.pp_anim_density_setting_cb_level,'value');
        setting.levelquiver.param=get(h.pp_anim_density_setting_edit_level,'userdata');
        setting.saveanim.check=get(h.pp_anim_density_setting_cb_saveanim,'value');
        setting.gifset.param=get(h.pp_anim_density_setting_edit_gifset,'userdata');
        setting.level.eta_check=get(h.pp_anim_density_cb_level_eta,'value');
        setting.level.phi_check=get(h.pp_anim_density_cb_level_phi,'value');
        setting.level.quiver_check=get(h.pp_anim_density_cb_level_quiver,'value');
        setting.level.quiver_check=get(h.pp_anim_density_cb_level_quiver,'value');
        setting.level.extremecrest=get(h.pp_anim_density_cb_level_extremeCrest,'value');
        setting.level.extremetrough=get(h.pp_anim_density_cb_level_extremeTrough,'value');
        setting.bathy.cb=get(h.pp_anim_density_cb_bathy,'value');
        setting.bathy.scale=get(h.pp_anim_density_edit_bathyScale,'userdata');
        
        
        colorfig=cellstr(get(h.pp_anim_density_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.pp_anim_density_setting_popup_colormap,'value'));
        axesfig=h.pp_anim_density_axes;
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
         if setting.var>2
             
             if ID_model_phiform==1
                 if ~isfield(simuldata.output,'phi')
                     set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                     uicontrol(h.pp_anim_density_popup_var);
                     return;
                 end
                 if length(simuldata.output.phi(:,1,1))==1
                     set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                     uicontrol(h.pp_anim_density_popup_var);
                     return;
                 end
             else
                 if ~isfield(simuldata.output,'u')
                     set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                     uicontrol(h.pp_anim_density_popup_var);
                     return;
                 end
                 if length(simuldata.output.u(:,1,1))==1
                     set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                     uicontrol(h.pp_anim_density_popup_var);
                     return;
                 end
             end
        end
        
        if setting.bathy.cb==1 && (isempty(setting.bathy.scale)||length(setting.bathy.scale)>1||setting.bathy.scale<0)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a scale for the bathymetry plot')
            uicontrol(h.pp_anim_density_edit_bathyScale);
            return;   
        end
        
        if setting.view.check==1
            if isempty(setting.view.param)||length(setting.view.param)>2||length(setting.view.param)<1
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_anim_density_setting_edit_view);
                return;
            end
            if length(setting.view.param)==1 && (setting.view.param>3||setting.view.param<2)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.pp_anim_density_setting_edit_view);
                return;
            end
        end
        if setting.clim.check==1
            if isempty(setting.clim.param)||length(setting.clim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a colorbar axes limit [zmin;zmax]'])
                uicontrol(h.pp_anim_density_setting_edit_clim);
                return;
            end
        end
        
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.pp_anim_density_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(1)<X(1) || setting.xlim.param(1)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_anim_density_setting_edit_xlim);
                return; 
            end
            
            if setting.xlim.param(2)<X(1) || setting.xlim.param(2)>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_anim_density_setting_edit_xlim);
                return; 
            end
            
        end
        if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.pp_anim_density_setting_edit_ylim);
                return;
            end
            
            if setting.ylim.param(1)<Y(1) || setting.ylim.param(1)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_anim_density_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Y(1) || setting.ylim.param(2)>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_anim_density_setting_edit_ylim);
                return; 
            end
        end
        
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a time limit [tmin;tmax]'])
                uicontrol(h.pp_anim_density_setting_edit_tlim);
                return;
            end
           
            if  any(setting.tlim.param<T(1)-0.0001) || any(setting.tlim.param>T(end)+0.0001)
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in the interval ['...
                ,num2str(T(1)),';',num2str(T(end)),'] (s)'])
            uicontrol(h.pp_anim_density_setting_edit_tlim);
            return;
            end
        end
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_anim_density_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_anim_density_setting_edit_coarse);
                return;
            end
        end
        
        if setting.ampliref.check==1
            if length(setting.ampliref.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference'])
                uicontrol(h.pp_anim_density_setting_edit_ampliref);
                return;
            end
            if  setting.ampliref.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify an amplitude reference that is > 0'])
                uicontrol(h.pp_anim_density_setting_edit_ampliref);
                return;
            end
        end
        
        if setting.levelquiver.check==1
            if isempty(setting.levelquiver.param)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify level/quiver option'])
                uicontrol(h.pp_anim_density_setting_edit_level);
                return;
            end
        end
        
        if setting.saveanim.check==1
            if isempty(setting.gifset.param)||length(setting.gifset.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify gif parameters: [delay factor;number of loop]'])
                uicontrol(h.pp_anim_density_setting_edit_gifset);
                return;
            end
        end
        
        if setting.clim.check==0
            X=simuldata.output.X;Y=simuldata.output.Y;
            T=simuldata.output.time;
            if setting.tlim.check==1
                indt1=funC_closest(T,setting.tlim.param(1));
                indt2=funC_closest(T,setting.tlim.param(2));
            else
                indt1=1;indt2=length(T);
            end
            
            if setting.xlim.check==1
                indx1=funC_closest(X,setting.xlim.param(1));
                indx2=funC_closest(X,setting.xlim.param(2));
            else
                indx1=1;indx2=length(X);
            end
            
            if setting.ylim.check==1
                indy1=funC_closest(Y,setting.ylim.param(1));
                indy2=funC_closest(Y,setting.ylim.param(2));
            else
                indy1=1;indy2=length(Y);
            end
            
            if setting.var==1
                var=simuldata.output.eta(indt1:indt2,indy1:indy2,indx1:indx2);
                cclim=[min(min(min(var))) max(max(max(var)))];
                set(h.pp_anim_density_setting_edit_clim,'userdata',cclim)
                set(h.pp_anim_density_setting_edit_clim,'string',...
                    [num2str(roundn(min(min(min(var))),-3)),';',...
                    num2str(roundn(max(max(max(var))),-3))])
            elseif setting.var==3
%                 var=simuldata.output.phi(indt1:indt2,indy1:indy2,indx1:indx2);
%                 cclim=[min(min(min(var))) max(max(max(var)))];
%                 set(h.pp_anim_density_setting_edit_clim,'userdata',cclim)
%                 set(h.pp_anim_density_setting_edit_clim,'string',...
%                     [num2str(roundn(min(min(min(var))),-3)),';',...
%                     num2str(roundn(max(max(max(var))),-3))])
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
        funP_anim_density(simuldata,setting,axesfig,h);
        
    end

    function callback_pp_plot_anim_line(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        if get(h.pp_anim_line_pause,'userdata')==1
            set(h.monitorbox,'foregroundcolor','r','string','>> Please press stop button.')
            return;
        end
        
        simuldata=get(h.pp_project_load_data,'userdata');
        setting.var=get(h.pp_anim_line_popup_var,'value');
        
        setting.zlim.check=get(h.pp_anim_line_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_anim_line_setting_edit_zlim,'userdata');
        setting.spatlim.check=get(h.pp_anim_line_setting_cb_spatlim,'value');
        setting.spatlim.param=get(h.pp_anim_line_setting_edit_spatlim,'userdata');
        setting.tlim.check=get(h.pp_anim_line_setting_cb_tlim,'value');
        setting.tlim.param=get(h.pp_anim_line_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.pp_anim_line_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_anim_line_setting_edit_coarse,'userdata');
        setting.saveanim.check=get(h.pp_anim_line_setting_cb_saveanim,'value');
        setting.gifset.param=get(h.pp_anim_line_setting_edit_gifset,'userdata');
        axesfig=h.pp_anim_line_axes;
        spatsnap.x_check=get(h.pp_anim_line_x_cb,'value');
        spatsnap.x_snap=get(h.pp_anim_line_x_edit,'userdata');
        spatsnap.y_check=get(h.pp_anim_line_y_cb,'value');
        spatsnap.y_snap=get(h.pp_anim_line_y_edit,'userdata');
        setting.MTCcheck=get(h.pp_anim_line_MTC_cb,'value');
        setting.MTTcheck=get(h.pp_anim_line_MTT_cb,'value');
        setting.bathy.cb=get(h.pp_anim_line_cb_bathy,'value');
        setting.bathy.scale=get(h.pp_anim_line_edit_bathyScale,'userdata');
        
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
        if setting.var>1
            if ID_model_phiform==1
                if ~isfield(simuldata.output,'phi')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_anim_line_popup_var);
                    return;
                end
                if length(simuldata.output.phi(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_anim_line_popup_var);
                    return;
                end
            else
                if ~isfield(simuldata.output,'u')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_anim_line_popup_var);
                    return;
                end
                if length(simuldata.output.u(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_anim_line_popup_var);
                    return;
                end
            end
        end
        
         X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if spatsnap.x_check==0 && spatsnap.y_check==0
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x or y axis')
            uicontrol(h.pp_anim_line_x_edit);
            return;
        end
        
        if spatsnap.x_check==1
            if isempty(spatsnap.x_snap)||length(spatsnap.x_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in x axis')
                uicontrol(h.pp_anim_line_x_edit);
                return;
            end
            
             if spatsnap.x_snap<X(1) || spatsnap.x_snap>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                uicontrol(h.pp_anim_line_x_edit);
                return; 
            end
            
        end
        if spatsnap.y_check==1
            if isempty(spatsnap.y_snap)||length(spatsnap.y_snap)~=1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a position in y axis')
                uicontrol(h.pp_anim_line_y_edit);
                return;
            end
            
            if spatsnap.y_snap<Y(1) || spatsnap.y_snap>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                uicontrol(h.pp_anim_line_y_edit);
                return; 
            end
        end
        
         if setting.bathy.cb==1 && (isempty(setting.bathy.scale)||length(setting.bathy.scale)>1||setting.bathy.scale<0)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a scale for the bathymetry plot')
            uicontrol(h.pp_anim_line_edit_bathyScale);
            return;   
        end
        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_anim_line_setting_edit_zlim);
                return;
            end
        end
        
        if setting.spatlim.check==1
            if spatsnap.y_check==1
                
                          
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2 || setting.spatlim.param(1)>setting.spatlim.param(2)
                    set(h.monitorbox,'foregroundcolor','r','string','>>Specify x axis limit [xmin;xmax]')
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<X(1) || setting.spatlim.param(1)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<X(1) || setting.spatlim.param(2)>X(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
            else
                
                if isempty(setting.spatlim.param)||length(setting.spatlim.param)~=2 || setting.spatlim.param(1)>setting.spatlim.param(2)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(1)<Y(1) || setting.spatlim.param(1)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
                
                if setting.spatlim.param(2)<Y(1) || setting.spatlim.param(2)>Y(end)
                    set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
                    uicontrol(h.pp_anim_line_setting_edit_spatlim);
                    return;
                end
            end
        end
        
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time axis limit [tmin;tmax]'])
                uicontrol(h.pp_lanim_line_setting_edit_tlim);
                return;
            end
            
            if any(setting.tlim.param<T(1)-0.0001) || any(setting.tlim.param>T(end)+0.0001)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time limit in the interval [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
                uicontrol(h.pp_anim_line_setting_edit_tlim);
                return;
            end

            
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_anim_line_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_anim_line_setting_edit_coarse);
                return;
            end
        end
        
        if setting.saveanim.check==1
            if isempty(setting.gifset.param)||length(setting.gifset.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify gif parameters: [delay factor;number of loop]'])
                uicontrol(h.pp_anim_line_setting_edit_gifset);
                return;
            end
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        funP_anim_line(spatsnap,simuldata,setting,axesfig,h);
    end

    function callback_resume(hObj,eventdata,axesnow,href)
        if get(href,'userdata')==1;
            uiresume(gcbf)
            set(href,'userdata',0);
        end
    end

    function callback_pause(hObj,eventdata)
        set(hObj,'userdata',1);
    end

    function callback_stop(hObj,eventdata,href)
        if get(href,'userdata')==1;
            uiresume(gcbf);
            set(href,'userdata',0);
        end
        
        set(hObj,'userdata',1);
    end



    function callback_zoom_in(hObj,~,h)
        if get(hObj,'value')
            bdown = get(h.fig,'WindowButtonDownFcn');
            bmove = get(h.fig,'WindowButtonMotionFcn');
            
            if ~isempty(bdown) %if not empty
                %bdown == 'localModeWindowButtonDownFcn'
                if iscell(bdown)
                    temp = func2str(bdown{1,1});
                    if ~strcmpi(temp(1,1:9),'localMode')
                        set(h.fig,'WindowButtonDownFcn','');
                    end
                else
                    set(h.fig,'WindowButtonDownFcn','');
                end
            end
            
            if ~isempty(bmove) %if not empty
                %bdown == 'localModeWindowButtonDownFcn'
                if iscell(bmove)
                    temp = func2str(bmove{1,1});
                    if ~strcmpi(temp(1,1:9),'localMode')
                        set(h.fig,'WindowButtonMotionFcn','');
                    end
                else
                    set(h.fig,'WindowButtonMotionFcn','');
                end
            end
            
            datacursormode off;
            zoom on;
        else
            zoom off;
        end
    end

    function callback_zoom_out(hObj,~,h)
        if get(hObj,'value')
            gca;
            zoom off; datacursormode off;
            zoom out;
        end
    end

    function callback_rotate(hObj,~,h)
        if get(hObj,'value')==1
            rotate3d on;
        else
            rotate3d off;
        end
    end


    function callback_pan(hObj,~,h)
        if get(hObj,'value')==1
            pan on;
        else
            pan off;
        end
    end

    function callback_datacursor(hObj,~,h)
        zoom off;
        if get(hObj,'value')==1
            datacursormode on;
        else
            datacursormode off;
        end
    end

    function callback_preview_wave_spectrum(hObj,~,h)
        callback_preview_wave(h.preview.wave_popup_var,[],h);
    end

    function callback_pp_validation_load_meas_data(hObj,~,h)
       projdir=get(h.pp_project_popup_projdir,'string');
       savename=get(h.pp_project_edit_name,'string');
       workdir=[projdir,'\',savename,'\'];
       [file_name,directory]=uigetfile([workdir,'/','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
        set(h.monitorbox,'foregroundcolor','k','string',['>>Uploading data...'])
        if directory~=0
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
          if isstruct(my_data)==0
             myBuoy_X=my_data(1,2:end); 
             myBuoy_Y=my_data(2,2:end); 
             myBuoy_T=my_data(3:end,1);
             myBuoy_data=my_data(3:end,2:end);
             if length(myBuoy_X)~=length(myBuoy_Y)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return;
             end
             if length(myBuoy_data(1,:))~=length(myBuoy_X)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return; 
             end
             if length(myBuoy_data(:,1))~=length(myBuoy_T)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return; 
             end
             NX=length(myBuoy_X);
             if NX~=0
                 if NX~=1
                     for i=1:NX-1
                         Buoy_X{i}=[num2str(myBuoy_X(i)),';',];
                         Buoy_Y{i}=[num2str(myBuoy_Y(i)),';',];
                     end
                 end
                 Buoy_X{NX}=num2str(myBuoy_X(NX));
                 Buoy_Y{NX}=num2str(myBuoy_Y(NX));
                 set(h.pp_validation_buoy_x_edit,'String',[Buoy_X{1:NX}]);
                 set(h.pp_validation_buoy_x_edit,'userdata',myBuoy_X);
                 set(h.pp_validation_buoy_y_edit,'String',[Buoy_Y{1:NX}]);
                 set(h.pp_validation_buoy_y_edit,'userdata',myBuoy_Y);
                 
                 simuldata=get(h.pp_project_load_data,'userdata');
                 if isfield(simuldata,'output')
                 T=simuldata.output.time;
                 Tmin=min(T(1),myBuoy_T(1));
                 Tmax=min(T(end),myBuoy_T(end));
                 else
                 Tmin=myBuoy_T(1);
                 Tmax=myBuoy_T(end);   
                 end
                 set(h.pp_validation_buoy_time_edit,'String',...
                     [num2str(Tmin),';',num2str(Tmax)],...
                     'userdata',[Tmin;Tmax]);
             end
             set(h.pp_validation_buoy_number_popup,'string',{'1'},'value',1);
               
            set(hObj,'Userdata',my_data); 
            set(h.monitorbox,'foregroundcolor','k','string',['>>',file_name,' has been loaded'])
          else
            set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
          end
        else
          set(h.monitorbox,'foregroundcolor','k','string',['>>']) 
        end
    end

    function callback_pp_validation_buoy_cb_timeshift(hObj,~,h)
     Id=get(hObj,'value');
     if Id==1
        set(h.pp_validation_buoy_cb_timeshift_def,'enable','on','value',1)
        set(h.pp_validation_buoy_cb_timeshift_dt,'enable','on','value',0)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','off')
     else
        set(h.pp_validation_buoy_cb_timeshift_def,'enable','off','value',0)
        set(h.pp_validation_buoy_cb_timeshift_dt,'enable','off','value',0)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','off')
     end
    end
    function callback_pp_validation_buoy_cb_timeshift_dt(hObj,~,h)
      Id=get(hObj,'value');
     if Id==1
        set(h.pp_validation_buoy_cb_timeshift_def,'enable','on','value',0)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','on')
     else
        set(h.pp_validation_buoy_cb_timeshift_def,'enable','on','value',1)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','off')
     end   
    end

   function callback_pp_validation_buoy_edit_timeshift_dt(hObj,~,h)
     param=str2num(get(hObj,'string'));
     set(hObj,'userdata',param)
     set(h.monitorbox,'foregroundcolor','k','string',['>>']) 
     
       if isempty(param)||length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a number how many dt to be shifted.'])
            uicontrol(hObj);
            return;
        end
   end

    function callback_pp_validation_buoy_cb_timeshift_def(hObj,~,h)
      Id=get(hObj,'value');
     if Id==1
        set(h.pp_validation_buoy_cb_timeshift_dt,'enable','on','value',0)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','off')
     else
        set(h.pp_validation_buoy_cb_timeshift_dt,'enable','on','value',1)
        set(h.pp_validation_buoy_edit_timeshift_dt,'enable','on')
     end     
    end

    function callback_pp_validation_buoy_edit_xy(hObj,~,h)
      param=str2num(get(hObj,'string'));
      set(hObj,'userdata',param)
     set(h.monitorbox,'foregroundcolor','k','string',['>>']) 
     
        if isempty(param)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify buoy positions'])
            uicontrol(hObj);
            return;
        end  
    end

    function callback_pp_validation_buoy_plot(hObj,~,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        measdata=get(h.pp_validation_buoy_measdata_button,'userdata');
        setting.var=get(h.pp_validation_buoy_popup_var,'value');
        setting.cb.signal=get(h.pp_validation_buoy_cb_signal,'value');
        setting.cb.spectrum=get(h.pp_validation_buoy_cb_spectrum,'value');
        setting.cb.spectrum_var=get(h.pp_validation_buoy_popup_spectrum_var,'value');
        setting.cb.timeshift=get(h.pp_validation_buoy_cb_timeshift,'value');
        setting.cb.timeshift_def=get(h.pp_validation_buoy_cb_timeshift_def,'value');
        setting.cb.timeshift_dt=get(h.pp_validation_buoy_cb_timeshift_dt,'value');
        setting.cb.timeshift_dt_edit=get(h.pp_validation_buoy_edit_timeshift_dt,'userdata');
        
        setting.zlim.check=get(h.pp_validation_buoy_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_validation_buoy_setting_edit_zlim,'userdata');
        setting.horznlim.check=get(h.pp_validation_buoy_setting_cb_horznlim,'value');
        setting.horznlim.param=get(h.pp_validation_buoy_setting_edit_horznlim,'userdata');
        setting.coarse.check=get(h.pp_validation_buoy_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_validation_buoy_setting_edit_coarse,'userdata');
        setting.spsmooth.check=get(h.pp_validation_buoy_setting_cb_spsmooth,'value');
        setting.spsmooth.param=get(h.pp_validation_buoy_setting_edit_spsmooth,'userdata');
        
        setting.savefig.check=get(h.pp_validation_buoy_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_validation_buoy_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_validation_buoy_setting_popup_savefig,'value'));
        axesfig=h.pp_validation_buoy_axes;
        spatsnap.x_snap=get(h.pp_validation_buoy_x_edit,'userdata');
        spatsnap.y_snap=get(h.pp_validation_buoy_y_edit,'userdata');
        spatsnap.tinterv=get(h.pp_validation_buoy_time_edit,'userdata');
             
        setting.savedat=get(h.pp_validation_buoy_savedata,'value');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
          
        if setting.var==2
            if ID_model_phiform==1
                if ~isfield(simuldata.output,'phi')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_validation_buoy_popup_var);
                    return;
                end
                if length(simuldata.output.phi(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave potential data')
                    uicontrol(h.pp_validation_buoy_popup_var);
                    return;
                end
            else
                if ~isfield(simuldata.output,'u')
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_validation_buoy_popup_var);
                    return;
                end
                if length(simuldata.output.u(:,1,1))==1
                    set(h.monitorbox,'foregroundcolor','r','string','>> There is no wave velocity data')
                    uicontrol(h.pp_validation_buoy_popup_var);
                    return;
                end
            end
        end
        
        if isempty(measdata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify measurement data')
            uicontrol(h.pp_validation_buoy_measdata_button);
            return;
        end
        
        X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        
        if isempty(spatsnap.x_snap)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify buoy positions in x axis, x:[x1;x2;x2]')
            uicontrol(h.pp_validation_buoy_x_edit);
            return;
        end
        
        if any(spatsnap.x_snap<X(1)) || any(spatsnap.x_snap>X(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
            uicontrol(h.pp_validation_buoy_x_edit);
            return;
        end
        
        if isempty(spatsnap.y_snap)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify buoy positions in y axis, y:[y1;y2;y3]')
            uicontrol(h.pp_validation_buoy_y_edit);
            return;
        end
        
        if any(spatsnap.y_snap<Y(1)) || any(spatsnap.y_snap>Y(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
            uicontrol(h.pp_validation_buoy_y_edit);
            return;
        end
        
        if length(spatsnap.x_snap)~=length(spatsnap.y_snap)
            set(h.monitorbox,'foregroundcolor','r','string','>> Number of Buoy positions in x and y axis must be the same')          
            uicontrol(h.pp_validation_buoy_x_edit);
            return;    
        end
        
        if length(spatsnap.x_snap)>length(measdata(1,2:end))
            set(h.monitorbox,'foregroundcolor','r','string','>> Number of Buoy positions in x and y axis is larger than in the measurement data')          
            uicontrol(h.pp_validation_buoy_x_edit);
            return;     
        end
        
        
        if length(spatsnap.tinterv)~=2 || spatsnap.tinterv(1)>=spatsnap.tinterv(2)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a time interval, t:[tmin;tmax]')
            uicontrol(h.pp_validation_buoy_time_edit);
            return;
        end
        
        if any(spatsnap.tinterv<T(1)) || any(spatsnap.tinterv>T(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time interval in [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
            uicontrol(h.pp_validation_buoy_time_edit);
            return;
        end
        
        if setting.cb.timeshift_dt==1 && length(setting.cb.timeshift_dt_edit)~=1
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a factor of time shifting')
            uicontrol(h.pp_validation_buoy_edit_timeshift_dt);
            return;
        end

        
        if setting.zlim.check==1
            if isempty(setting.zlim.param)||length(setting.zlim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a vertical axes limit [zmin;zmax]'])
                uicontrol(h.pp_validation_buoy_setting_edit_zlim);
                return;
            end
        end
        
        if setting.horznlim.check==1
            if isempty(setting.horznlim.param)||length(setting.horznlim.param)~=2 || setting.horznlim.param(1)>=setting.horznlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a horizontal axes limit, [min;max]'])
                uicontrol(h.pp_validation_buoy_setting_edit_horznlim);
                return;
            end
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.pp_line_statistic_setting_edit_coarse);
                return;
            end
        end
        
        
        if setting.spsmooth.check==1
            if length(setting.spsmooth.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a smooth factor for spectrum'])
                uicontrol(h.pp_validation_buoy_setting_cb_spsmooth);
                return;
            end
        end
        
        Nx=length(spatsnap.x_snap);
        for ii=1:Nx
            strBuoy{ii}=num2str(ii);
        end
        strBuoy{Nx}=num2str(Nx);
        set(h.pp_validation_buoy_number_popup,'string',strBuoy,'value',1);
        set(h.pp_validation_buoy_X1disp_edit,'string',num2str(spatsnap.x_snap(1)));
        set(h.pp_validation_buoy_Y1disp_edit,'string',num2str(spatsnap.y_snap(1)));
       
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
    
        funP_validation_plot_buoy(spatsnap,simuldata,measdata,setting,axesfig,h.fig,1);
    end

    function callback_pp_validation_buoy_plot_popup(hObj,~,h)
    Idplot=get(hObj,'value');
    
    set(h.monitorbox,'foregroundcolor','k','string','>>')
        simuldata=get(h.pp_project_load_data,'userdata');
        measdata=get(h.pp_validation_buoy_measdata_button,'userdata');
        setting.var=get(h.pp_validation_buoy_popup_var,'value');
        setting.cb.signal=get(h.pp_validation_buoy_cb_signal,'value');
        setting.cb.spectrum=get(h.pp_validation_buoy_cb_spectrum,'value');
        setting.cb.spectrum_var=get(h.pp_validation_buoy_popup_spectrum_var,'value');
        setting.cb.timeshift=get(h.pp_validation_buoy_cb_timeshift,'value');
        setting.cb.timeshift_def=get(h.pp_validation_buoy_cb_timeshift_def,'value');
        setting.cb.timeshift_dt=get(h.pp_validation_buoy_cb_timeshift_dt,'value');
        setting.cb.timeshift_dt_edit=get(h.pp_validation_buoy_edit_timeshift_dt,'userdata');
        
        setting.zlim.check=get(h.pp_validation_buoy_setting_cb_zlim,'value');
        setting.zlim.param=get(h.pp_validation_buoy_setting_edit_zlim,'userdata');
        setting.horznlim.check=get(h.pp_validation_buoy_setting_cb_horznlim,'value');
        setting.horznlim.param=get(h.pp_validation_buoy_setting_edit_horznlim,'userdata');
        setting.coarse.check=get(h.pp_validation_buoy_setting_cb_coarse,'value');
        setting.coarse.param=get(h.pp_validation_buoy_setting_edit_coarse,'userdata');
        setting.spsmooth.check=get(h.pp_validation_buoy_setting_cb_spsmooth,'value');
        setting.spsmooth.param=get(h.pp_validation_buoy_setting_edit_spsmooth,'userdata');
        
        setting.savefig.check=get(h.pp_validation_buoy_setting_cb_savefig,'value');
        format=cellstr(get(h.pp_validation_buoy_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.pp_validation_buoy_setting_popup_savefig,'value'));
        axesfig=h.pp_validation_buoy_axes;
        spatsnap.x_snap=get(h.pp_validation_buoy_x_edit,'userdata');
        spatsnap.y_snap=get(h.pp_validation_buoy_y_edit,'userdata');
        spatsnap.tinterv=get(h.pp_validation_buoy_time_edit,'userdata');
             
        setting.savedat=get(h.pp_validation_buoy_savedata,'value');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        
        set(h.pp_validation_buoy_X1disp_edit,'string',num2str(spatsnap.x_snap(Idplot)));
        set(h.pp_validation_buoy_Y1disp_edit,'string',num2str(spatsnap.y_snap(Idplot)));
        
    funP_validation_plot_buoy(spatsnap,simuldata,measdata,setting,axesfig,h.fig,Idplot);    
    end


    function callback_pp_validation_quant_load_meas_data(hObj,~,h)
       projdir=get(h.pp_project_popup_projdir,'string');
       savename=get(h.pp_project_edit_name,'string');
       workdir=[projdir,'\',savename,'\'];
       [file_name,directory]=uigetfile([workdir,'/','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
        set(h.monitorbox,'foregroundcolor','k','string',['>>Uploading data...'])
        if directory~=0
            temp=load([directory,file_name]);
            if isstruct(temp)
                namevar = fieldnames(temp);
                my_data=temp.(namevar{1});
            else
                my_data=temp;
            end
            clearvars temp;
          if isstruct(my_data)==0
             myBuoy_X=my_data(1,2:end); 
             myBuoy_Y=my_data(2,2:end); 
             myBuoy_T=my_data(3:end,1);
             myBuoy_data=my_data(3:end,2:end);
             if length(myBuoy_X)~=length(myBuoy_Y)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return;
             end
             if length(myBuoy_data(1,:))~=length(myBuoy_X)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return; 
             end
             if length(myBuoy_data(:,1))~=length(myBuoy_T)
                 set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
                 return; 
             end
             NX=length(myBuoy_X);
             if NX~=0
                 if NX~=1
                     for i=1:NX-1
                         Buoy_X{i}=[num2str(myBuoy_X(i)),';',];
                         Buoy_Y{i}=[num2str(myBuoy_Y(i)),';',];
                     end
                 end
                 Buoy_X{NX}=num2str(myBuoy_X(NX));
                 Buoy_Y{NX}=num2str(myBuoy_Y(NX));
                 set(h.pp_validation_quant_buoy_x_edit,'String',[Buoy_X{1:NX}]);
                 set(h.pp_validation_quant_buoy_x_edit,'userdata',myBuoy_X);
                 set(h.pp_validation_quant_buoy_y_edit,'String',[Buoy_Y{1:NX}]);
                 set(h.pp_validation_quant_buoy_y_edit,'userdata',myBuoy_Y);
                 
                 simuldata=get(h.pp_project_load_data,'userdata');
                 if isfield(simuldata,'output')
                 T=simuldata.output.time;
                 Tmin=min(T(1),myBuoy_T(1));
                 Tmax=min(T(end),myBuoy_T(end));
                 else
                 Tmin=myBuoy_T(1);
                 Tmax=myBuoy_T(end);   
                 end
                 set(h.pp_validation_quant_buoy_time_edit,'String',...
                     [num2str(Tmin),';',num2str(Tmax)],...
                     'userdata',[Tmin;Tmax]);
             end
               
            set(hObj,'Userdata',my_data); 
            set(h.monitorbox,'foregroundcolor','k','string',['>>',file_name,' has been loaded'])
          else
            set(h.monitorbox,'foregroundcolor','r','string',['>>wrong input file!'])    
          end
        else
          set(h.monitorbox,'foregroundcolor','k','string',['>>']) 
        end
    end

    function callback_pp_validation_quant_buoy_cb_timeshift(hObj,~,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.pp_validation_quant_buoy_cb_timeshift_def,'enable','on','value',1)
            set(h.pp_validation_quant_buoy_cb_timeshift_dt,'enable','on','value',0)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','off')
        else
            set(h.pp_validation_quant_buoy_cb_timeshift_def,'enable','off','value',0)
            set(h.pp_validation_quant_buoy_cb_timeshift_dt,'enable','off','value',0)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','off')
        end
    end
    function callback_pp_validation_quant_buoy_cb_timeshift_dt(hObj,~,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.pp_validation_quant_buoy_cb_timeshift_def,'enable','on','value',0)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','on')
        else
            set(h.pp_validation_quant_buoy_cb_timeshift_def,'enable','on','value',1)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','off')
        end
    end

    function callback_pp_validation_quant_buoy_edit_timeshift_dt(hObj,~,h)
        param=str2num(get(hObj,'string'));
        if isempty(param)||length(param)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a number how many dt to be shifted.'])
            uicontrol(hObj);
            return;
        else
            set(hObj,'userdata',param)
            set(h.monitorbox,'foregroundcolor','k','string',['>>'])
        end
    end

    function callback_pp_validation_quant_buoy_cb_timeshift_def(hObj,~,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.pp_validation_quant_buoy_cb_timeshift_dt,'enable','on','value',0)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','off')
        else
            set(h.pp_validation_quant_buoy_cb_timeshift_dt,'enable','on','value',1)
            set(h.pp_validation_quant_buoy_edit_timeshift_dt,'enable','on')
        end
    end

    function callback_pp_validation_quant_buoy_calc(hObj,~,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        buoy.x=get(h.pp_validation_quant_buoy_x_edit,'userdata');
        buoy.y=get(h.pp_validation_quant_buoy_y_edit,'userdata');
        buoy.tinterv=get(h.pp_validation_quant_buoy_time_edit,'userdata');
        buoy.savedata=get(h.pp_validation_quant_buoy_cb_savedata,'value');
        setting.cb.timeshift=get(h.pp_validation_quant_buoy_cb_timeshift,'value');
        setting.cb.timeshift_def=get(h.pp_validation_quant_buoy_cb_timeshift_def,'value');
        setting.cb.timeshift_dt=get(h.pp_validation_quant_buoy_cb_timeshift_dt,'value');
        setting.cb.timeshift_dt_edit=get(h.pp_validation_quant_buoy_edit_timeshift_dt,'userdata');
        
        simuldata=get(h.pp_project_load_data,'userdata');
        filename=get(h.pp_project_edit_name,'string');
        projdir=get(h.pp_project_popup_projdir,'string');
        measdata=get(h.pp_validation_quant_buoy_measdata_button,'userdata');
        
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.pp_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.pp_project_popup_projdir);
            return;
        end
        
        if isempty(simuldata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a simulation data')
            uicontrol(h.pp_project_load_data);
            return;
        end
        
          
        if isempty(measdata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify measurement data')
            uicontrol(h.pp_validation_quant_buoy_measdata_button);
            return;
        end
        
         X=simuldata.output.X;Y=simuldata.output.Y;
        T=simuldata.output.time;
        
        if isempty(buoy.x)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify position of buoy(s) in x axis'])
            uicontrol(h.pp_validation_quant_buoy_x_edit);
            return;
        end
        
        if any(buoy.x<X(1)) || any(buoy.x>X(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in x axis in the interval [',num2str(X(1)),';',num2str(X(end)),'] (m)']);
            uicontrol(h.pp_validation_quant_buoy_x_edit);
            return;
        end
        
        if isempty(buoy.y)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify position of buoy(s) in y axis'])
            uicontrol(h.pp_validation_quant_buoy_y_edit);
            return;
        end
        
        if any(buoy.y<Y(1)) || any(buoy.y>Y(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a position or multiple position in y axis in the interval [',num2str(Y(1)),';',num2str(Y(end)),'] (m)']);
            uicontrol(h.pp_validation_quant_buoy_y_edit);
            return;
        end
        
        if length(buoy.x)~=length(buoy.y)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Number of input in x and y should be same'])
            uicontrol(h.pp_validation_quant_buoy_x_edit);
            return;
        end
        
         if length(buoy.x)>length(measdata(1,2:end))
            set(h.monitorbox,'foregroundcolor','r','string','>> Number of Buoy positions in x and y axis is larger than in the measurement data')          
            uicontrol(h.pp_validation_buoy_x_edit);
            return;     
        end
       
        if length(buoy.tinterv)~=2 || buoy.tinterv(1)>=buoy.tinterv(2)
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a time interval (t start;  t end)'])
            uicontrol(h.pp_validation_quant_buoy_time_edit);
            return;
        end
        
        if any(buoy.tinterv<T(1)) || any(buoy.tinterv>T(end))
            set(h.monitorbox,'foregroundcolor','r','string',['>>Specify time interval in [',num2str(T(1)),';',num2str(T(end)),'] (s)']);
            uicontrol(h.pp_validation_quant_buoy_time_edit);
            return;
        end
        
        if setting.cb.timeshift_dt==1 && length(setting.cb.timeshift_dt_edit)~=1
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a factor of time shifting')
            uicontrol(h.pp_validation_quant_buoy_edit_timeshift_dt);
            return;
        end
        
        simuldata.Proj.path=h.pathnow;
        simuldata.Proj.projdir=projdir;
        simuldata.Proj.savename=filename;
        simuldata.Proj.workdir=[projdir,'\',filename,'\'];
        if ~isdir(simuldata.Proj.workdir)
            mkdir(simuldata.Proj.workdir);
        end
        
         funP_validation_quant_buoy(simuldata,measdata,buoy,setting,h);   
    end

    function callback_preview_spatial(hObj,eventdata,h)
        Id=get(hObj,'value');
        
        preproc=get(hObj,'userdata');
        input=preproc(1).input;dom=preproc(1).dom;
        set(h.monitorbox,'foregroundcolor','k','string',['>>']);
         if Id==4 && get(h.bathymetry_checkbox_friction,'value')==0
            set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
            set(hObj,'value',1);
        elseif Id==5 && strcmpi(input.wall.option,'No')
            set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
            set(hObj,'value',1);
        elseif Id==6 && strcmpi(input.wall.option,'No') 
            set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
            set(hObj,'value',1);
        elseif Id==6 && strcmpi(input.wall.option,'Yes') && dom.wall.NInfl==0
            set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
            set(hObj,'value',1);
        elseif Id==7 && strcmpi(input.wall.option,'No')
            set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
             set(hObj,'value',1);
        elseif Id==7 && strcmpi(input.wall.option,'Yes') && dom.wall.NInfl==0 
             set(h.monitorbox,'foregroundcolor','r','string',['>>The option is not available']);
             set(hObj,'value',1);
         end
        
        PreProcSpatialView(h,preproc);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function func_save_guistate(GUIinput)
        if ~isdir(GUIinput.proj.workdir)
            fc=fix(clock);
            GUIinput.proj.createddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
            GUIinput.proj.createdtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];
            mkdir(GUIinput.proj.workdir);
        else
            if ~isfield(GUIinput.proj,'createddate')
                fc=fix(clock);
                GUIinput.proj.createddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
                GUIinput.proj.createdtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];
                
            end
        end
        
        fc=fix(clock);
        GUIinput.proj.modifieddate=[num2str(fc(3)),'/',num2str(fc(2)),'/',num2str(fc(1))];
        GUIinput.proj.modifiedtime=[ num2str(fc(4)),':',num2str(fc(5)),':',num2str(fc(6))];
        
        
        if GUIinput.wave_bdy_assim.checkVal==2
        GUIinput.wave_bdy_assim.assimdata=[];
        GUIinput.wave_bdy_assim.assimdata_phi=[];
        end
        save([GUIinput.proj.workdir,'\abproj_',GUIinput.proj.name,'.mat'], 'GUIinput')
        
        if isempty(GUIinput.proj.projhist)
            GUIinput.proj.projhist = {GUIinput.proj.name, strcat(GUIinput.proj.createddate, ' - '...
                , GUIinput.proj.createdtime), strcat(GUIinput.proj.modifieddate, ' - ', ...
                GUIinput.proj.modifiedtime), GUIinput.proj.workdir, GUIinput.proj.note};
        else
            if ~isempty(GUIinput.proj.projhist{1,1})
                flagNew=0;
                for ii=1:length(GUIinput.proj.projhist(:,1))
                    if strcmpi(GUIinput.proj.projhist(ii,1),GUIinput.proj.name)==1 && strcmpi(GUIinput.proj.projhist(ii,4),GUIinput.proj.workdir)==1
                        temp = GUIinput.proj.projhist{ii,2};
                        GUIinput.proj.projhist(ii,:) =[];
                        GUIinput.proj.projhist= [{GUIinput.proj.name, temp, strcat(GUIinput.proj.modifieddate, ' - ', ...
                            GUIinput.proj.modifiedtime), GUIinput.proj.workdir, GUIinput.proj.note};GUIinput.proj.projhist];
                        clear temp
                        flagNew=1;
                        break;
                    end
                end
                if flagNew==0
                    GUIinput.proj.projhist = [{GUIinput.proj.name, strcat(GUIinput.proj.createddate, ' - '...
                        , GUIinput.proj.createdtime), strcat(GUIinput.proj.modifieddate, ' - ', ...
                        GUIinput.proj.modifiedtime), GUIinput.proj.workdir, GUIinput.proj.note};GUIinput.proj.projhist];
                end
            else
                GUIinput.proj.projhist = {GUIinput.proj.name, strcat(GUIinput.proj.createddate, ' - '...
                    , GUIinput.proj.createdtime), strcat(GUIinput.proj.modifieddate, ' - ', ...
                    GUIinput.proj.modifiedtime), GUIinput.proj.workdir, GUIinput.proj.note};
            end
        end
        
        if length(GUIinput.proj.projhist(:,1))>10
            GUIinput.proj.projhist = GUIinput.proj.projhist(1:10,:);
        end
        project = GUIinput.proj.projhist;
        userfold=getenv('LOCALAPPDATA');
        projhistfile=[userfold,'\projhistAB.mat'];
        save(projhistfile,'project');
        
    end
%%%%%%%%%%%%%%%%%%%%% function for initialization gui handles%%%%%%%%%%%%%%%%%%%%%%
    function initialization_gui(h)
        if h.flagOpenProj==0
            %%%%%%%%%%%%%Input panel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(h.project_edit_name,'string',h.projectname);
            set(h.project_popup_projdir,'string',h.projectdirectory);
            set(h.project_edit_note,'string',h.usernote);
            set(h.model_break_edit_initiation,'string','0.8','userdata',0.8);
            set(h.model_break_edit_termination,'string','0.2','userdata',0.2);
            set(h.model_break_edit_Tstar,'string','0.5','userdata',0.5);
            set(h.waveinput_influx_checkbox_ramp,'value',1);
            set(h.waveinput_influx_edit_ramp,'userdata',2,'string','2');
            set(h.waveinput_bdy_assim_popup,'value',1);
            callback_bdyassim_popup(h.waveinput_bdy_assim_popup,[],h);
            set(h.options_checkbox_intflow,'value',0,'enable','on');
            set(h.options_edit_interval_init,'userdata',[],'string','','enable','off');
            set(h.options_edit_interval_end,'userdata',[],'string','','enable','off');
            set(h.options_edit_step,'userdata',[],'string','','enable','off');
            set(h.options_checkbox_default,'value',1);
            set(h.options_edit_partition,'userdata',[],'string','','enable','off');
            set(h.options_checkbox_combine,'value',1);
            set(h.options_edit_ode_tol,'userdata',0.001,'string','0.001')
        else
            set(h.monitorbox,'foregroundcolor','k','string','>>loading a project]')
            reset_handles_gui(h);
            input_handles_gui(h);
            set(h.monitorbox,'foregroundcolor','k','string','>>')
        end
        global postproc_flag
        postproc_flag=0;
    end

    function reset_handles_postproc(h)
        set(h.pp_project_popup_projdir,'string','--')
        set(h.pp_project_edit_name,'string','');
        set(h.pp_project_edit_data,'string','')
        set(h.pp_project_load_data,'userdata',[]);
        set(h.pp_density_profile_popup_var,'value',1);
        set(h.pp_density_profile_time_edit, 'string', '', 'userdata', '');
        set(h.pp_density_profile_cb_level_eta,'value', 0, 'enable', 'on');
        set(h.pp_density_profile_cb_level_phi,'value', 0, 'enable', 'on');
        set(h.pp_density_profile_cb_level_quiver,'enable','off','value',0);
        set(h.pp_density_profile_cb_level_extremeCrest,'value',0,'enable','on');
        
        set(h.pp_density_profile_setting_cb_view, 'value', 0);
        set(h.pp_density_profile_setting_edit_view,'string', '','userdata', '','enable','off');
        set(h.pp_density_profile_setting_cb_clim, 'value', 0);
        set(h.pp_density_profile_setting_edit_clim,'string', '','userdata', '','enable','off');
        set(h.pp_density_profile_setting_cb_xlim,'value',0);
        set(h.pp_density_profile_setting_edit_xlim,'string', '','userdata', '','enable','off');
        set(h.pp_density_profile_setting_cb_ylim,'value',0);
        set(h.pp_density_profile_setting_edit_ylim,'string','','userdata','','enable','off');
        set(h.pp_density_profile_setting_cb_coarse, 'value',0);
        set(h.pp_density_profile_setting_edit_coarse,'string', '','userdata', '','enable','off');
        set(h.pp_density_profile_setting_cb_level, 'value', 0);
        set(h.pp_density_profile_setting_edit_level  , 'string','','userdata','','enable','off');
        set(h.pp_density_profile_setting_cb_savefig, 'value',0);
        set(h.pp_density_profile_setting_popup_savefig, 'value', 1, 'enable', 'off');
        set(h.pp_density_profile_setting_popup_colormap, 'value', 1);
        
        cla(h.pp_density_profile_axes,'reset');
        
        set(h.pp_density_signal_popup_var, 'value', 1);
        set(h.pp_density_signal_x_cb, 'value', 0);
        set(h.pp_density_signal_x_edit, 'string','','userdata','','enable','off');
        set(h.pp_density_signal_y_cb, 'value', 0);
        set(h.pp_density_signal_y_edit, 'string','','userdata','','enable','off');
        
        set(h.pp_density_signal_setting_cb_view, 'value', 0);
        set(h.pp_density_signal_setting_edit_view, 'string','','userdata','','enable','off');
        set(h.pp_density_signal_setting_cb_clim, 'value', 0);
        set(h.pp_density_signal_setting_edit_clim, 'string','','userdata','','enable','off');
        set(h.pp_density_signal_setting_cb_spatlim, 'value', 0);
        set(h.pp_density_signal_setting_edit_spatlim , 'string','','userdata','','enable','off');
        set(h.pp_density_signal_setting_cb_tlim, 'value', 0);
        set(h.pp_density_signal_setting_edit_tlim , 'string','','userdata','','enable','off');
        
        set(h.pp_density_signal_setting_cb_coarse, 'value', 0);
        set(h.pp_density_signal_setting_edit_coarse  , 'string','','userdata','','enable','off');
        set(h.pp_density_signal_setting_cb_savefig, 'value', 0);
        set(h.pp_density_signal_setting_popup_savefig, 'value', 1,'enable', 'off');
        set(h.pp_density_signal_setting_popup_colormap, 'value', 1);
        
        cla(h.pp_density_signal_axes,'reset');
        
%         set(h.pp_density_statistic_hs_cb, 'value', 0);
%         set(h.pp_density_statistic_mtc_cb, 'value', 0);
%         set(h.pp_density_statistic_mtt_cb, 'value', 0);
%         set(h.pp_density_statistic_mwl_cb, 'value', 0);
        
        set(h.pp_density_statistic_setting_cb_view, 'value', 0);
        set(h.pp_density_statistic_setting_edit_view,'string', '','userdata', '','enable','off');
        set(h.pp_density_statistic_setting_cb_clim, 'value', 0);
        set(h.pp_density_statistic_setting_edit_clim,'string', '','userdata', '','enable','off');
        set(h.pp_density_statistic_setting_cb_xlim,'value',0);
        set(h.pp_density_statistic_setting_edit_xlim,'string', '','userdata', '','enable','off');
        set(h.pp_density_statistic_setting_cb_ylim,'value',0);
        set(h.pp_density_statistic_setting_edit_ylim,'string','','userdata','','enable','off');
        set(h.pp_density_statistic_setting_cb_tlim, 'value', 0);
        set(h.pp_density_statistic_setting_edit_tlim , 'string','','userdata','','enable','off');
        set(h.pp_density_statistic_setting_cb_coarse, 'value',0);
        set(h.pp_density_statistic_setting_edit_coarse,'string', '','userdata', '','enable','off');
        set(h.pp_density_statistic_setting_cb_savefig, 'value',0);
        set(h.pp_density_statistic_setting_popup_savefig, 'value', 1, 'enable', 'off');
        set(h.pp_density_statistic_setting_popup_colormap, 'value', 1);
        
        cla(h.pp_density_statistic_axes);
        
        set(h.pp_line_profile_popup_var, 'value', 1);
        set(h.pp_line_profile_time_edit, 'string', '', 'userdata', '');
        
        set(h.pp_line_profile_x_cb, 'value', 0);
        set(h.pp_line_profile_x_edit,'string', '','userdata', '','enable','off');
        set(h.pp_line_profile_y_cb, 'value', 0);
        set(h.pp_line_profile_y_edit, 'string', '','userdata', '','enable','off');
        
        set(h.pp_line_profile_setting_cb_zlim, 'value', 0);
        set(h.pp_line_profile_setting_edit_zlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_profile_setting_cb_spatlim, 'value', 0);
        set(h.pp_line_profile_setting_edit_spatlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_profile_setting_cb_coarse, 'value', 0);
        set(h.pp_line_profile_setting_edit_coarse, 'string', '','userdata', '','enable','off');
        set(h.pp_line_profile_setting_cb_savefig, 'value', 0);
        set(h.pp_line_profile_setting_popup_savefig, 'value', 1, 'enable', 'off');
        
        cla(h.pp_line_profile_axes);
        
        set(h.pp_line_buoy_popup_var, 'value', 1);
        set(h.pp_line_buoy_x_edit,'string', '','userdata', '');
        set(h.pp_line_buoy_y_edit,'string', '','userdata', '');
        set(h.pp_line_buoy_cb_signal, 'value', 1, 'enable', 'on');
        set(h.pp_line_buoy_cb_spectrum, 'value', 0, 'enable', 'off');
        set(h.pp_line_buoy_popup_spectrum_var, 'value', 1, 'enable', 'off');
        set(h.pp_line_buoy_setting_cb_zlim, 'value', 0);
        set(h.pp_line_buoy_setting_edit_zlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_buoy_setting_cb_horznlim, 'value', 0);
        set(h.pp_line_buoy_setting_edit_horznlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_buoy_setting_cb_coarse, 'value', 0);
        set(h.pp_line_buoy_setting_edit_coarse, 'string', '','userdata', '','enable','off');
        set(h.pp_line_buoy_setting_cb_savefig, 'value', 0);
        set(h.pp_line_buoy_setting_popup_savefig, 'value', 1,'enable','off');
        
        cla(h.pp_line_buoy_axes);
        
        set(h.pp_line_statistic_popup_var, 'value', 1);
        set(h.pp_line_statistic_x_cb, 'value', 0);
        set(h.pp_line_statistic_x_edit,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_statistic_y_cb, 'value', 0);
        set(h.pp_line_statistic_y_edit,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_statistic_setting_cb_zlim, 'value', 0);
        set(h.pp_line_statistic_setting_edit_zlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_statistic_setting_cb_spatlim, 'value', 0);
        set(h.pp_line_statistic_setting_edit_spatlim, 'string', '','userdata', '','enable','off');
        set(h.pp_line_statistic_setting_cb_tlim, 'value', 0);
        set(h.pp_line_statistic_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_statistic_setting_cb_coarse, 'value', 0);
        set(h.pp_line_statistic_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_statistic_setting_edit_coarse, 'string', '','userdata', '','enable','off');
        set(h.pp_line_statistic_setting_cb_savefig, 'value', 0);
        set(h.pp_line_statistic_setting_popup_savefig, 'value', 1,'enable','off');
        cla(h.pp_line_statistic_axes);
        
        set(h.pp_line_ham_mom_popup_var, 'value', 1);
        set(h.pp_line_ham_mom_setting_cb_xlim, 'value', 0);
        set(h.pp_line_ham_mom_setting_edit_xlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_ham_mom_setting_cb_ylim, 'value', 0);
        set(h.pp_line_ham_mom_setting_edit_ylim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_ham_mom_setting_cb_tlim, 'value', 0);
        set(h.pp_line_ham_mom_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_ham_mom_setting_cb_zlim, 'value', 0);
        set(h.pp_line_ham_mom_setting_edit_zlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_ham_mom_setting_cb_coarse, 'value', 0);
        set(h.pp_line_ham_mom_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_ham_mom_setting_cb_savefig, 'value', 0);
        set(h.pp_line_ham_mom_setting_popup_savefig, 'value', 1,'enable','off');
        
        cla(h.pp_line_ham_mom_axes);
        
        set(h.pp_line_breaking_popup_var, 'value', 1);
        set(h.pp_line_breaking_setting_cb_xlim, 'value', 0);
        set(h.pp_line_breaking_setting_edit_xlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_breaking_setting_cb_ylim, 'value', 0);
        set(h.pp_line_breaking_setting_edit_ylim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_breaking_setting_cb_tlim, 'value', 0);
        set(h.pp_line_breaking_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_breaking_setting_cb_zlim, 'value', 0);
        set(h.pp_line_breaking_setting_edit_zlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_breaking_setting_cb_coarse, 'value', 0);
        set(h.pp_line_breaking_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_breaking_setting_cb_savefig, 'value', 0);
        set(h.pp_line_breaking_setting_popup_savefig, 'value', 1,'enable','off');
        
        cla(h.pp_line_breaking_axes);
        
        
        set(h.pp_line_MTAA_popup_var, 'value', 1);
        set(h.pp_line_MTAA_setting_cb_xlim, 'value', 0);
        set(h.pp_line_MTAA_setting_edit_xlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_MTAA_setting_cb_ylim, 'value', 0);
        set(h.pp_line_MTAA_setting_edit_ylim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_MTAA_setting_cb_tlim, 'value', 0);
        set(h.pp_line_MTAA_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_MTAA_setting_cb_zlim, 'value', 0);
        set(h.pp_line_MTAA_setting_edit_zlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_MTAA_setting_cb_coarse, 'value', 0);
        set(h.pp_line_MTAA_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_MTAA_setting_cb_savefig, 'value', 0);
        set(h.pp_line_MTAA_setting_popup_savefig, 'value', 1,'enable','off');
        
        cla(h.pp_line_MTAA_axes);
        
        set(h.pp_line_Extreme_popup_plotopt, 'value', 1);
        set(h.pp_line_Extreme_setting_cb_xlim, 'value', 0);
        set(h.pp_line_Extreme_setting_edit_xlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_Extreme_setting_cb_ylim, 'value', 0);
        set(h.pp_line_Extreme_setting_edit_ylim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_Extreme_setting_cb_tlim, 'value', 0);
        set(h.pp_line_Extreme_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_Extreme_setting_cb_zlim, 'value', 0);
        set(h.pp_line_Extreme_setting_edit_zlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_Extreme_setting_cb_coarse, 'value', 0);
        set(h.pp_line_Extreme_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_line_Extreme_setting_cb_savefig, 'value', 0);
        set(h.pp_line_Extreme_setting_popup_savefig, 'value', 1,'enable','off');
        
        cla(h.pp_line_Extreme_axes);
        
        
        set(h.pp_anim_density_popup_var, 'value', 1);
        set(h.pp_anim_density_cb_level_eta, 'value', 0);
        set(h.pp_anim_density_cb_level_phi, 'value', 0);
        set(h.pp_anim_density_cb_level_quiver, 'value', 0, 'enable', 'off');
        set(h.pp_anim_density_setting_cb_view, 'value', 0);
        set(h.pp_anim_density_setting_edit_view,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_clim, 'value', 0);
        set(h.pp_anim_density_setting_edit_clim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_xlim, 'value', 0);
        set(h.pp_anim_density_setting_edit_xlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_ylim, 'value', 0);
        set(h.pp_anim_density_setting_edit_ylim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_tlim , 'value', 0);
        set(h.pp_anim_density_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_coarse, 'value', 0);
        set(h.pp_anim_density_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_cb_level, 'value', 0);
        set(h.pp_anim_density_setting_edit_level  , 'string','','userdata','','enable','off');
        
        set(h.pp_anim_density_setting_cb_saveanim, 'value', 0);
        set(h.pp_anim_density_setting_edit_gifset,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_density_setting_popup_colormap, 'value', 1,'enable','on');
        
        cla(h.pp_anim_density_axes);
        
        set(h.pp_anim_line_popup_var, 'value', 1);
        set(h.pp_anim_line_x_cb, 'value', 0);
        set(h.pp_anim_line_x_edit,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_y_cb, 'value', 0);
        set(h.pp_anim_line_y_edit,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_setting_cb_zlim, 'value', 0);
        set(h.pp_anim_line_setting_edit_zlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_setting_cb_spatlim, 'value', 0);
        set(h.pp_anim_line_setting_edit_spatlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_setting_cb_tlim, 'value', 0);
        set(h.pp_anim_line_setting_edit_tlim,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_setting_cb_coarse, 'value', 0);
        set(h.pp_anim_line_setting_edit_coarse,'string', '','userdata', '', 'enable', 'off');
        set(h.pp_anim_line_setting_cb_saveanim, 'value', 0);
        set(h.pp_anim_line_setting_edit_gifset,'string', '','userdata', '', 'enable', 'off');
        
        cla(h.pp_anim_line_axes);
        
        set(h.pp_validation_buoy_popup_var,'value',1);
        set(h.pp_validation_buoy_measdata_button,'userdata',[]);
        set(h.pp_validation_buoy_x_edit,'string','','userdata',[]);
        set(h.pp_validation_buoy_y_edit,'string','','userdata',[]);
        set(h.pp_validation_buoy_time_edit,'string','','userdata',[]);
        set(h.pp_validation_buoy_cb_timeshift,'value',0);
        set(h.pp_validation_buoy_cb_timeshift_def,'value',0,'enable','off')
        set(h.pp_validation_buoy_cb_timeshift_dt,'value',0,'enable','off')
        set(h.pp_validation_buoy_edit_timeshift_dt,'string','','userdata',[],'enable','off')
        set(h.pp_validation_buoy_cb_signal,'value',1,'enable','on')
        set(h.pp_validation_buoy_cb_spectrum,'value',0,'enable','off');
        set(h.pp_validation_buoy_popup_spectrum_var,'value',1,'enable','off')
        set(h.pp_validation_buoy_number_popup,'String',{'1'},'value',1)
        set(h.pp_validation_buoy_X1disp_edit,'string','')
        set(h.pp_validation_buoy_Y1disp_edit,'string','')
        set(h.pp_validation_buoy_setting_cb_zlim, 'value', 0);
        set(h.pp_validation_buoy_setting_edit_zlim, 'string', '','userdata', '','enable','off');
        set(h.pp_validation_buoy_setting_cb_horznlim, 'value', 0);
        set(h.pp_validation_buoy_setting_edit_horznlim, 'string', '','userdata', '','enable','off');
        set(h.pp_validation_buoy_setting_cb_coarse, 'value', 0);
        set(h.pp_validation_buoy_setting_edit_coarse, 'string', '','userdata', '','enable','off');
        set(h.pp_validation_buoy_setting_cb_savefig, 'value', 0);
        set(h.pp_validation_buoy_setting_popup_savefig, 'value', 1,'enable','off');

        cla(h.pp_validation_buoy_axes);
        
        
        
    end

    function reset_handles_gui(h)
        set(h.project_popup_projdir,'string','--')
        set(h.project_edit_name,'string','');
        set(h.project_edit_note,'string','');
        set(h.model_popup_dynamic,'value',1);
        set(h.model_popup_dispersion,'value',1);
        set(h.model_popup_breaking,'value',2);
        set(h.model_break_edit_initiation,'string','0.8','userdata',0.8);
        set(h.model_break_edit_termination,'string','0.2','userdata',0.2);
        set(h.model_break_edit_Tstar,'string','0.5','userdata',0.5);
        callback_breaking(h.model_popup_breaking,[],h);
        set(h.model_current_cb,'value',0);
        set(h.model_current_ux_edit,'userdata',[],'string','');
        set(h.model_current_uy_edit,'userdata',[],'string','');
        callback_current_cb(h.model_current_cb,[],h);
        
        set(h.spatial_edit_xmin,'userdata',[],'string','');
        set(h.spatial_edit_xmax,'userdata',[],'string','');
        set(h.spatial_edit_ymin,'userdata',[],'string','');
        set(h.spatial_edit_ymax,'userdata',[],'string','');
        set(h.spatial_edit_dx,'userdata',[],'string','');
        set(h.spatial_edit_dy,'userdata',[],'string','');
        set(h.cutfrac_k,'userdata',2,'string',2);
        set(h.wall_popup,'value',1);
        callback_wall_popup(h.wall_popup,[],h);
        
        set(h.damping_popup,'value',1);
        set(h.damping_table,'data',[],'userdata',[],'visible','off');
        callback_damping_popup(h.damping_popup,[],h);
        
        set(h.bathymetry_popup_type,'value',1);
        set(h.bathymetry_edit_depth,'userdata',[],'string','');
        set(h.bathymetry_edit_mindepth,'userdata',[],'string','');
        set(h.bathymetry_edit_maxdepth,'userdata',[],'string','');
        set(h.bathymetry_edit_slope,'userdata',[],'string','');
        set(h.bathymetry_edit_startslope,'userdata',[],'string','');
        set(h.bathymetry_edit_maxdepthshore,'userdata',[],'string','');
        set(h.bathymetry_edit_slopeshore,'userdata',[],'string','');
        set(h.bathymetry_edit_shoreposition,'userdata',[],'string','');
        set(h.bathymetry_button_load ,'userdata',[]);
        set(h.bathymetry_popup_interpdepth,'value',1);
        callback_bathy_popup(h.bathymetry_popup_type,[],h);%%it is called to asjust again all handles
        set(h.bathymetry_checkbox_friction,'value',1);
        set(h.bathymetry_table_friction,'data',[],'visible','off');
        
        
        set(h.waveinput_ivp_popup_type,'value',1);
        set(h.waveinput_ivp_edit_A,'userdata',[],'string','');
        set(h.waveinput_ivp_edit_centre_position,'userdata',[],'string','');
        set(h.waveinput_ivp_edit_sd,'userdata',[],'string','');
        set(h.waveinput_ivp_button_load,'userdata',[]);
        set(h.waveinput_ivp_edit_filename,'string','');
        callback_ivp_popup(h.waveinput_ivp_popup_type,[],h);
        
        set(h.waveinput_influx_popup,'value',1);
        set(h.waveinput_influx_proptable,'data',[],'userdata',[],'visible','off');
        set(h.waveinput_influx_methodtable,'data',[],'visible','off');
        set(h.waveinput_influx_checkbox_ramp,'value',1);
        set(h.waveinput_influx_edit_ramp,'userdata',2,'string','2');
        set(h.waveinput_influx_checkbox_nonlinadj,'value',0,'enable','off');
        set(h.waveinput_influx_edit_nonlinadj,'userdata',[],'string','','enable','off');
        callback_influx_popup(h.waveinput_influx_popup,[],h)
        
        set(h.waveinput_bdy_assim_popup,'value',1);
        set(h.waveinput_bdy_assim_shape_popup,'value',1);
        set(h.waveinput_bdy_assim_shape_popup,'userdata',[]);
        set(h.waveinput_bdy_assim_edit_R1,'userdata',[]);
        set(h.waveinput_bdy_assim_edit_xc,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_edit_yc,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_edit_smooth,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_button_load,'userdata',[]);
        if ID_model_phiform==1
        set(h.waveinput_bdy_assim_button_load_phi,'userdata',[]);
        set(h.waveinput_bdy_assim_cb_phi ,'value',0);
        else
        set(h.waveinput_bdy_assim_button_load_v,'userdata',[]);    
        set(h.waveinput_bdy_assim_button_load_u,'userdata',[]);
        set(h.waveinput_bdy_assim_cb_vel ,'value',0);    
        end
        set(h.waveinput_bdy_assim_propdir_popup,'value',1);
        set(h.waveinput_bdy_assim_cb_nonlinAdj,'value',0);
        set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'userdata',[],'string','');
        
        set(h.waveinput_bdy_assim_time_edit_interval_init,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_time_edit_interval_end,'userdata',[],'string','');
        set(h.waveinput_bdy_assim_time_edit_step,'userdata',[],'string','');
       if ID_model_phiform==1
        callback_bdyassim_nonlinadj(h.waveinput_bdy_assim_cb_phi,[],h);
        callback_bdyaasim_cb_phi(h.waveinput_bdy_assim_cb_phi,[],h);
       else
        callback_bdyassim_nonlinadj(h.waveinput_bdy_assim_cb_vel,[],h);
        callback_bdyaasim_cb_vel(h.waveinput_bdy_assim_cb_vel,[],h);  
       end
        callback_bdyassim_popup(h.waveinput_bdy_assim_popup,[],h);
        
        
        
        set(h.time_edit_interval_init,'userdata',[],'string','');
        set(h.time_edit_interval_end,'userdata',[],'string','');
        set(h.time_edit_step,'userdata',[],'string','');
        set(h.options_outputvar_popup,'value',1);
        set(h.options_checkbox_intflow,'value',0,'enable','on');
        set(h.options_edit_interval_init,'userdata',[],'string','','enable','off');
        set(h.options_edit_interval_end,'userdata',[],'string','','enable','off');
        set(h.options_edit_step,'userdata',[],'string','','enable','off');
        set(h.options_checkbox_default,'value',1);
        set(h.options_edit_partition,'userdata',[],'string','','enable','off');
        callback_ode_partition(h.options_checkbox_default,[],h)
        set(h.options_checkbox_combine,'value',1);
        set(h.options_edit_ode_tol,'userdata',0.001,'string','0.001')
        set(h.options_odesolv_popup,'value', 1,'visible','on');
        cla(h.preview.dispersion_axes1);
        cla(h.preview.dispersion_axes2);
        cla(findall(h.panel.prev_spatial,'type','axes'));
        delete(findall(h.panel.prev_wave,'type','axes'));
        set(h.preview.log,'string','')
    end

    function input_handles_gui(h)
        
        GUIinput=h.GUIinput;
        set(h.project_edit_name,'string',GUIinput.proj.name);
        set(h.pp_project_edit_name,'string',GUIinput.proj.name);
        set(h.pp_project_edit_name,'string',GUIinput.proj.name);
        set(h.IF_project_edit_name,'string',GUIinput.proj.name)
        set(h.project_edit_note,'string',GUIinput.proj.note);
        if h.flagOpenProj==1
            set(h.project_popup_projdir,'string',GUIinput.proj.projdir);
            set(h.pp_project_popup_projdir,'string',GUIinput.proj.projdir);
            set(h.IF_project_popup_projdir,'string',GUIinput.proj.projdir);
        else
            set(h.project_popup_projdir,'string','--');
            set(h.pp_project_popup_projdir,'string','--');
            set(h.IF_project_popup_projdir,'string','--');
        end
        
        set(h.model_popup_dynamic,'value',GUIinput.modeldynVal);
        set(h.model_popup_dispersion,'value',GUIinput.modeldispVal);
        set(h.model_popup_breaking,'value',GUIinput.modelbreak.checkVal);
        set(h.model_break_edit_initiation,'string',...
            num2str(GUIinput.modelbreak.initiation),'userdata',...
            GUIinput.modelbreak.initiation);
        set(h.model_break_edit_termination,'string',...
            num2str(GUIinput.modelbreak.termination),'userdata',...
            GUIinput.modelbreak.termination);
        set(h.model_break_edit_Tstar,'string',num2str(GUIinput.modelbreak.Tstar),...
            'userdata',GUIinput.modelbreak.Tstar);
        callback_breaking(h.model_popup_breaking,[],h);
        if isfield(GUIinput,'modelcurrent_cb')
            set(h.model_current_cb,'value',GUIinput.modelcurrent.check);
            set(h.model_current_ux_edit,'userdata',GUIinput.modelcurrent.ux,'string',...
                num2str(GUIinput.modelcurrent.ux));
            set(h.model_current_uy_edit,'userdata',GUIinput.modelcurrent.uy,'string',...
                num2str(GUIinput.modelcurrent.uy));
            callback_current_cb(h.model_current_cb,[],h);
        end
        
        set(h.spatial_edit_xmin,'userdata',GUIinput.spatialinterv.xmin,...
            'string',num2str(GUIinput.spatialinterv.xmin));
        set(h.spatial_edit_xmax,'userdata',GUIinput.spatialinterv.xmax,...
            'string',num2str(GUIinput.spatialinterv.xmax));
        set(h.spatial_edit_ymin,'userdata',GUIinput.spatialinterv.ymin,...
            'string',num2str(GUIinput.spatialinterv.ymin));
        set(h.spatial_edit_ymax,'userdata',GUIinput.spatialinterv.ymax,...
            'string',num2str(GUIinput.spatialinterv.ymax));
        
        if ~isfield(GUIinput,'spatialgrid')
            Dx=(GUIinput.spatialinterv.xmax-GUIinput.spatialinterv.xmin)...
                ./(2^GUIinput.spatialdxdy.px);
            Dy=(GUIinput.spatialinterv.ymax-GUIinput.spatialinterv.ymin)...
                ./(2^GUIinput.spatialdxdy.py);
            set(h.spatial_edit_dx,'userdata',Dx,...
                'string',num2str(Dx));
            set(h.spatial_edit_dy,'userdata',Dy,...
                'string',num2str(Dy));
        else
            set(h.spatial_edit_dx,'userdata',GUIinput.spatialgrid.dx,...
                'string',num2str(GUIinput.spatialgrid.dx));
            set(h.spatial_edit_dy,'userdata',GUIinput.spatialgrid.dy,...
                'string',num2str(GUIinput.spatialgrid.dy));
        end
        
        set(h.cutfrac_k,'userdata',GUIinput.fourier_cutfrac.k,...
            'string',GUIinput.fourier_cutfrac.k);
        set(h.wall_popup,'value',GUIinput.spatialwall.checkVal);
        if strcmpi(h.moduleRestrict,'Flat_NoWall')
            if GUIinput.spatialwall.checkVal==1
                set(h.wall_popup,'value',2);
                set(h.monitorbox,'String','>> The wall is not available for this licence. ','foregroundcolor','k');
            end
        end
      
        
        tableInputdat=GUIinput.spatialwall.data;
        for ii=1:length(tableInputdat(:,1))
            for jj=1:length(tableInputdat(1,:))
                if isempty(tableInputdat{ii,jj})
                   tableInputdat(ii,jj)={0}; 
                end
            end
        end
        
        if length(tableInputdat(1,:))==11
            temp(:,1:8)=tableInputdat(:,1:8);
            for ii=1:length(tableInputdat(:,1))
                        temp(ii,9)={'Energy truncation'};
            end
            temp(:,10:12)=tableInputdat(:,9:11);
            tableInputdat=temp;
        end
        
        for ii=1:length(tableInputdat(:,1))
            if strcmpi(tableInputdat(ii,9),'Localization')
            tableInputdat(ii,9)={'Energy truncation'};
            end
        end
        
        set(h.wall_table,'userdata',GUIinput.spatialwall.userdata,...
            'data',tableInputdat);
        
        callback_wall_popup(h.wall_popup,[],h);
        
        if GUIinput.spatialwall.checkVal == 1
            set(h.wall_table,'visible','on')
        end
        
        set(h.damping_popup,'value',GUIinput.spatialdamp.checkVal);
        NdF=length(GUIinput.spatialdamp.data);
        datTab=GUIinput.spatialdamp.data;
       
        if NdF==5
        datTab=    GUIinput.spatialdamp.data;
        datTab(1,6)={'-'};
        if length(datTab(:,1))>1
           datTab(2:end,6)={0}; 
        end
        end
        set(h.damping_table,'data',datTab,...
            'userdata',GUIinput.spatialdamp.userdata);
        callback_damping_popup(h.damping_popup,[],h);
        if GUIinput.spatialdamp.checkVal==1
            set(h.damping_table,'visible','on')
        end
        
        set(h.bathymetry_popup_type,'value',GUIinput.spatialbathy.typeVal);
        if strcmpi(h.moduleRestrict,'Flat_NoWall')
            if GUIinput.spatialbathy.typeVal>1
                set(h.bathymetry_popup_type,'value',1);
                set(h.monitorbox,'String','>> Only Flat bottom is available for this licence. ','foregroundcolor','k');
            end
        end
        
        set(h.bathymetry_edit_depth,'userdata',GUIinput.spatialbathy.depth,...
            'string',num2str(GUIinput.spatialbathy.depth));
        set(h.bathymetry_edit_mindepth,'userdata',GUIinput.spatialbathy.xymindepth,...
            'string',num2str(GUIinput.spatialbathy.xymindepth));
        set(h.bathymetry_edit_maxdepth,'userdata',GUIinput.spatialbathy.xymaxdepth,...
            'string',num2str(GUIinput.spatialbathy.xymaxdepth));
        set(h.bathymetry_edit_slope,'userdata',GUIinput.spatialbathy.slope,...
            'string',num2str(GUIinput.spatialbathy.slope));
        set(h.bathymetry_edit_startslope,'userdata',GUIinput.spatialbathy.startslope,...
            'string',num2str(GUIinput.spatialbathy.startslope));
        set(h.bathymetry_edit_maxdepthshore,'userdata',GUIinput.spatialbathy.maxdepthshore,...
            'string',num2str(GUIinput.spatialbathy.maxdepthshore));
        if isfield(GUIinput.spatialbathy,'mindepthshore')
          set(h.bathymetry_edit_mindepthshore,'userdata',GUIinput.spatialbathy.mindepthshore,...
            'string',num2str(GUIinput.spatialbathy.mindepthshore));
        end
          
        set(h.bathymetry_edit_slopeshore,'userdata',GUIinput.spatialbathy.slopeshore,...
            'string',num2str(GUIinput.spatialbathy.slopeshore));
        set(h.bathymetry_edit_shoreposition,'userdata',GUIinput.spatialbathy.shoreposition,...
            'string',num2str(GUIinput.spatialbathy.shoreposition));
        set(h.bathymetry_button_load ,'userdata',GUIinput.spatialbathy.userdata);
        if ~isempty(GUIinput.spatialbathy.userdata)
        set(h.bathymetry_edit_filename,'string','data loaded');
        end
        set(h.bathymetry_popup_interpdepth,'value',GUIinput.spatialbathy.interp-1);
       % nunu
        if isfield(GUIinput.spatialbathy,'xymiddepth')
            set(h.bathymetry_edit_middepth,'userdata',GUIinput.spatialbathy.xymiddepth,...
                'string',num2str(GUIinput.spatialbathy.xymiddepth));
        else
            set(h.bathymetry_edit_middepth,'userdata',0.5*(GUIinput.spatialbathy.xymindepth+GUIinput.spatialbathy.xymaxdepth),...
                'string',num2str(0.5*(GUIinput.spatialbathy.xymindepth+GUIinput.spatialbathy.xymaxdepth)));
        end
        callback_bathy_popup_interpdepth(h.bathymetry_popup_interpdepth,[],h);
        set(h.bathymetry_checkbox_friction,'value',GUIinput.spatialfriction.check);
        callback_bathy_popup(h.bathymetry_popup_type,[],h);%%it is called to adjust again all handles
        
        if GUIinput.spatialfriction.check==1
            set(h.bathymetry_table_friction,'visible','on');
            set(h.bathymetry_button_addrow,'visible','on');
            set(h.bathymetry_button_deleterow,'visible','on');
            
            if length(GUIinput.spatialfriction.data(1,:))==5
                set(h.bathymetry_table_friction,'data',[]);
            else
                set(h.bathymetry_table_friction,'data',GUIinput.spatialfriction.data);
                set(h.bathymetry_table_friction,'userdata',GUIinput.spatialfriction.userdata);
            end
        end
        set(h.waveinput_ivp_popup_type,'value',GUIinput.wave_ivp.typeVal);
        set(h.waveinput_ivp_edit_A,'userdata',GUIinput.wave_ivp.A,'string',...
            num2str(GUIinput.wave_ivp.A));
        set(h.waveinput_ivp_edit_centre_position,'userdata',GUIinput.wave_ivp.centerposition,...
            'string',num2str(GUIinput.wave_ivp.centerposition));
        set(h.waveinput_ivp_edit_sd,'userdata',GUIinput.wave_ivp.stdev,...
            'string',num2str(GUIinput.wave_ivp.stdev));
        set(h.waveinput_ivp_button_load,'userdata',GUIinput.wave_ivp.userdata);
        if isfield(GUIinput.wave_ivp.userdata,'name')
            set(h.waveinput_ivp_edit_filename,'string',GUIinput.wave_ivp.userdata.name);
        end
        callback_ivp_popup(h.waveinput_ivp_popup_type,[],h);
        
        set(h.waveinput_influx_popup,'value',GUIinput.wave_influx.typeVal);
        
        
        propertiesdata=GUIinput.wave_influx.propertiesdata;
        
        Ninf=length(GUIinput.wave_influx.propertiesdata(:,1));
        for ii=1:Ninf
            if strcmpi(propertiesdata(ii,1),'User-defined')
                propertiesdata(ii,1)={'User-defined (signal)'};
            end
            if strcmp(propertiesdata(ii,1),'Jonswap')
                    propertiesdata(ii,1)={'JONSWAP'};
            end
        end
        
        if length(propertiesdata(1,:))==8
            temp=propertiesdata;
            propertiesdata=cell(Ninf,7);
            for ii=1:Ninf
                if strcmp(temp(ii,1),'Harmonic') &&  strcmp(temp(ii,7),'Uniform')
                propertiesdata(ii,1:6)=temp(ii,1:6);
                propeiesdata(ii,7)={0};
                else
                propertiesdata(ii,1:6)=temp(ii,1:6);
                propertiesdata(ii,7)=temp(ii,8);
                end
            end
        end
        
        if length(GUIinput.wave_influx.methoddata(1,:))==11
            for ii=1:Ninf
                if strcmpi(propertiesdata(ii,1),'User-defined (signal)')
                    if  strcmp(GUIinput.wave_influx.methoddata(ii,2),'Horizontal')
                        Xinfl=GUIinput.wave_influx.userdata(ii).spatial;
                        wave_influx_userdata(ii).inflX=Xinfl;
                        wave_influx_userdata(ii).inflY=cell2mat(GUIinput.wave_influx.methoddata(ii,3)).*...
                            ones(size(Xinfl));
                    else
                        Yinfl=GUIinput.wave_influx.userdata(ii).spatial;
                        wave_influx_userdata(ii).inflX=cell2mat(GUIinput.wave_influx.methoddata(ii,6)).*...
                            ones(size(Yinfl));
                        wave_influx_userdata(ii).inflY=Yinfl;
                    end
                    wave_influx_userdata(ii).time=GUIinput.wave_influx.userdata(ii).time;
                    wave_influx_userdata(ii).eta=GUIinput.wave_influx.userdata(ii).eta;
                else
                    wave_influx_userdata(ii).inflX=[];  wave_influx_userdata(ii).inflY=[];
                    wave_influx_userdata(ii).time=[];  wave_influx_userdata(ii).eta=[];
                end
            end
             set(h.waveinput_influx_proptable,'data',propertiesdata,...
            'userdata',wave_influx_userdata);
        else
             set(h.waveinput_influx_proptable,'data',propertiesdata,...
            'userdata',GUIinput.wave_influx.userdata);
        end
        
       
        
        if length(GUIinput.wave_influx.methoddata(1,:))==11
            Ninf=length(GUIinput.wave_influx.methoddata(:,1));
            methoddata=cell(Ninf,14);
            methoddata(1:Ninf,1)=GUIinput.wave_influx.methoddata(:,1);
            methoddata(:,12:14)=GUIinput.wave_influx.methoddata(:,9:11);
            for ii=1:Ninf
                if strcmp(GUIinput.wave_influx.methoddata(ii,1),'Line')
                 methoddata(ii,1)={'Point'};   
                end
                 methoddata(ii,2)={'Straight'};
                 methoddata(ii,7:10)={'-'};
                 if  strcmp(GUIinput.wave_influx.methoddata(ii,2),'Horizontal')
                 methoddata(ii,3)=GUIinput.wave_influx.methoddata(ii,4);
                 methoddata(ii,4)=GUIinput.wave_influx.methoddata(ii,3);
                 methoddata(ii,5)=GUIinput.wave_influx.methoddata(ii,5);
                 methoddata(ii,6)=GUIinput.wave_influx.methoddata(ii,3);
                 else
                 methoddata(ii,3)=GUIinput.wave_influx.methoddata(ii,6);
                 methoddata(ii,4)=GUIinput.wave_influx.methoddata(ii,7);
                 methoddata(ii,5)=GUIinput.wave_influx.methoddata(ii,6);
                 methoddata(ii,6)=GUIinput.wave_influx.methoddata(ii,8);    
                 end
            end
        else
            methoddata=GUIinput.wave_influx.methoddata;
        end
        
      
        
        set(h.waveinput_influx_methodtable,'data',methoddata);
        set(h.waveinput_influx_checkbox_ramp,'value',GUIinput.wave_influx.rampcheck);
        set(h.waveinput_influx_edit_ramp,'userdata',GUIinput.wave_influx.rampfactor,...
            'string',num2str(GUIinput.wave_influx.rampfactor));
        if isfield(GUIinput.wave_influx,'ramplinecheck')
            set(h.waveinput_influxline_checkbox_ramp,'value',GUIinput.wave_influx.ramplinecheck);
            set(h.waveinput_influxline_edit_ramp,'userdata',GUIinput.wave_influx.ramplinefactor,...
                'string',num2str(GUIinput.wave_influx.ramplinefactor));
        end
        set(h.waveinput_influx_checkbox_nonlinadj,'value',GUIinput.wave_influx.nonlinadjcheck);
        set(h.waveinput_influx_edit_nonlinadj,'userdata',GUIinput.wave_influx.nonlinadjfactor,...
            'string',num2str(GUIinput.wave_influx.nonlinadjfactor));
        callback_influx_popup(h.waveinput_influx_popup,[],h)
        
        if isfield(GUIinput,'wave_bdy_assim')
            set(h.waveinput_bdy_assim_popup,'value',GUIinput.wave_bdy_assim.checkVal);
            callback_bdyassim_popup(h.waveinput_bdy_assim_popup,[],h);
            if isfield(GUIinput.wave_bdy_assim,'shapeCheck');
                set(h.waveinput_bdy_assim_shape_popup,'value',GUIinput.wave_bdy_assim.shapeCheck,...
                    'userdata',GUIinput.wave_bdy_assim.shapeUserdata);
                callback_bdyassim_shape(h.waveinput_bdy_assim_shape_popup,[],h);
                set(h.waveinput_bdy_assim_edit_R1,'userdata',GUIinput.wave_bdy_assim.R1,...
                    'string',num2str(GUIinput.wave_bdy_assim.R1))
                set(h.waveinput_bdy_assim_edit_xc,'userdata',GUIinput.wave_bdy_assim.xc,...
                    'string',num2str(GUIinput.wave_bdy_assim.xc))
                set(h.waveinput_bdy_assim_edit_yc,'userdata',GUIinput.wave_bdy_assim.yc,...
                    'string',num2str(GUIinput.wave_bdy_assim.yc))
                set(h.waveinput_bdy_assim_edit_smooth,'userdata',GUIinput.wave_bdy_assim.smoothfact,...
                    'string',num2str(GUIinput.wave_bdy_assim.smoothfact));
                set(h.waveinput_bdy_assim_button_load,'userdata',GUIinput.wave_bdy_assim.assimdata);
                if ~isempty(GUIinput.wave_bdy_assim.assimdata)
                    set(h.waveinput_bdy_assim_edit_data,'string','data loaded');
                end
                set(h.waveinput_bdy_assim_propdir_popup,'value',GUIinput.wave_bdy_assim.propdir);
                
                if ID_model_phiform==1
                if isfield(GUIinput.wave_bdy_assim,'cb_phi')
                    set(h.waveinput_bdy_assim_cb_phi,'value',GUIinput.wave_bdy_assim.cb_phi);
                    set(h.waveinput_bdy_assim_button_load_phi,'userdata',GUIinput.wave_bdy_assim.assimdata_phi);
                    if ~isempty(GUIinput.wave_bdy_assim.assimdata_phi)
                        set(h.waveinput_bdy_assim_edit_data_phi,'string','data loaded');
                    end
                end
                callback_bdyaasim_cb_phi(h.waveinput_bdy_assim_cb_phi,[],h);
                else
                    if isfield(GUIinput.wave_bdy_assim,'cb_vel')
                        set(h.waveinput_bdy_assim_cb_vel,'value',GUIinput.wave_bdy_assim.cb_vel);
                        set(h.waveinput_bdy_assim_button_load_u,'userdata',GUIinput.wave_bdy_assim.assimdata_u);
                        set(h.waveinput_bdy_assim_button_load_v,'userdata',GUIinput.wave_bdy_assim.assimdata_v);
                        
                        if ~isempty(GUIinput.wave_bdy_assim.assimdata_u)
                            set(h.waveinput_bdy_assim_edit_data_u,'string','data loaded');
                        end
                    else
                        set(h.waveinput_bdy_assim_cb_vel,'value',0);
                        set(h.waveinput_bdy_assim_button_load_u,'userdata',[]);
                        set(h.waveinput_bdy_assim_button_load_v,'userdata',[]);
                    end
                 callback_bdyaasim_cb_vel(h.waveinput_bdy_assim_cb_vel,[],h);   
                end
                
                if isfield(GUIinput.wave_bdy_assim,'cb_nonlinAdj')
                    set(h.waveinput_bdy_assim_cb_nonlinAdj,'value',GUIinput.wave_bdy_assim.cb_nonlinAdj);
                    if isfield(GUIinput.wave_bdy_assim,'nonlinAdj_distance')
                        set(h.waveinput_bdy_assim_edit_nonlinAdj_distance,'userdata',...
                            GUIinput.wave_bdy_assim.nonlinAdj_distance,'string',...
                            num2str(GUIinput.wave_bdy_assim.nonlinAdj_distance));
                        set(h.waveinput_bdy_assim_edit_nonlinAdj_smooth,'userdata',...
                            GUIinput.wave_bdy_assim.nonlinAdj_smooth,'string',...
                            num2str(GUIinput.wave_bdy_assim.nonlinAdj_smooth));
                    end
                end
                
                if ID_model_phiform==1
                callback_bdyassim_nonlinadj(h.waveinput_bdy_assim_cb_phi,[],h)
                else
                callback_bdyassim_nonlinadj(h.waveinput_bdy_assim_cb_vel,[],h)    
                end
                
                set(h.waveinput_bdy_assim_time_edit_interval_init,'userdata',GUIinput.wave_bdy_assim.tinit,...
                    'string',num2str(GUIinput.wave_bdy_assim.tinit));
                set(h.waveinput_bdy_assim_time_edit_interval_end,'userdata',GUIinput.wave_bdy_assim.tend,...
                    'string',num2str(GUIinput.wave_bdy_assim.tend));
                set(h.waveinput_bdy_assim_time_edit_step,'userdata',GUIinput.wave_bdy_assim.dt,...
                    'string',num2str(GUIinput.wave_bdy_assim.dt));
            end
        else
            set(h.waveinput_bdy_assim_popup,'value',1);
            callback_bdyassim_popup(h.waveinput_bdy_assim_popup,[],h);
        end
        
        set(h.time_edit_interval_init,'userdata',GUIinput.time.t_start,...
            'string',num2str(GUIinput.time.t_start));
        set(h.time_edit_interval_end,'userdata',GUIinput.time.t_end,...
            'string',num2str(GUIinput.time.t_end));
        set(h.time_edit_step,'userdata',GUIinput.time.dt,...
            'string',num2str(GUIinput.time.dt));
        if isfield(GUIinput,'option_outputvar')
            set(h.options_outputvar_popup,'value',GUIinput.option_outputvar)
        end
        set(h.options_checkbox_intflow,'value',GUIinput.option_intflow.check,...
            'enable','on');
        set(h.options_edit_interval_init,'userdata',GUIinput.option_intflow.tinit,...
            'string',num2str(GUIinput.option_intflow.tinit),'enable','off');
        set(h.options_edit_interval_end,'userdata',GUIinput.option_intflow.tend,...
            'string',num2str(GUIinput.option_intflow.tend),'enable','off');
        set(h.options_edit_step,'userdata',GUIinput.option_intflow.dt_fact,...
            'string',num2str(GUIinput.option_intflow.dt_fact),'enable','off');
        set(h.options_checkbox_default,'value',GUIinput.option_partition.check_def);
        set(h.options_edit_partition,'userdata',GUIinput.option_partition.total,...
            'string',num2str(GUIinput.option_partition.total),'enable','off');
        callback_ode_partition(h.options_checkbox_default,[],h)
        set(h.options_checkbox_combine,'value',GUIinput.option_partition.check_combine);
        
        if isfield(GUIinput,'option_ode_tol')
            set(h.options_edit_ode_tol,'userdata', GUIinput.option_ode_tol,'string',...
                num2str( GUIinput.option_ode_tol));
        else
            set(h.options_edit_ode_tol,'userdata', 0.001,'string',...
                num2str(0.001));
        end
        
        if isfield(GUIinput,'option_ode_sol')
            set(h.options_odesolv_popup,'value', GUIinput.option_ode_sol,...
                'visible','on','enable','on');
        else
            set(h.options_odesolv_popup,'value', 1,'visible','on','enable','on');
        end
        
        if isfield(GUIinput,'option_mc')
            if isfield(GUIinput.option_mc,'numrun_check')
                set(h.options_montecarlo_cb_Nsimul,'value',GUIinput.option_mc.numrun_check);
                set(h.options_montecarlo_edit_Nsimul,'userdata',num2str(GUIinput.option_mc_numrun_edit));
            else
                set(h.options_montecarlo_cb_Nsimul,'value',0);
                set(h.options_montecarlo_edit_Nsimul,'userdata',[],'string','','enable','off');
            end
            callback_mc_N_run(h.options_montecarlo_cb_Nsimul,[],h);
            
            set(h.options_montecarlo_cb_numset,'value',GUIinput.option_mc.numset_check);
            callback_mc_cb_numset(h.options_montecarlo_cb_numset,[],h);
            set(h.options_montecarlo_edit_px,'userdata',GUIinput.option_mc.numset_px);
            set(h.options_montecarlo_edit_py,'userdata',GUIinput.option_mc.numset_py);
            
            set(h.options_montecarlo_cb_waveinput,'value',GUIinput.option_mc.waveinput_check);
            callback_mc_cb_waveinput(h.options_montecarlo_cb_waveinput,[],h);
            set(h.options_montecarlo_pb_waveinput,'userdata',GUIinput.option_mc.waveinput_userdata);
            callback_mc_popup_waveinput(h.options_montecarlo_popup_waveinput,[],h);
            
            
            set(h.options_montecarlo_cb_check,'value',GUIinput.option_mc.check);
            callback_mc_cb_check(h.options_montecarlo_cb_check,[],h);
            
            
            
            
            if GUIinput.option_mc.numset_check==1
                Ndat=length(GUIinput.option_mc.numset_px);
                if Ndat>1
                    for i=1:Ndat-1
                        strPx{i}=[num2str(GUIinput.option_mc.numset_px(i)),';',];
                    end
                    strPx{Ndat}=num2str(GUIinput.option_mc.numset_px(Ndat));
                else
                    strPx{Ndat}=num2str(GUIinput.option_mc.numset_px(1));
                end
                set(h.options_montecarlo_edit_px,'string',[strPx{1:end}]);
                
                Ndat=length(GUIinput.option_mc.numset_py);
                if Ndat>1
                    for i=1:Ndat-1
                        strPy{i}=[num2str(GUIinput.option_mc.numset_py(i)),';',];
                    end
                    strPy{Ndat}=num2str(GUIinput.option_mc.numset_py(Ndat));
                else
                    strPy{Ndat}=num2str(GUIinput.option_mc.numset_py(1));
                end
                set(h.options_montecarlo_edit_py,'string',[strPy{1:end}]);
            end
        end
      guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
  
    end
    
    
    function callback_IF_projectname(hObj,eventdata,h)
      projname=get(hObj,'String');
        projdir=get(h.IF_project_popup_projdir,'string');
        
        if isempty(projdir)
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Specify a working directory');
            uicontrol(h.IF_project_edit_name);
            return;
        end
        
        if exist([projdir,'/',projname],'dir')
            set(h.monitorbox,'foregroundcolor','r','string',...
                '>>Warning: Project exists already, it will be overwritten');
            uicontrol(h.IF_project_edit_name);
            flag_warn_workdir=1;
        else
            flag_warn_workdir=0;
            set(h.monitorbox,'foregroundcolor','k','string','>>');
        end  
    end

    function callback_IF_project_button_projdir(hObj,eventdata,h)
        pathnow=h.pathnow;
        workingdir = uigetdir(pathnow,'browse a directory');
        if workingdir ~= 0
            set(h.IF_project_popup_projdir,'String',workingdir);
        else
            set(h.monitorbox,'String',['>>Warning: No directory is loaded'],'foregroundcolor','k');
            uicontrol(h.IF_project_button_projdir)
            return;
        end
    end


    function callback_IF_loaddata(hObj,eventdata,h)
        projdir=get(h.IF_project_popup_projdir,'string');
        
        if isempty(projdir)||strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify a project directory for post-processing');
            return;
        end
        
        
        
        [file_name,directory]=uigetfile([projdir,'\','*.txt; *.dat; *.mat; *.asc'],'Load Data file');
        
        if directory~=0
            set(h.monitorbox,'foregroundcolor','k','string','>>loading...');
            pause(0.001)
            my_data=load([directory,file_name]);
            my_data.Load=1;
            if isempty(my_data)
                set(h.IF_project_edit_data,'string','');
                set(h.monitorbox,'foregroundcolor','r','string','>>No data loaded');
            elseif ~isfield(my_data,'dom')||~isfield(my_data,'outkinematic')
                set(h.IF_project_edit_data,'string','')
                set(h.monitorbox,'foregroundcolor','r','string','>>Wrong input file');
            else
                reset_handles_postprocIF(h);
                set(h.IF_project_popup_projdir,'string',projdir);
                set(h.IF_project_edit_name,'string',my_data.Proj.savename);
                set(h.IF_project_edit_data,'string',file_name);
                set(h.IF_project_edit_note,'string',my_data.Proj.usernote);
                set(h.monitorbox,'foregroundcolor','k','string','>>Data loaded');
                set(h.IF_project_load_data,'userdata',my_data);
                input_handles_postprocIF(h,my_data);
            end
            
        end
    end

    function reset_handles_postprocIF(h)
        set(h.IF_project_popup_projdir,'string','--')
        set(h.IF_project_edit_note,'string','');
        set(h.IF_project_edit_data,'string','','userdata',[]);
        set(h.IF_project_edit_name,'string','')
        set(h.IF_calc_x_edit_interval_init,'string','','userdata',[]);
        set(h.IF_calc_x_edit_interval_end,'string','','userdata',[]);
        set(h.IF_calc_x_edit_step,'string','','userdata',[]);
        set(h.IF_calc_y_edit_interval_init,'string','','userdata',[]);
        set(h.IF_calc_y_edit_interval_end,'string','','userdata',[]);
        set(h.IF_calc_y_edit_step,'string','','userdata',[]);
        set(h.IF_calc_z_edit_interval_init,'string','','userdata',[],'enable','on');
        set(h.IF_calc_z_edit_interval_end,'string','','userdata',[],'enable','on');
        set(h.IF_calc_z_edit_step,'string','','userdata',[],'enable','on');
        set(h.IF_cal_z_checkbox_equidistant,'value',1);
        set(h.IF_cal_z_checkbox_nonequidistant,'value',0);
        set(h.IF_calc_z_edit_input,'string','','userdata',[],'enable','off')
        set(h.IF_calc_IDdata_edit,'userdata',1,'string',1)
    
    end

    function input_handles_postprocIF(h,inputdat)
     set(h.IF_calc_time_edit_interval_init,'userdata',inputdat.outkinematic.time(1), ...
          'string',num2str(inputdat.outkinematic.time(1)));
     set(h.IF_calc_time_edit_interval_end,'userdata',inputdat.outkinematic.time(end), ...
          'string',num2str(inputdat.outkinematic.time(end)));
     set(h.IF_calc_time_edit_step,'userdata',1,'string',num2str(1));
     Xi=inputdat.dom.X(1)+inputdat.dom.fbl.l;Xf=inputdat.dom.X(end)-inputdat.dom.fbl.r;
     Yi=inputdat.dom.Y(1)+inputdat.dom.fbl.b;Yf=inputdat.dom.Y(end)-inputdat.dom.fbl.t;
     Zi=min(min(inputdat.dom.bathy.profile));Zf=max(max(max(inputdat.outkinematic.eta)));  
     set(h.IF_calc_x_edit_interval_init,'userdata',Xi,'string',num2str(Xi));
     set(h.IF_calc_x_edit_interval_end,'userdata',Xf,'string',num2str(Xf));
     set(h.IF_calc_x_edit_step,'userdata',1,'string','1');
     set(h.IF_calc_y_edit_interval_init,'userdata',Yi,'string',num2str(Yi));
     set(h.IF_calc_y_edit_interval_end,'userdata',Yf,'string',num2str(Yf));
     set(h.IF_calc_y_edit_step,'userdata',1,'string','1')
     set(h.IF_calc_z_edit_interval_init,'userdata',Zi,'string',num2str(Zi));
     set(h.IF_calc_z_edit_interval_end,'userdata',Zf,'string',num2str(Zf));
     set(h.IF_calc_z_edit_step,'userdata',3,'string','3')
     set(h.IF_calc_IDdata_edit,'userdata',1,'string',1)
      
    end

    function callback_IF_edit(hObj,eventdata,h)
    val=str2num(get(hObj,'string')); 
    set(hObj,'userdata',val);
    end

    function callback_IF_cb_Zequidistant(hObj,eventadata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.IF_calc_z_edit_interval_init,'enable','on');
            set(h.IF_calc_z_edit_interval_end,'enable','on');
            set(h.IF_calc_z_edit_step,'enable','on');
            set(h.IF_cal_z_checkbox_nonequidistant,'value',0);
            set(h.IF_calc_z_edit_input,'enable','off')
        else
            set(h.IF_calc_z_edit_interval_init,'enable','off');
            set(h.IF_calc_z_edit_interval_end,'enable','off');
            set(h.IF_calc_z_edit_step,'enable','off');
            set(h.IF_cal_z_checkbox_nonequidistant,'value',1);
            set(h.IF_calc_z_edit_input,'enable','on')
        end
    end
    function callback_IF_cb_Znonequidistant(hObj,eventadata,h)
        Id=get(hObj,'value');
        if Id==0
            set(h.IF_calc_z_edit_interval_init,'enable','on');
            set(h.IF_calc_z_edit_interval_end,'enable','on');
            set(h.IF_calc_z_edit_step,'enable','on');
            set(h.IF_cal_z_checkbox_equidistant,'value',1);
            set(h.IF_calc_z_edit_input,'enable','off')
        else
            set(h.IF_calc_z_edit_interval_init,'enable','off');
            set(h.IF_calc_z_edit_interval_end,'enable','off');
            set(h.IF_calc_z_edit_step,'enable','off');
            set(h.IF_cal_z_checkbox_equidistant,'value',0);
            set(h.IF_calc_z_edit_input,'enable','on')
        end
    end

    function callback_IF_cb_calc(hObj,eventdata,h)
       Id=get(hObj,'value');
       if Id==1
          set(h.IF_calc_time_edit_interval_init,'enable','on');
          set(h.IF_calc_time_edit_interval_end,'enable','on');
          set(h.IF_calc_time_edit_step,'enable','on');
          set(h.IF_calc_x_edit_interval_init,'enable','on');
          set(h.IF_calc_x_edit_interval_end,'enable','on');
          set(h.IF_calc_x_edit_step,'enable','on');
          set(h.IF_calc_y_edit_interval_init,'enable','on');
          set(h.IF_calc_y_edit_interval_end,'enable','on');
          set(h.IF_calc_y_edit_step,'enable','on');
          set(h.IF_cal_z_checkbox_equidistant,'enable','on');
          set(h.IF_cal_z_checkbox_nonequidistant,'enable','on');
     
          if get(h.IF_cal_z_checkbox_equidistant,'value')==1
          set(h.IF_calc_z_edit_interval_init,'enable','on');
          set(h.IF_calc_z_edit_interval_end,'enable','on');
          set(h.IF_calc_z_edit_step,'enable','on');
          set(h.IF_calc_z_edit_input,'enable','off');
          else
          set(h.IF_calc_z_edit_interval_init,'enable','off');
          set(h.IF_calc_z_edit_interval_end,'enable','off');
          set(h.IF_calc_z_edit_step,'enable','off');    
          set(h.IF_calc_z_edit_input,'enable','on');
          end
          set(h.IF_calc_IDdata_edit,'enable','on');
          set(h.IF_calc_pushbutton,'enable','on');
          set(h.IF_calc_load_push,'enable','off')
          set(h.IF_calc_load_checkbox,'value',0)
       else
          set(h.IF_calc_time_edit_interval_init,'enable','off');
          set(h.IF_calc_time_edit_interval_end,'enable','off');
          set(h.IF_calc_time_edit_step,'enable','off');
          set(h.IF_calc_x_edit_interval_init,'enable','off');
          set(h.IF_calc_x_edit_interval_end,'enable','off');
          set(h.IF_calc_x_edit_step,'enable','off');
          set(h.IF_calc_y_edit_interval_init,'enable','off');
          set(h.IF_calc_y_edit_interval_end,'enable','off');
          set(h.IF_calc_y_edit_step,'enable','off');
          set(h.IF_cal_z_checkbox_equidistant,'enable','off');
          set(h.IF_cal_z_checkbox_nonequidistant,'enable','off');
     
          set(h.IF_calc_z_edit_interval_init,'enable','off');
          set(h.IF_calc_z_edit_interval_end,'enable','off');
          set(h.IF_calc_z_edit_step,'enable','off');  
          set(h.IF_calc_z_edit_input,'enable','off');
         
          set(h.IF_calc_IDdata_edit,'enable','off');
          set(h.IF_calc_pushbutton,'enable','off');
          set(h.IF_calc_load_push,'enable','on')
          set(h.IF_calc_load_checkbox,'value',1)
       end
    end


function callback_IF_cb_calc_load(hObj,eventdata,h)
       Id=get(hObj,'value');
       if Id==0
          set(h.IF_calc_time_edit_interval_init,'enable','on');
          set(h.IF_calc_time_edit_interval_end,'enable','on');
          set(h.IF_calc_time_edit_step,'enable','on');
          set(h.IF_calc_x_edit_interval_init,'enable','on');
          set(h.IF_calc_x_edit_interval_end,'enable','on');
          set(h.IF_calc_x_edit_step,'enable','on');
          set(h.IF_calc_y_edit_interval_init,'enable','on');
          set(h.IF_calc_y_edit_interval_end,'enable','on');
          set(h.IF_calc_y_edit_step,'enable','on');
          set(h.IF_cal_z_checkbox_equidistant,'enable','on');
          if get(h.IF_cal_z_checkbox_equidistant,'value')==1
          set(h.IF_calc_z_edit_interval_init,'enable','on');
          set(h.IF_calc_z_edit_interval_end,'enable','on');
          set(h.IF_calc_z_edit_step,'enable','on');
          set(h.IF_calc_z_edit_input,'enable','off');
          else
          set(h.IF_calc_z_edit_interval_init,'enable','off');
          set(h.IF_calc_z_edit_interval_end,'enable','off');
          set(h.IF_calc_z_edit_step,'enable','off');    
          set(h.IF_calc_z_edit_input,'enable','on');
          end
          set(h.IF_calc_IDdata_edit,'enable','on');
          set(h.IF_calc_pushbutton,'enable','on');
          set(h.IF_calc_load_push,'enable','off')
          set(h.IF_calc_checkbox,'value',1)
       else
          set(h.IF_calc_time_edit_interval_init,'enable','off');
          set(h.IF_calc_time_edit_interval_end,'enable','off');
          set(h.IF_calc_time_edit_step,'enable','off');
          set(h.IF_calc_x_edit_interval_init,'enable','off');
          set(h.IF_calc_x_edit_interval_end,'enable','off');
          set(h.IF_calc_x_edit_step,'enable','off');
          set(h.IF_calc_y_edit_interval_init,'enable','off');
          set(h.IF_calc_y_edit_interval_end,'enable','off');
          set(h.IF_calc_y_edit_step,'enable','off');
          set(h.IF_cal_z_checkbox_equidistant,'enable','off');
          set(h.IF_cal_z_checkbox_nonequidistant,'enable','off');
     
          set(h.IF_calc_z_edit_interval_init,'enable','off');
          set(h.IF_calc_z_edit_interval_end,'enable','off');
          set(h.IF_calc_z_edit_step,'enable','off');  
          set(h.IF_calc_z_edit_input,'enable','off');
        
          set(h.IF_calc_IDdata_edit,'enable','off');
          set(h.IF_calc_pushbutton,'enable','off');
          set(h.IF_calc_load_push,'enable','on')
           set(h.IF_calc_checkbox,'value',0)
       end
    end
    

    function callback_IF_pb_calculation(hObj,eventdata,h)
        inputIF.time.init=get(h.IF_calc_time_edit_interval_init,'userdata');
        inputIF.time.end=get(h.IF_calc_time_edit_interval_end,'userdata');
        inputIF.time.step=get(h.IF_calc_time_edit_step,'userdata');
        inputIF.x.init=get(h.IF_calc_x_edit_interval_init,'userdata');
        inputIF.x.end=get(h.IF_calc_x_edit_interval_end,'userdata');
        inputIF.x.step=get(h.IF_calc_y_edit_step,'userdata');
        inputIF.y.init=get(h.IF_calc_y_edit_interval_init,'userdata');
        inputIF.y.end=get(h.IF_calc_y_edit_interval_end,'userdata');
        inputIF.y.step=get(h.IF_calc_y_edit_step,'userdata');
        inputIF.z.init=get(h.IF_calc_z_edit_interval_init,'userdata');
        inputIF.z.end=get(h.IF_calc_z_edit_interval_end,'userdata');
        inputIF.z.Nz=get(h.IF_calc_z_edit_step,'userdata');
        inputIF.z.input=get(h.IF_calc_z_edit_input,'userdata');
        inputIF.z.Idequidistant=get(h.IF_cal_z_checkbox_equidistant,'value');
        inputIF.data=get(h.IF_project_load_data,'userdata');
        inputIF.fileId=get(h.IF_calc_IDdata_edit,'userdata');
        % assignin('base','inputIF',inputIF)
        if isempty(inputIF.time.init) || inputIF.time.init<inputIF.data.outkinematic.time(1)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify an initial time!')
            uicontrol(h.IF_calc_time_edit_interval_init);
            return;
        end
        if isempty(inputIF.time.end) || inputIF.time.end>inputIF.data.outkinematic.time(end)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a final time!')
            uicontrol(h.IF_calc_time_edit_interval_end);
            return;
        end
        if inputIF.time.init==inputIF.time.end
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify an initial & final time or load simulation data !')
            uicontrol(h.IF_calc_time_edit_interval_init);
            return;
        end
        
        if isempty(inputIF.time.step) || inputIF.time.step<1
             set(h.monitorbox,'foregroundcolor','r','string','>> Specify a time step!')
            uicontrol(h.IF_calc_time_edit_step);
            return;
        end
        
         if isempty(inputIF.x.init) || inputIF.x.init<inputIF.data.dom.X(1)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify an initial x!')
            uicontrol(h.IF_calc_x_edit_interval_init);
            return;
        end
        if isempty(inputIF.x.end) || inputIF.x.end>inputIF.data.dom.X(end)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a final x!')
            uicontrol(h.IF_calc_x_edit_interval_end);
            return;
        end
        if isempty(inputIF.x.step) || inputIF.x.step<1
             set(h.monitorbox,'foregroundcolor','r','string','>> Specify a x step!')
            uicontrol(h.IF_calc_x_edit_step);
            return;
        end
        
        if isempty(inputIF.y.init) || inputIF.y.init<inputIF.data.dom.Y(1)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify an initial y!')
            uicontrol(h.IF_calc_y_edit_interval_init);
            return;
        end
        if isempty(inputIF.y.end) || inputIF.y.end>inputIF.data.dom.Y(end)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a final y!')
            uicontrol(h.IF_calc_y_edit_interval_end);
            return;
        end
        if isempty(inputIF.y.step) || inputIF.y.step<1
             set(h.monitorbox,'foregroundcolor','r','string','>> Specify a y step!')
            uicontrol(h.IF_calc_y_edit_step);
            return;
        end
        if inputIF.z.Idequidistant==1
            if isempty(inputIF.z.init) || inputIF.z.init<min(min(inputIF.data.dom.bathy.profile))
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify an initial z, larger than maximum depth!')
                uicontrol(h.IF_calc_z_edit_interval_init);
                return;
            end
            maxEta=max(max(max(inputIF.data.outkinematic.eta)));
            if isempty(inputIF.z.end) || inputIF.z.end>maxEta
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a final z, less than or equal to maximum elevation (',num2str(maxEta),')[m]')
                uicontrol(h.IF_calc_z_edit_interval_end);
                return;
            end
            if isempty(inputIF.z.Nz) || inputIF.z.Nz<1
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify a number of z discretization!')
                uicontrol(h.IF_calc_z_edit_step);
                return;
            end
        else
             if isempty(inputIF.z.input) 
                set(h.monitorbox,'foregroundcolor','r','string','>> Specify z levels [m]')
                uicontrol(h.IF_calc_z_edit_input);
                return;
            end
            
        end
        
        if isempty(inputIF.fileId) || inputIF.fileId<0
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a ID data file for the output filename.')
            uicontrol(h.IF_calc_IDdata_edit);
            return;
        end
        [IFcalc]=InternalFlowCalc2D(inputIF);
        set(hObj,'userdata',IFcalc);
%         assignin('base','IFcalc',IFcalc)
    end

    function callback_IF_pb_calc_load(hObj,eventdata,h)
        projdir=get(h.IF_project_popup_projdir,'string');
        
        if isempty(projdir)||strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>>Specify a project directory for post-processing');
            return;
        end
        
        
        [file_name,directory]=uigetfile([projdir,'\','*.mat'],'Load Data file');
        
        if directory~=0
            set(h.monitorbox,'foregroundcolor','k','string','>>loading...');
            pause(0.001)
            my_data=load([directory,file_name]);
            my_data.Load=1;
           % assignin('base','my_data',my_data);
            if isempty(my_data)
                set(h.IF_calc_load_text,'string','no data');
                set(h.monitorbox,'foregroundcolor','r','string','>>No data loaded');
            elseif ~isfield(my_data,'IFcalc')
                set(h.IF_calc_load_text,'string','no data');
                set(h.monitorbox,'foregroundcolor','r','string','>>Wrong input file');
            else
                set(hObj,'userdata',my_data);
                set(h.monitorbox,'foregroundcolor','k','string','>>Data loaded');
                set(h.IF_calc_load_text,'string',[file_name]); 
            end
            
        end
    end

    function callback_IF_pop_density_axesopt(hObj,eventdata,h)
       Id=get(hObj,'value');
       if Id==1
       set(h.IF_plot_density_var1_text1,'string','at time:') ;
       set(h.IF_plot_density_var1_text2,'string','[s]') ;
       set(h.IF_plot_density_var2_text1,'string','at z:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;
       elseif Id==2
       set(h.IF_plot_density_var1_text1,'string','at time:') ;
       set(h.IF_plot_density_var1_text2,'string','[s]') ;
       set(h.IF_plot_density_var2_text1,'string','at y:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;    
       elseif Id==3
       set(h.IF_plot_density_var1_text1,'string','at time:') ;
       set(h.IF_plot_density_var1_text2,'string','[s]') ;
       set(h.IF_plot_density_var2_text1,'string','at x:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;    
       elseif Id==4
       set(h.IF_plot_density_var1_text1,'string','at y:') ;
       set(h.IF_plot_density_var1_text2,'string','[m]') ;
       set(h.IF_plot_density_var2_text1,'string','at z:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;    
        elseif Id==5
       set(h.IF_plot_density_var1_text1,'string','at x:') ;
       set(h.IF_plot_density_var1_text2,'string','[m]') ;
       set(h.IF_plot_density_var2_text1,'string','at z:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;   
       elseif Id==6
       set(h.IF_plot_density_var1_text1,'string','at x:') ;
       set(h.IF_plot_density_var1_text2,'string','[m]') ;
       set(h.IF_plot_density_var2_text1,'string','at y:') ;
       set(h.IF_plot_density_var2_text2,'string','[m]') ;   
       end
    end

    function callback_IF_edit_density(hObj,eventdata,h)
       param=str2num(get(hObj,'string'));
       set(hObj,'userdata',param)
       set(h.monitorbox,'foregroundcolor','k','string',['>>'])
       if isempty(param) || length(param)~=1
           set(h.monitorbox,'foregroundcolor','r','string',['>> Specify an input value'])
           uicontrol(hObj);
           return;
       end       
    end


    function callback_pp_setting_cb_level1(hObj,eventdata,h,href1,href2)
       Id=get(hObj,'value');
       Idref=get(href1,'value');
       if Id==1 && Idref==1
           set(href2,'enable','on');
       else
           set(href2,'enable','off');
       end
           
    end
        
    function callback_IF_button_density(hObj,eventdata,h)
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        IdCalcdat=get(h.IF_calc_checkbox,'value');
        if  IdCalcdat==1
            propdata=get(h.IF_calc_pushbutton,'userdata');
        else
            tempdata=get(h.IF_calc_load_push,'userdata');
            propdata=tempdata.IFcalc;
        end
        inputvar.propertiesId=get(h.IF_plot_density_prop_popup,'value');
        inputvar.axesoptId=get(h.IF_plot_density_axesopt_popup,'value');
        inputvar.var1=get(h.IF_plot_density_var1_edit,'userdata');
        inputvar.var2=get(h.IF_plot_density_var2_edit,'userdata');
        inputvar.levelCb=get(h.IF_plot_density_profile_cb_level,'value');
      %  inputvar.quiverCb=get(h.IF_plot_density_profile_cb_quiver,'value');
        
        setting.view.check=get(h.IF_plot_density_profile_setting_cb_view,'value');
        setting.view.param=get(h.IF_plot_density_profile_setting_edit_view,'userdata');
        setting.clim.check=get(h.IF_plot_density_profile_setting_cb_clim,'value');
        setting.clim.param=get(h.IF_plot_density_profile_setting_edit_clim,'userdata');
        setting.xlim.check=get(h.IF_plot_density_profile_setting_cb_xlim,'value');
        setting.xlim.param=get(h.IF_plot_density_profile_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.IF_plot_density_profile_setting_cb_ylim,'value');
        setting.ylim.param=get(h.IF_plot_density_profile_setting_edit_ylim,'userdata');
        setting.coarse.check=get(h.IF_plot_density_profile_setting_cb_coarse,'value');
        setting.coarse.param=get(h.IF_plot_density_profile_setting_edit_coarse,'userdata');
        setting.level.check=get(h.IF_plot_density_profile_setting_cb_level,'value');
        setting.level.param=get(h.IF_plot_density_profile_setting_edit_level,'userdata');
        setting.savefig.check=get(h.IF_plot_density_profile_setting_cb_savefig,'value');
        format=cellstr(get(h.IF_plot_density_profile_setting_popup_savefig,'string'));
        setting.savefig.format=format(get(h.IF_plot_density_profile_setting_popup_savefig,'value'));
        
        colorfig=cellstr(get(h.IF_plot_density_profile_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.IF_plot_density_profile_setting_popup_colormap,'value'));
        
        filename=get(h.IF_project_edit_name,'string');
        projdir=get(h.IF_project_popup_projdir,'string');
        
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.IF_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.IF_project_popup_projdir);
            return;
        end
        
        if isempty(propdata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a calculated internal properties data')
            uicontrol(h.IF_calc_pushbutton);
            return;
        end
     
        X=propdata.x;Y=propdata.y;Z=propdata.z;
        T=propdata.time;
        if isempty(inputvar.var1) || length(inputvar.var1)~=1
          set(h.monitorbox,'foregroundcolor','r','string',['>> Specify an input value'])
            uicontrol(h.IF_plot_density_var1_edit);
            return;       
        end
        if isempty(inputvar.var2) || length(inputvar.var2)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify an input value'])
            uicontrol(h.IF_plot_density_var2_edit);
            return;     
        end
        
        if inputvar.axesoptId==1
            if inputvar.var1<T(1)  || inputvar.var1>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in interval (',num2str(T(1)),';',num2str(T(end)),')']);
                uicontrol(h.IF_plot_density_var1_edit);
                return;
            end
             if inputvar.var2<Z(1)  || inputvar.var2>Z(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a Z value in interval (',num2str(Z(1)),';',num2str(Z(end)),')']);
                uicontrol(h.IF_plot_density_var2_edit);
                return;
            end
        elseif inputvar.axesoptId==2
            if inputvar.var1<T(1)  || inputvar.var1>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in interval (',num2str(T(1)),';',num2str(T(end)),')']);
                uicontrol(h.IF_plot_density_var1_edit);
                return;
            end
             if inputvar.var2<Y(1)  || inputvar.var2>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a Y value in interval (',num2str(Y(1)),';',num2str(Y(end)),')']);
                uicontrol(h.IF_plot_density_var2_edit);
                return;
             end
          elseif inputvar.axesoptId==3
            if inputvar.var1<T(1)  || inputvar.var1>T(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in interval (',num2str(T(1)),';',num2str(T(end)),')']);
                uicontrol(h.IF_plot_density_var1_edit);
                return;
            end
             if inputvar.var2<X(1)  || inputvar.var2>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a X value in interval (',num2str(X(1)),';',num2str(X(end)),')']);
                uicontrol(h.IF_plot_density_var2_edit);
                return;
            end  
        end
        setting.workdir=[projdir,'\',filename,'\'];
        if ~isdir(setting.workdir)
            mkdir(setting.workdir);
        end
        fun_plotting_density_IF2D(h,inputvar,propdata,setting);
    end


    function callback_IF_pop_animdensity_axesopt(hObj,eventdata,h)
        Id=get(hObj,'value');
        if Id==1
            set(h.IF_anim_density_var2_text1,'string','at z:') ;
            set(h.IF_anim_density_var2_text2,'string','[m]') ;
        elseif Id==2
            set(h.IF_anim_density_var2_text1,'string','at y:') ;
            set(h.IF_anim_density_var2_text2,'string','[m]') ;
        elseif Id==3
            set(h.IF_anim_density_var2_text1,'string','at x:') ;
            set(h.IF_anim_density_var2_text2,'string','[m]') ;
        end
    end

    function calback_IF_anim_plot(hObj,eventdata,h)
        if get(h.pp_anim_pause,'userdata')==1
            set(h.monitorbox,'foregroundcolor','r','string','>> Please press stop button.')
            return;
        end
        
        set(h.monitorbox,'foregroundcolor','k','string','>>')
        
        IdCalcdat=get(h.IF_calc_checkbox,'value');
        if  IdCalcdat==1
            propdata=get(h.IF_calc_pushbutton,'userdata');
        else
            tempdata=get(h.IF_calc_load_push,'userdata');
            propdata=tempdata.IFcalc;
        end
        inputvar.propertiesId=get(h.IF_anim_density_prop_popup,'value');
        inputvar.axesoptId=get(h.IF_anim_density_axesopt_popup,'value');
        inputvar.var2=get(h.IF_anim_density_var2_edit,'userdata');
        inputvar.levelCb=get(h.IF_anim_density_profile_cb_level,'value');

        setting.view.check=get(h.IF_anim_density_setting_cb_view,'value');
        setting.view.param=get(h.IF_anim_density_setting_edit_view,'userdata');
        setting.clim.check=get(h.IF_anim_density_setting_cb_clim,'value');
        setting.clim.param=get(h.IF_anim_density_setting_edit_clim,'userdata');
        setting.xlim.check=get(h.IF_anim_density_setting_cb_xlim,'value');
        setting.xlim.param=get(h.IF_anim_density_setting_edit_xlim,'userdata');
        setting.ylim.check=get(h.IF_anim_density_setting_cb_ylim,'value');
        setting.ylim.param=get(h.IF_anim_density_setting_edit_ylim,'userdata');
        setting.tlim.check=get(h.IF_anim_density_setting_cb_tlim,'value');
        setting.tlim.param=get(h.IF_anim_density_setting_edit_tlim,'userdata');
        setting.coarse.check=get(h.IF_anim_density_setting_cb_coarse,'value');
        setting.coarse.param=get(h.IF_anim_density_setting_edit_coarse,'userdata');
        setting.level.check=get(h.IF_anim_density_setting_cb_level,'value');
        setting.level.param=get(h.IF_anim_density_setting_edit_level,'userdata');
        setting.saveanim.check=get(h.IF_anim_density_setting_cb_saveanim,'value');
        setting.gifset.param=get(h.IF_anim_density_setting_edit_gifset,'userdata');
        
        colorfig=cellstr(get(h.IF_anim_density_setting_popup_colormap,'string'));
        setting.colormap=colorfig(get(h.IF_anim_density_setting_popup_colormap,'value'));
        
        filename=get(h.IF_project_edit_name,'string');
        projdir=get(h.IF_project_popup_projdir,'string');
        
        if isempty(filename)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project name!')
            uicontrol(h.IF_project_edit_name);
            return;
        end
        
        if isempty(projdir) || strcmpi(projdir,'--')
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a project directory!')
            uicontrol(h.IF_project_popup_projdir);
            return;
        end
        
        if isempty(propdata)
            set(h.monitorbox,'foregroundcolor','r','string','>> Specify a calculated internal properties data')
            uicontrol(h.IF_calc_pushbutton);
            return;
        end
     
        X=propdata.x;Y=propdata.y;Z=propdata.z;
        T=propdata.time;
        
        if isempty(inputvar.var2) || length(inputvar.var2)~=1
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify an input value'])
            uicontrol(h.IF_anim_density_var2_edit);
            return;     
        end
        
        if inputvar.axesoptId==1
             if inputvar.var2<Z(1)  || inputvar.var2>Z(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a Z value in interval (',num2str(Z(1)),';',num2str(Z(end)),')']);
                uicontrol(h.IF_anim_density_var2_edit);
                return;
            end
        elseif inputvar.axesoptId==2
             if inputvar.var2<Y(1)  || inputvar.var2>Y(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a Y value in interval (',num2str(Y(1)),';',num2str(Y(end)),')']);
                uicontrol(h.IF_anim_density_var2_edit);
                return;
             end
          elseif inputvar.axesoptId==3
             if inputvar.var2<X(1)  || inputvar.var2>X(end)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a X value in interval (',num2str(X(1)),';',num2str(X(end)),')']);
                uicontrol(h.IF_anim_density_var2_edit);
                return;
            end  
        end
        
         if setting.view.check==1
            if isempty(setting.view.param)||length(setting.view.param)>2||length(setting.view.param)<1
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.IF_anim_density_setting_edit_view);
                return;
            end
            if length(setting.view.param)==1 && (setting.view.param>3||setting.view.param<2)
                set(h.monitorbox,'foregroundcolor','r','string',['>> Specify view: 2 or 3',...
                    ' or [az;el], az is azimuth and el is vertical elevation of the view point'])
                uicontrol(h.IF_anim_density_setting_edit_view);
                return;
            end
        end
        if setting.clim.check==1
            if isempty(setting.clim.param)||length(setting.clim.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a colorbar axes limit [zmin;zmax]'])
                uicontrol(h.IF_anim_density_setting_edit_clim);
                return;
            end
        end
        
        if setting.xlim.check==1
            if isempty(setting.xlim.param)||length(setting.xlim.param)~=2 || setting.xlim.param(1)>setting.xlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x axis limit [xmin;xmax]'])
                uicontrol(h.IF_anim_density_setting_edit_xlim);
                return;
            end
            if inputvar.axesoptId==1 ||  inputvar.axesoptId==2
              Xvar1=X(1);XvarN=X(end);
             elseif inputvar.axesoptId==3
              Xvar1=Y(1);XvarN=Y(end);  
            end
            
            if setting.xlim.param(1)<Xvar1 || setting.xlim.param(1)>XvarN
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(Xvar1),';',num2str(XvarN),'] (m)']);
                uicontrol(h.IF_anim_density_setting_edit_xlim);
                return;
            end
            
            if setting.xlim.param(2)<Xvar1 || setting.xlim.param(2)>XvarN
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify x lim in the interval [',num2str(Xvar1),';',num2str(XvarN),'] (m)']);
                uicontrol(h.IF_anim_density_setting_edit_xlim);
                return;
            end
            
        end
        
         if setting.ylim.check==1
           if isempty(setting.ylim.param)||length(setting.ylim.param)~=2 || setting.ylim.param(1)>setting.ylim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y axis limit [ymin;ymax]'])
                uicontrol(h.IF_anim_density_setting_edit_ylim);
                return;
           end
            
            if inputvar.axesoptId==1 
              Yvar1=Y(1);YvarN=Y(end);
             elseif inputvar.axesoptId==2 || inputvar.axesoptId==3
              Yvar1=Z(1);YvarN=Z(end);  
            end
            
            if setting.ylim.param(1)<Yvar1 || setting.ylim.param(1)>YvarN
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Yvar1),';',num2str(YvarN),'] (m)']);
                uicontrol(h.IF_anim_density_setting_edit_ylim);
                return; 
            end
            
            if setting.ylim.param(2)<Yvar1 || setting.ylim.param(2)>YvarN
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify y lim in the interval [',num2str(Yvar1),';',num2str(YvarN),'] (m)']);
                uicontrol(h.IF_anim_density_setting_edit_ylim);
                return; 
            end
        end
        
        if setting.tlim.check==1
            if isempty(setting.tlim.param)||length(setting.tlim.param)~=2 || setting.tlim.param(1)>setting.tlim.param(2)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a time limit [tmin;tmax]'])
                uicontrol(h.IF_anim_density_setting_edit_tlim);
                return;
            end
           
            if  any(setting.tlim.param<T(1)-0.0001) || any(setting.tlim.param>T(end)+0.0001)
            set(h.monitorbox,'foregroundcolor','r','string',['>> Specify a time in the interval ['...
                ,num2str(T(1)),';',num2str(T(end)),'] (s)'])
            uicontrol(h.IF_anim_density_setting_edit_tlim);
            return;
            end
        end
        
        if setting.coarse.check==1
            if length(setting.coarse.param)~=1
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor'])
                uicontrol(h.IF_anim_density_setting_edit_coarse);
                return;
            end
            if  setting.coarse.param<=0
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify a coarse factor that is > 0'])
                uicontrol(h.IF_anim_density_setting_edit_coarse);
                return;
            end
        end
        if setting.level.check==1
            if isempty(setting.level.param)
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify level/quiver option'])
                uicontrol(h.IF_anim_density_setting_edit_level);
                return;
            end
        end
        if setting.saveanim.check==1
            if isempty(setting.gifset.param)||length(setting.gifset.param)~=2
                set(h.monitorbox,'foregroundcolor','r','string',['>>Specify gif parameters: [delay factor;number of loop]'])
                uicontrol(h.IF_anim_density_setting_edit_gifset);
                return;
            end
        end
         setting.workdir=[projdir,'\',filename,'\'];
        if ~isdir(setting.workdir)
            mkdir(setting.workdir);
        end
        funIF_anim_density(h,inputvar,propdata,setting);
    end
%%%%%%%%%%%%%%%%%%%%%functions for layout%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function menuSelectFcn(~,~,h)
        % get selected node
               
        nodes = mtree.getSelectedNodes;
        if isempty(nodes), return; end
        n = nodes(1);
        
        if strcmp(n,'com.mathworks.hg.peer.UITreeNode:Project') == 1
            menupanel.SelectedChild = 1;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Model') == 1
            menupanel.SelectedChild = 2;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Spatial') == 1
            menupanel.SelectedChild = 3;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Domain Area') == 1
            menupanel.SelectedChild = 4;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Boundaries') == 1
            menupanel.SelectedChild = 5;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Wall') == 1
            menupanel.SelectedChild = 6;
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Fourier Bdry.') == 1
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            menupanel.SelectedChild = 7;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Bottom') == 1
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            menupanel.SelectedChild = 8;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Wave Input') == 1
            menupanel.SelectedChild = 9;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Initial Condition') == 1
            menupanel.SelectedChild = 10;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Influx') == 1
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            menupanel.SelectedChild = 11;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Boundary Assimilation') == 1
            menupanel.SelectedChild = 12;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Time') == 1
            menupanel.SelectedChild = 13;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Options') == 1
            menupanel.SelectedChild = 14;            
        end
    end

    function mySelectFcn_preview(~,~)
        nodes = mtree_prev.getSelectedNodes;
        if isempty(nodes), return; end
        n = nodes(1);
        
        %         % only consider a leaf node (skip folders)
        %         if ~n.isLeaf, return; end
        
        if strcmp(n,'com.mathworks.hg.peer.UITreeNode:Dispersion') == 1
            menupanel.SelectedChild = 15;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Spatial') == 1
            menupanel.SelectedChild = 16;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Wave-input') == 1
            menupanel.SelectedChild = 17;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Log file') == 1
            menupanel.SelectedChild = 18;
        else
            %menupanel.SelectedChild = 1;
        end
    end

    function mySelectFcn_postproc(~,~,h)
        % get selected node
        nodes = mtree_postproc.getSelectedNodes;
        if isempty(nodes), return; end
        n = nodes(1);
      
        %         % only consider a leaf node (skip folders)
        %         if ~n.isLeaf, return; end
        
        if strcmp(n,'com.mathworks.hg.peer.UITreeNode:Project') == 1
            menupanel.SelectedChild = 19;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Plotting') == 1
            menupanel.SelectedChild = 20;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Density plots') == 1
            menupanel.SelectedChild = 21;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Profile') == 1
            menupanel.SelectedChild = 22;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Signal') == 1
            menupanel.SelectedChild = 23;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Statistic') == 1
            menupanel.SelectedChild = 24;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Line plots') == 1
            menupanel.SelectedChild = 25;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Profile ') == 1
            menupanel.SelectedChild = 26;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Buoy ') == 1
            menupanel.SelectedChild = 27;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Statistic ') == 1
            menupanel.SelectedChild = 28;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Hamiltonian & Momentum') == 1
            menupanel.SelectedChild = 29;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Breaking events') == 1
            menupanel.SelectedChild = 30;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Max. Avg. Temp. Area Amplitude') == 1
            menupanel.SelectedChild = 31;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Extreme events') == 1
            menupanel.SelectedChild = 32;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Quantitative (Buoy)') == 1
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            menupanel.SelectedChild = 33;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Animation') == 1
            menupanel.SelectedChild = 34;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Density') == 1
            menupanel.SelectedChild = 35;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Line') == 1
            menupanel.SelectedChild = 36;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Validation') == 1
            menupanel.SelectedChild = 37;
       elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Buoy') == 1
            menupanel.SelectedChild = 38;
       elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Quantitative') == 1
            guiAB2dwaveResizeFcn([],[],h);%resizing table for fixing bug (table disappear) in Matlab 2017a
            menupanel.SelectedChild = 39;
        end
    end

    function mySelectFcn_IF(~,~,h)
        % get selected node
        nodes = mtree_IF.getSelectedNodes;
        if isempty(nodes), return; end
        n = nodes(1);
        if strcmp(n,'com.mathworks.hg.peer.UITreeNode:Project') == 1
            menupanel.SelectedChild = 40;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Calculation') == 1
            menupanel.SelectedChild = 41;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Plotting') == 1
            menupanel.SelectedChild = 42;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Density plots') == 1
            menupanel.SelectedChild = 43;
        elseif strcmp(n,'com.mathworks.hg.peer.UITreeNode:Animation') == 1
            menupanel.SelectedChild = 44;    
        end
    end
end