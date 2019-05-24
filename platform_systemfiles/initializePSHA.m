function handles=initializePSHA(handles)

delete(findall(handles.FIGSeismicHazard,'type','legend'));
delete(findall(handles.FIGSeismicHazard,'tag','Colorbar'));
delete(findall(handles.ax1,'type','line'));
delete(findall(handles.ax1,'type','patch'));
delete(findall(handles.ax1,'tag','satext'));
delete(findall(handles.ax1,'tag','gmap'));
delete(findall(handles.ax1,'tag','patchtxt'));
delete(findall(handles.ax1,'type','text'));
delete(findall(handles.ax2,'type','line'));
handles.switchmode.CData=handles.form1;

if ~isfield(handles,'sys')
    load pshatoolbox_RealValues opt
    handles.opt=opt;
    handles.Boundary_check.Enable='off'; handles.Boundary_check.Value = 0;
    handles.Layers_check.Enable  ='off'; handles.Layers_check.Value   = 0;
    handles.FIGSeismicHazard.Name='Seismic Hazard';
    handles.site_menu.Value=1;
    handles.site_menu.String={''};
    handles.site_menu_psda.Value=1;
    handles.site_menu_psda.String={''};
    handles.IM_select.Value=1;
    handles.IM_select.String={''};
else
    handles = rmfield(handles,{'sys','model'});
end

if handles.opt.ellipsoid.Code==0
    xlabel(handles.ax1,'X (km)','fontsize',10,'tag','xlabel')
    ylabel(handles.ax1,'Y (km)','fontsize',10,'tag','ylabel')
else
    xlabel(handles.ax1,'Longitude (�)','fontsize',10,'tag','xlabel')
    ylabel(handles.ax1,'Latitude (�)','fontsize',10,'tag','ylabel')
end

%--- initialization ---------------------------------------------------
set(handles.ax1,'nextplot','add','box','on','dataaspectratio',[1 1 1],'XGrid','off','YGrid','off');
handles.MapOptionsButton.Enable='inactive';
handles.SeismicHazardTools.Enable = 'off';
handles.Edit.Enable           = 'off';
handles.PSHAMenu.Enable       = 'off';
handles.DSHAMenu.Enable       = 'off';
handles.HazardMap.Enable      = 'off';
handles.po_grid.Enable        = 'on';  
handles.po_sources.Enable     = 'off';
handles.SourceLabels.Enable   = 'off';
handles.po_sourcemesh.Enable  = 'off';
handles.po_region.String      = ' ';
handles.po_region.Enable      = 'off';
handles.po_region.Visible     = 'off';
handles.Enable_Deaggregation.Checked = 'off';
handles.Run_shakefield.Enable = 'off';
handles.ComputeKmeans.Enable  = 'off';
handles.shakepanel.Visible    = 'off';
handles.po_grid.Value         = 0;
handles.po_sources.Value      = 1;
handles.SourceLabels.Value    = 0;
handles.po_sourcemesh.Value   = 0;
handles.po_googleearth.Value  = 1;
handles.PSHApannel.Visible='off';
handles.DSHApannel.Visible='off';

DATA=load('pshatoolbox_emptyGEmap.mat');
handles.ax1.XLim = DATA.XLIM;
handles.ax1.YLim = DATA.YLIM;
gmap             = image(DATA.xx,DATA.yy,DATA.cc, 'Parent', handles.ax1,'Tag','gmap','Visible','on');
handles.ax1.YDir = 'normal';
uistack(gmap,'bottom');
handles.ax1.Layer='top';

%--- AX2 properties -----------------------------------------------------
set(handles.ax2,...
    'box','on',...
    'xscale','log',...
    'yscale','log')
xlabel(handles.ax2,'IM','fontsize',10)
ylabel(handles.ax2,'Mean Rate of Exceedance','fontsize',10)

%% platform data structure
handles.im             = [];
handles.MRE            = [];
handles.MREPCE         = cell(0,0);
handles.deagg          = [];
handles.shakefield     = [];
handles.scenarios      = cell(0,9);
handles.tessel         = {[],zeros(0,2)};
handles.Zscenario      = 1;
handles.krate          = []; %k-mean cluster rate
handles.kY             = []; %k-mean cluster scenarios
handles.TriScatterd    = cell(0,1);
handles.h.p            = cell(0,4);
handles.h.t            = cell(0,2);
handles.h.shape        = [];
handles.h.lmax         = 50;
handles.h.Vs30         = [];
handles.site_selection = [];
handles.plothazmethod  = '';
handles.HazOptions     = struct('mod',1,'avg',[1 0 0 50 0],'sbh',[1 0 0 1],'dbt',[1 0 0 50],'map',[457 1],'pce',[0 1 50],'rnd',1);
handles.isREGULAR      = [];
handles.isPCE          = [];
handles.platformMode   = 1;
handles.TT             = [];
handles.Y              = [];
handles.hdist          = [];
handles.L              = [];
handles.optkm          = [];
plot_sites_PSHA(handles);
