function plot_sites_PSHA(handles)

delete(findall(handles.ax1,'tag','siteplot'));

if isempty(handles.h.p)
    handles.po_sites.Enable='off';
    handles.po_sites.Value=0;
    handles.runMRE.Enable='off';
    handles.runMRD.Enable='off';
    handles.ShakeField.Enable='off';
    handles.ExploreDeagg.Enable='off';
    return
end

site   = cell2mat(handles.h.p(:,2:3));
Nsites = length(site(:,1));
c = repmat([0      0.44706      0.74118],Nsites,1);

switch handles.po_sites.Value
    case 0, VIS='off';
    case 1, VIS='on';
end
scatter(handles.ax1,site(:,2),site(:,1),[],c,'filled','markeredgecolor','k','tag','siteplot','ButtonDownFcn',{@site_click_PSHA;handles;1},'visible',VIS);
handles.po_sites.Enable='on';

if isfield(handles,'model')
    switch handles.model(1).id
        case 'USGS_NHSM_2008'
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='off';
            handles.SeismicHazardTools.Enable='off';
            handles.runMRE.Enable='on';
            handles.runMRD.Enable='off';
        case 'USGS_NHSM_2014'
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='off';
            handles.SeismicHazardTools.Enable='off';
            handles.runMRE.Enable='on';
            handles.runMRD.Enable='off';
        otherwise
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='on';
            handles.SeismicHazardTools.Enable='on';
            handles.runMRE.Enable='on';
            handles.runMRD.Enable='on';
    end
else
    handles.Shakefield.Enable='off';
    handles.Run.Enable='off';
end

if size(handles.h.p,1)>0
    handles.site_menu.String = handles.h.p(:,1);
    handles.site_menu_psda.String = handles.h.p(:,1);
    if isempty(handles.site_selection)
        handles.site_menu.Value  = 1;
        handles.site_menu_psda.Value  = 1;
    else
        handles.site_menu.Value       = handles.site_selection(1);
        handles.site_menu_psda.Value  = handles.site_selection(1);
    end
    ch=findall(handles.ax1,'tag','siteplot');
    uistack(ch, 'top');
end


