function plot_hazardmap_PSHA(handles)

MapOptions = handles.HazOptions;

% This piece of code prevents the patches to be generated for paths elements
t      = handles.h.t;
ispath = regexp(t(:,1),'path');
for i=1:length(ispath)
    if isempty(ispath{i}),ispath{i}=0;end
end
ispath=cell2mat(ispath);
t(ispath==1,:)=[];
if size(t,1)==0
    return
end

% ------------------------------------------------------------------------

delete(findall(handles.FIGSeismicHazard,'Tag','Colorbar'));
ch=findall(handles.ax1,'tag','gmap'); uistack(ch,'bottom');
nt = size(handles.h.t,1);
for i=1:nt
    ch = findall(handles.ax1,'tag',num2str(i));
    delete(ch);
end

ch=findall(handles.ax1,'Tag','satext');delete(ch);
val      = handles.IM_select.Value;

lambda   = nansum(handles.MRE(:,:,val,:,:),4);
lambda   = permute(lambda,[5,1,2,3,4]);

switch MapOptions.mod
    case 1 %weithed average of hazard
        if MapOptions.avg(1)==1 % default weights
            weights = handles.sys.WEIGHT(:,4);
            lambda  = prod(bsxfun(@power,lambda,weights),1);
        end
        
        if MapOptions.avg(2)==1
            weights  = MapOptions.rnd;
            lambda  = prod(bsxfun(@power,lambda,weights),1);
        end
        
        if MapOptions.avg(3)==1
            weights  = [];
            Per      = MapOptions.avg(4);
            lambda   = nansum(handles.MRE(:,:,val,:,:),4);
            lambda   = permute(lambda,[5,1,2,3,4]);
            lambda   = prctile(lambda,Per,1);
        end
        
    case 2 %single branch
        ptr      = MapOptions.sbh(1);
        lambda   = nansum(handles.MRE(:,:,val,:,:),4);
        lambda   = permute(lambda,[5,1,2,3,4]);
        lambda   = lambda(ptr,:,:);
end
lambda   = permute(lambda,[2,3,1]);
im       = handles.opt.im;


pall={'parula','autumn','bone','colorcube','cool','copper','flag','gray','hot','hsv','jet','lines','pink','prism','spring','summer','white','winter'};
pallet = pall{MapOptions.map(2)};
set(handles.FIGSeismicHazard,'colormap',feval(pallet));
hazard   = 1/MapOptions.map(1);
logh     = log(hazard);
logHIM   = log(lambda);


umax = -inf;
umin = inf;
label = handles.h.p(:,1);
p     = cell2mat(handles.h.p(:,[3,2]));

if ~isempty(handles.h.shape)
    active = vertcat(handles.h.shape.active);
    active = find(active);
end
ii     = 1;
handles.SaText = text(NaN,NaN,'','parent',handles.ax1,'Tag','satext');

% interpolation oh hazard curves
Nsites = size(handles.h.p,1);
v   = zeros(Nsites,1);
for j=1:Nsites
    xxx =logHIM(j,:);
    if max(xxx)<logh || length(xxx)<2
        v(j)  = 0;
    else
        v(j)  = interp1(xxx(~isinf(xxx)),im(~isinf(xxx)),logh,'linear');
    end
end
handles.uIM=v;
% -------------

for i=1:size(t,1)
    vptr=regexp(label,t{i});
    for j=1:length(vptr)
        if isempty(vptr{j})
            vptr{j}=0;
        end
    end
    vptr = find(cell2mat(vptr));
    x    = p(vptr,1);
    y    = p(vptr,2);
    u    = v(vptr);
    F =  scatteredInterpolant(x,y,u,'linear','none');
    handles.TriScatterd{i}=F;
    if regexp(t{i,1},'grid')
        conn = t{i,2};
        gps  = [x,y];
        in   = 1:size(gps,1);
    end
    
    if regexp(t{i,1},'shape')
        ind = active(ii);
        xl=handles.h.shape(ind).Lon';
        yl=handles.h.shape(ind).Lat';
        faces = handles.h.shape(ind).faces;
        pv =[xl,yl];
        gps  = zeros(0,2);
        conn = zeros(0,3);
        offset=0;
        for count=1:size(faces,1)
            indface = faces(count,:);
            indface(isnan(indface))=[];
            pvface = pv(indface,:);
            [gps_i,conn_i]=triangulate_maps(pvface,[x,y]);
            gps = [gps;gps_i]; %#ok<AGROW>
            conn = [conn;conn_i+offset]; %#ok<AGROW>
            offset=size(gps,1);
        end
        u = F(gps(:,1),gps(:,2));
        in = inpolygon(gps(:,1),gps(:,2),xl,yl);
        ii=ii+1;
    end
    handles.shading(i) = patch(...
        'parent',handles.ax1,...
        'vertices',gps,...
        'faces',conn,...
        'facevertexcdata',u,...
        'facecol','interp',...
        'edgecol','none',...
        'linewidth',0.5,...
        'facealpha',0.7,...
        'Tag',num2str(i),...
        'ButtonDownFcn',{@site_click_PSHA,handles,2},...
        'visible','on');
    uistack(handles.shading(i),'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    
    uin = u;
    uin = uin(in);
    umin = min(min(uin),umin);
    umax = max(max(uin),umax);
    
end

delete(findall(handles.ax1,'tag','siteplot'));
caxis([umin umax])
switch handles.po_sites.Value
    case 0, VIS='off';
    case 1, VIS='on';
end
scatter(handles.ax1,p(:,1),p(:,2),[],v,'filled','markeredgecolor','k','tag','siteplot','ButtonDownFcn',{@site_click_PSHA;handles;1},'visible',VIS);
handles.colorbar=colorbar('peer',handles.ax1,'location','eastoutside','position',[0.94 0.16 0.02 0.65],'ylim',[umin,umax]);
set(get(handles.colorbar,'Title'),'String',handles.IM_select.String{handles.IM_select.Value})
handles.ax1.ButtonDownFcn={@clear_satxt;handles};

% restores map
h = findall(handles.ax1,'tag','gmap');
if ~isempty(h)
    h.ButtonDownFcn={@clear_satxt;handles};
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
end

function clear_satxt(hObject,eventdata,handles) %#ok<*INUSL,*INUSD>
ch=findall(handles.ax1,'Tag','satext');delete(ch);

