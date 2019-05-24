function figurecompare_uimenu(~, ~,ax)

fig = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1000);


if isempty(fig)
    fig=figure(1000);
    copyobj(ax,fig);
    % hold on
    set(gca,'ActivePositionProperty','outerposition')
    set(gca,'Units','normalized')
    set(gca,'OuterPosition',[0 0 1 1])
    set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
    XL=get(gca,'xlim');
    YL=get(gca,'ylim');
    
    ch = findall(gcf,'tag','redballs');delete(ch);
    ch = findall(gcf,'tag','nanvalues');delete(ch);
    ch = get(gca,'title'); delete(ch);
    
    if strcmp(get(gca,'xscale'),'linear') && strcmp(get(gca,'yscale'),'linear')
        f = uimenu('Label','akZoom');
        uimenu(f,'Label','Zoom X Y','Callback',{@switchzoom,{0,XL,YL}});
        uimenu(f,'Label','Zoom X','Callback',{@switchzoom,{1,XL,YL}});
        uimenu(f,'Label','Zoom Y','Callback',{@switchzoom,{2,XL,YL}});
    end
    legend('Model 1')
else
    newlines = findall(ax,'type','line');
    copyobj(newlines,findobj(1000,'type','axes'));
%     oldlines = findall(fig,'type','line');
%     Nlines   = length(oldlines);
%     co       = get(gca,'colororder');
    L=findall(fig,'type','legend');
    Str=L.String;
    
    Str{end+1}=['Model ',num2str(length(Str)+1)];
    delete(L);
    nax=findall(fig,'type','axes');
    legend(nax,Str);
    
    % this brings up Figure 1000
    figure(1000);
    if length(Str)==2
        f = uimenu('Label','Data');
        uimenu(f,'Label','Flip','Callback',{@dataflip,fig});
    end
    
    h  = findall(fig,'type','line');
    co = [         0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
    for i=1:length(h)
        ind = mod(i,7);
        if ind==0,ind=7;end
        h(end+1-i).Color=co(ind,:);
    end
end

function dataflip(~,~,fid)

linea = findall(fid,'type','line');
uistack(linea(end),'top');



