function plot_PSDA_pce(handles)

if ~isfield(handles,'sys'), return;end
if ~isfield(handles,'lambdaPCE'),return;end

delete(findall(handles.ax1,'type','line'));
handles.ax1.NextPlot='add';


param    = handles.PCEOptions;
site_ptr = handles.pop_site.Value;
d = handles.d;
model_ptr = param(5);
lambdaPCE2 = nansum(handles.lambdaPCE(:,site_ptr,:,:,model_ptr),4);
lambdaPCE2 = permute(lambdaPCE2,[1 3 2]);
lambdaPCE2(lambdaPCE2<0)=nan;

% plot all scenarios
% notnan =~isnan(sum(lambdaPCE2,2));
notnan = true(size(lambdaPCE2,1),1);

switch param(4)
    case 0
        plot(handles.ax1,d,lambdaPCE2(notnan,:),'-','color',[1 1 1]*0.8,'visible','off','HandleVisibility','off');
    case 1
        plot(handles.ax1,nan,nan,'-','color',[1 1 1]*0.8,'visible','off','HandleVisibility','on');
        plot(handles.ax1,d,lambdaPCE2(notnan,:),'-','color',[1 1 1]*0.8,'visible','on' ,'HandleVisibility','off');
        Leg1 = 'PCE simulations';
end


% plot median
handles.ax1.ColorOrderIndex=1;
if param(1)
    y = exp(nanmean(log(lambdaPCE2),1));
    plot(handles.ax1,d,y,'.-');
    Leg2='Median';
end

% plot percentile
if param(2)
    Per = param(3);
    y   = prctile(lambdaPCE2,Per,1);
    plot(handles.ax1,d,y,'.-');
    Leg2=sprintf('Percentile %g',Per);
end

switch param(4)
    case 0, Leg=legend(Leg2);
    case 1, Leg=legend({Leg1,Leg2});
end
Leg.FontSize=8;
Leg.EdgeColor=[1 1 1];
Leg.Location='SouthWest';
Leg.Tag='hazardlegend';
set(handles.ax1,'layer','top')

