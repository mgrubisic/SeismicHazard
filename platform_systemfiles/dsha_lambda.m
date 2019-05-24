function[im,lambdaIS] = dsha_lambda(handles,~)

IMptr      = handles.pop_field.Value;
siteptr    = handles.site_menu_psda.Value;
Nsim       = str2double(handles.NumSim.String);
opt        = handles.opt;
im         = handles.opt.im;
Nim        = size(im,1);
IMs        = [opt.IM1;opt.IM2];
Nsites     = size(handles.h.p,1);
IND        = siteptr + Nsites*(IMptr-1);

%% EMPIRICAL HAZARD CURVE FROM DSHA ANALYSIS
lambdaIS   = zeros(Nim,1);
Nscen      = length(handles.shakefield);
Y_IS       = zeros(Nscen,Nsim);

Nz        = size(handles.L,1);
Nsim      = str2double(handles.NumSim.String);
for i=1:Nscen
    mu     = handles.shakefield(i).mulogIM(IND);
    Li     = handles.L(IND,:,handles.shakefield(i).Lptr);
    Zi     = normrnd(0,1,[Nz,Nsim]);
    Y_IS(i,:) = exp(mu+Li*Zi);
end
Y_IS    = Y_IS(:);
rate_IS = prod(cell2mat(handles.scenarios(:,10:12)),2);
rate_IS = repmat(rate_IS/Nsim,1,Nsim);
rate_IS = rate_IS(:);

for i=1:Nim
    lambdaIS(i)=sum(rate_IS(Y_IS>=im(i)));
end

%% kmean Clusters
if ~isempty(handles.krate)
    lambdaKM = zeros(Nim,1);
    kY       = handles.kY(:,IND);
    krate    = handles.krate;
    for i=1:Nim
        lambdaKM(i)=sum(krate(kY>=im(i)));
    end
end

%% SEISMIC HAZARD CURVE FROM PSHA ANALYSIS (FOR COMPARISON PURPOSES)
model_ptr  = str2double(unique(handles.scenarios(:,1)));
model      = handles.model(model_ptr);
[~,source_ptr] = intersect({model.source.label},unique(handles.scenarios(:,2)));
IM         = IMs(IMptr);
site       = cell2mat(handles.h.p(siteptr,2:4));
Vs30       = handles.h.Vs30(siteptr);
Nsource    = length(model.source);
MRE        = runhazard1(im,IM,site,Vs30,opt,model,Nsource,1);
lambdaPSHA = permute(MRE,[2,4,1,3]);
lambdaPSHA = nansum(lambdaPSHA(:,source_ptr),2);


%% plotting
delete(findall(handles.ax2,'Type','line'));
handles.ax2.ColorOrderIndex=1;
if ~isempty(handles.krate)
    plot(handles.ax2,im,lambdaPSHA ,'k--' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    plot(handles.ax2,im,lambdaIS   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    plot(handles.ax2,im,lambdaKM   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    LEG={'PSHA',sprintf('DSHA (IS %g)',length(rate_IS)),sprintf('DSHA (kM %g)',length(krate))};
else
    plot(handles.ax2,im,lambdaPSHA ,'k--' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),    
    plot(handles.ax2,im,lambdaIS   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    LEG={'PSHA',sprintf('DSHA (IS %g)',length(rate_IS))};
end
L=legend(handles.ax2,LEG);
L.Visible='on';
L.Location='northeast';
L.Box='off';

function[]=myfun(hObject, eventdata, handles) %#ok<INUSL>

H=datacursormode(handles.FIGSeismicHazard);
set(H,'enable','on','DisplayStyle','window','UpdateFcn',{@gethazarddata,handles.HazOptions.dbt(1)});
w = findobj('Tag','figpanel');
set(w,'Position',[ 409   485   150    60]);

function output_txt = gethazarddata(~,event_obj,dbt)

pos = get(event_obj,'Position');

if dbt==1
    output_txt = {...
        ['IM   : ',num2str(pos(1),4)],...
        ['Rate : ',num2str(pos(2),4)],...
        ['T    : ',num2str(1/pos(2),4),' years']};
end

if dbt==0
    output_txt = {...
        ['IM             : ',num2str(pos(1),4)],...
        ['P(IM>im|t) : ',num2str(pos(2),4)]};
end
