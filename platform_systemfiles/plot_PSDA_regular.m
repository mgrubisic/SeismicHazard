function plot_PSDA_regular(handles)

if ~isfield(handles,'sys'), return;end
if ~isfield(handles,'lambdaD'),return;end

delete(findall(handles.ax1,'type','line'));
handles.ax1.ColorOrderIndex=1;
handles.ax1.NextPlot='add';

site_ptr = handles.pop_site.Value;

REGOptions = handles.REGOptions;
d = handles.d;

modo = REGOptions.mode;


switch modo
    case 1 % averaged displacment hazard
        lambdaD = sum(handles.lambdaD(site_ptr,:,:,:),3);
        lambdaD = permute(lambdaD,[4 2 3 1]);
        
        if REGOptions.param(1)==1 % default weights
            handles.table.Data = main_psda(handles.T1,handles.T2,handles.T3);
            weights = cell2mat(handles.table.Data(:,end));
            lam2 = prod(bsxfun(@power,lambdaD,weights),1);
            plot(handles.ax1,d,lam2,'b.-');
            str     = 'Default Weights';
        end
        
        if REGOptions.param(2)==1 % random weights
            weights = handles.REGOptions.rnd;
            %Nt      = length(weights);
            handles.table.Data(:,4)=num2cell(weights);
            lam2 = prod(bsxfun(@power,lambdaD,weights),1);
            plot(handles.ax1,d,lam2,'b.-');
            str     = 'Random Weights';
        end
        
        if REGOptions.param(3)==1 % percentiles
            percentile = REGOptions.param(4);
            Nt  = size(lambdaD,1);
            handles.table.Data(:,4)=cell(Nt,1);
            lam2 = exp(prctile(log(lambdaD),percentile,1));
            plot(handles.ax1,d,lam2,'b.-');
            str  = sprintf('Percentile %g',percentile);
        end
        
        if REGOptions.param(5)==1 % displays source contribution
            Ny = size(lambdaD,1);
            c1 = [0.4660    0.6740    0.1880];
            c2 = [0.3010    0.7450    0.9330];
            mycolor = [linspace(c1(1),c2(1),Ny)',linspace(c1(2),c2(2),Ny)',linspace(c1(3),c2(3),Ny)'];
            handles.ax1.ColorOrder=[0 0.477 0.741; mycolor];
            lam1 = lambdaD;
            plot(handles.ax1,d,lam1,'-')
        end
        
        if REGOptions.param(5)==0 % legend
            Leg=legend(handles.ax1,str);
        else
            Nbranches = size(handles.table.Data,1);
            strb = cell(1,1+Nbranches);
            strb{1}=str;
            for i=1:Nbranches
                strb{i+1}=sprintf('Branch %g',i);
            end
            Leg=legend(handles.ax1,strb);
            Leg.FontSize=8;
            Leg.EdgeColor=[1 1 1];
            Leg.Location='SouthWest';
            Leg.Tag='hazardlegend';
        end
        
        switch handles.toggle1.Value
            case 0,Leg.Visible='off';
            case 1,Leg.Visible='on';
        end
        
        
    case 2 % single branch displacment
        branch_ptr = REGOptions.param(6);
        lambdaD = handles.lambdaD(site_ptr,:,:,branch_ptr);
        lambdaD = permute(lambdaD,[2 3 1]);
        
        lam2 = sum(lambdaD,2)';
        plot(handles.ax1,d,lam2,'b.-');
        str = sprintf('Branch %g',branch_ptr);
        
        if REGOptions.param(7)==1 % displays source contribution
            haz_ptr = handles.IJK(branch_ptr,1);
            NOTZERO = sum(lambdaD,1)>0;
            source_label = {handles.model(haz_ptr).source.label};
            str  = [str,source_label(NOTZERO)];
            
            Ny = sum(NOTZERO);
            c1 = [0.4660    0.6740    0.1880];
            c2 = [0.3010    0.7450    0.9330];
            mycolor = [linspace(c1(1),c2(1),Ny)',linspace(c1(2),c2(2),Ny)',linspace(c1(3),c2(3),Ny)'];
            handles.ax1.ColorOrder=[0 0.477 0.741; mycolor];
            lam1 = lambdaD(:,NOTZERO)';
            plot(handles.ax1,d,lam1,'-')
        end
        
        if REGOptions.param(8)==1 % displays source contribution
            branch_ptr = REGOptions.param(6);
            model_ptr = handles.IJK(branch_ptr,1);
            mechs     = {handles.model(model_ptr).source.mechanism};
            m1        = strcmp(mechs,'system');
            m2        = strcmp(mechs,'interface');
            m3        = strcmp(mechs,'intraslab');
            m4        = strcmp(mechs,'slab');
            m5        = strcmp(mechs,'crustal');
            m6        = strcmp(mechs,'fault');
            m7        = strcmp(mechs,'grid');
            lambdaD   = [nansum(lambdaD(:,m1),2) nansum(lambdaD(:,m2),2) nansum(lambdaD(:,m3),2) nansum(lambdaD(:,m4),2) nansum(lambdaD(:,m5),2) nansum(lambdaD(:,m6),2) nansum(lambdaD(:,m7),2)];
            NOTNAN    = (sum(lambdaD,1)>0);
            lam1      = lambdaD(:,NOTNAN)';
            plot(handles.ax1,d,lam1,'-')
            mechs = {'system','interface','intraslab','slab','crustal','fault','grid'};
            str = [{str},mechs(NOTNAN)];
        end
        
        Leg=legend(handles.ax1,str);
        Leg.FontSize=8;
        Leg.EdgeColor=[1 1 1];
        Leg.Location='SouthWest';
        Leg.Tag='hazardlegend';
        switch handles.toggle1.Value
            case 0,Leg.Visible='off';
            case 1,Leg.Visible='on';
        end
end

axis(handles.ax1,'auto')
cF   = get(0,'format');
format long g
if exist('lam1','var')
    data = num2cell([d;lam2;lam1]'); % branches
else
    data = num2cell([d;lam2]'); % average
end
c    = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',        {@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback',           {@figure2clipboard_uimenu,handles.ax1});
uimenu(c,'Label','Undock & compare','Callback', {@figurecompare_uimenu,handles.ax1});
set(handles.ax1,'uicontextmenu',c);
format(cF);

