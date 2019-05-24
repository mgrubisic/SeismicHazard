function[handles]=mat2psda(handles,varargin)

%% loads sys, opt, h, model
if nargin==3
    pathname=varargin{1};
    filename=varargin{2};
    if contains(filename,'.mat')
        load([pathname,filename],'sys','opt','h')
    elseif contains(filename,'.txt')
        [sys,opt,h]=loadPSHA([pathname,filename]);
    end
elseif nargin==5
    sys    = varargin{1};
    opt    = varargin{3};
    h      = varargin{4};
end

opt.MagDiscrete  = {'uniform',0.1};
models           = process_model(sys,opt);
isREGULAR        = horzcat(models.isregular);
handles.sys      = sys;
handles.model    = models(isREGULAR==1);
handles.modelpce = models(isREGULAR==0);
handles.opt      = opt;
handles.h        = h;
ME               = pshatoolbox_methods(5);

%% loads default models
func   = {ME.str};
ptrs   = sys.PTRS;

[d_default,Sadef,Ddef,SMLIB_default,default_reg,default_pce]=loadPSDADefaultParam(handles.model,handles.modelpce);

% d values
if ~isnan(ptrs(10,1))
    str = sys.DATA(ptrs(10,1):ptrs(10,2),:);
    str = regexp(str,'\ : ','split');
    handles.d      = eval(str{1}{2});
    handles.realSa = str2double(str{2}{2});
    handles.realD  = str2double(str{3}{2});
else
    handles.d      = d_default;
    handles.realSa = Sadef;
    handles.realD  = Ddef;
end

% ky and Ts
if ~isnan(ptrs(11,1))
    str = sys.DATA(ptrs(11,1):ptrs(11,2),:);
    str = regexp(str,'\ : ','split');
    
    % param 1
    fld = [str{1}{1},'_param'];
    val = lower(regexp(str{1}{2},'\ ','split'));
    val = struct(val{:});
    handles.(fld) = [str2double(val.mean),str2double(val.cov),str2double(val.samples)];
    
    %param2
    fld = [str{2}{1},'_param'];
    val = lower(regexp(str{2}{2},'\ ','split'));
    val = struct(val{:});
    handles.(fld) = [str2double(val.mean),str2double(val.cov),str2double(val.samples)];
else
    handles.Ts_param = default_reg.Ts_param;
    handles.ky_param = default_reg.ky_param;
end

% displacement model library
if ~isnan(ptrs(12,1))
    reservednames = {'avgHazard','integration','sampling'};
    str = sys.DATA(ptrs(12,1):ptrs(12,2),:);
    str = regexp(str,'\ handle ','split');
    
    Nmodels = length(str);
    SMLIB(1:Nmodels,1) = struct('id',[],'label',[],'str',[],'func',[],'isregular',[],'integrator',[],'primaryIM',[],'Safactor',[],'param',[]);
    for i=1:Nmodels
        stri = regexp(str{i}{2},'\ ','split');
        [~,C]=intersect(func,stri{1});
        SMLIB(i).id = str{i}{1};
        SMLIB(i).str        = ME(C).str;
        SMLIB(i).func       = ME(C).func;
        SMLIB(i).isregular  = ME(C).isregular;
        SMLIB(i).integrator = ME(C).integrator;
        SMLIB(i).primaryIM  = ME(C).primaryIM;
        SMLIB(i).Safactor   = ME(C).Safactor;
        if length(stri)>2
            fixparam = struct(stri{2:end});
            flds     = fields(fixparam);
            for j=1:length(flds)
                if ~ismember(flds{j},reservednames)
                    fixparam.(flds{j})=str2double(fixparam.(flds{j}));
                end
            end
            SMLIB(i).param   = fixparam;
        end
    end
    handles.sys.SMLIB=SMLIB;
else
    handles.sys.SMLIB = SMLIB_default;
end

% Displacement Models for regular PSDA Analysis
if ~isnan(ptrs(13,1))
    str = sys.DATA(ptrs(13,1):ptrs(13,2),:);
    Nmodels = size(str,1)-1;
    T3    = cell(Nmodels,6);
    newline = regexp(str{1},'\ ','split');
    slopeweights = str2double(newline(3:end)');
    slopeweights = slopeweights / sum(slopeweights);
    T3(:,6) = num2cell(slopeweights);
    
    str(1,:)=[];
    
    for j=1:Nmodels
        strj = regexp(str{j},'\ ','split');
        T3{j,1}=strj{1};
        modelassig = struct(strj{2:9});
        T3{j,2} = modelassig.interface;
        T3{j,3} = modelassig.intraslab;
        T3{j,4} = modelassig.crustal;
        T3{j,5} = modelassig.grid;
    end
    handles.T3=T3;
else
    handles.T3=default_reg.T3;
end

%% Displacement Models for PCE PSDA Analysis
if ~isnan(ptrs(14,1))
    str = sys.DATA(ptrs(14,1):ptrs(14,2),:);
    str = regexp(str,'\ ','split');
    
    Nmodels  = size(str,1);
    DataPCE  = cell(Nmodels,7);
    
    notPCE  = zeros(Nmodels,1);
    for j=1:Nmodels
        [~,Bj]=intersect(sys.BRANCH,str2double(str{j}(3:5)),'rows');
        modelassig   = struct(str{j}{6:end});
        model_ID     = models(Bj).id;
        if contains(model_ID,'(PCE)')
            DataPCE{j,1} = model_ID;
            DataPCE{j,2} = sprintf('%s, %s',modelassig.meanTs,modelassig.covTs);
            DataPCE{j,3} = sprintf('%s, %s',modelassig.meanky,modelassig.covky);
            DataPCE{j,4} = modelassig.interface;
            DataPCE{j,5} = modelassig.intraslab;
            DataPCE{j,6} = modelassig.crustal;
            DataPCE{j,7} = modelassig.grid;
        else
            fprintf('Warning: Hazard from model ID: "%s" is not PCE compatible\n',model_ID)
            notPCE(j)=true;
        end
    end
    DataPCE(notPCE==1,:)=[];

    isPCEGMM = ~horzcat(SMLIB.isregular);
    handles.tablePCE.ColumnFormat{4}={SMLIB(isPCEGMM).id};
    handles.tablePCE.ColumnFormat{5}={SMLIB(isPCEGMM).id};
    handles.tablePCE.ColumnFormat{6}={SMLIB(isPCEGMM).id};
    handles.tablePCE.ColumnFormat{7}={SMLIB(isPCEGMM).id};
    handles.tablePCE.Data = DataPCE;

else
    if ~isempty(handles.modelpce)
        handles.tablePCE.Data = default_pce;
    end
end

%% setup GUI for regular models
if any(isREGULAR)
    handles.pop_site.String=handles.h.p(:,1);
    handles.pop_site.Enable='on';
    handles.pop_site.Value=1;
    
    % Tables T1,T2
    pshaweights = sys.WEIGHT(isREGULAR==1,4);
    w1          = pshaweights/sum(pshaweights);
    id          = {handles.model.id}';
    handles.T1  = [id,num2cell(w1)];
    [Ts,~,dPTs] = trlognpdf_psda(handles.Ts_param);
    [ky,~,dPky] = trlognpdf_psda(handles.ky_param);
    Ts          = round(Ts*1e10)/1e10;
    [ind1,ind2] = meshgrid(1:length(Ts),1:length(ky));
    ind1        = ind1(:);
    ind2        = ind2(:);
    Nparam      = length(ind1);
    param_id    = cell(Nparam,1);
    for i=1:Nparam
        param_id{i}=sprintf('set%g',i);
    end
    handles.T2  = [param_id,num2cell([Ts(ind1),ky(ind2),dPTs(ind1).*dPky(ind2)])];
    
    [handles.table.Data,handles.IJK]=main_psda(handles.T1,handles.T2,handles.T3);
    handles.EditLogicTree.Enable='on';
end

%% two models
if ~isempty(handles.model) && ~isempty(handles.modelpce)
    handles.mode1.Value       = 1;
    handles.table.Enable      = 'on';
    handles.runRegular.Enable = 'on';
    handles.treebutton.Enable = 'on';
    handles.REG_DisplayOptions.Enable = 'on';
    handles.tablePCE.Enable = 'off';
    handles.runPCE.Enable   = 'inactive';
    handles.PCE_DisplayOptions.Enable='inactive';
end

%% regular yes, pce no
if ~isempty(handles.model) && isempty(handles.modelpce)
    handles.mode1.Value       = 1;
    handles.table.Enable      = 'on';
    handles.runRegular.Enable = 'on';
    handles.treebutton.Enable = 'on';
    handles.REG_DisplayOptions.Enable = 'on';
    handles.tablePCE.Enable = 'off';
    handles.runPCE.Enable   = 'inactive';
    handles.PCE_DisplayOptions.Enable='inactive';
end

%% regular yes, pce no
if isempty(handles.model) && ~isempty(handles.modelpce)
    handles.mode2.Value       = 1;
    handles.table.Enable      = 'off';
    handles.runRegular.Enable = 'inactive';
    handles.treebutton.Enable = 'inactive';
    handles.REG_DisplayOptions.Enable = 'inactive';
    handles.tablePCE.Enable = 'on';
    handles.runPCE.Enable   = 'on';
    handles.PCE_DisplayOptions.Enable='on';
end

