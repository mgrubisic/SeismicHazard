function[handles]=runPSDAAnalysis(handles)
delete(findall(handles.ax1,'type','line'));
drawnow
d        = handles.d;
Nd       = length(handles.d);
methods  = pshatoolbox_methods(5);
condTerm = handles.ConditionalIntegration.Checked;
label    = {methods.label}';

opt0             = handles.opt;
optL             = handles.opt;
[optL.IM,imethod] = PSDA_im_list(handles);
optL.im  = retrieve_im(opt0.im,opt0.IM,optL.IM);

handles.site_selection = 1:size(handles.h.p,1);
Nsites = length(handles.site_selection);
Nsources=0;
for i=1:length(handles.model)
    Nsources=max(Nsources,length(handles.model(i).source));
end

% Run Hazard Analyses
% -------------------------------------------------------------------------
handles.opt = optL;
if any(ismember(imethod,[1 2]))
    [handles.lambda,handles.deagg]=runlogictree2(handles);
end

if any(ismember(imethod,[3 4]))
    % Rathje R, Rathje F
end

if any(ismember(imethod,5))
    % PCE
end

handles.opt=opt0;

% -------------------------------------------------------------------------

Nbranches = size(handles.table.Data,1);
handles.lambdaD = zeros(Nsites,Nd,Nsources,Nbranches);

fprintf('\n');
spat  = 'Site %-17g | Branch %-3g of %-49g Runtime:  %-4.3f s\n';
t0 = tic;
fprintf('                               SLOPE DISPLACEMENT HAZARD \n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');


for site_ptr=1:Nsites
    for branch_ptr=1:Nbranches
        ti=tic;
        indT1  = handles.IJK(branch_ptr,1); % pointer to scenario and Tm value
        indT2  = handles.IJK(branch_ptr,2); % pointer to Ky and Ts values
        indT3  = handles.IJK(branch_ptr,3); % pointer to analyses models
        
        Ts     = handles.T2{indT2,2};
        ky     = handles.T2{indT2,3};
        B      = zeros(4,1);
        [~,B(1)] = intersect(label,handles.T3{indT3,2}); fun1 = methods(B(1)).func; % interface
        [~,B(2)] = intersect(label,handles.T3{indT3,3}); fun2 = methods(B(2)).func; % slab
        [~,B(3)] = intersect(label,handles.T3{indT3,4}); fun3 = methods(B(3)).func; % crustal
        [~,B(4)] = intersect(label,handles.T3{indT3,5}); fun4 = methods(B(4)).func; % grid, others
        
        ind = 1:Nsources;
        if any(ismember(imethod,[1 2]))
            for kkk=1:Nsources
                if isempty(handles.deagg{site_ptr,1,1,kkk,indT1})
                    ind(kkk)=0;
                end
            end
            ind(ind==0)=[];
        end
        
        % run sources
        for source_ptr=ind
            mechanism = handles.model(indT1).source(source_ptr).mechanism;
            
            switch mechanism
                case 'interface' , fun = fun1; Bs = B(1); param = handles.T3add{indT3,1};
                case 'intraslab' , fun = fun2; Bs = B(2); param = handles.T3add{indT3,2};
                case 'crustal'   , fun = fun3; Bs = B(3); param = handles.T3add{indT3,3};
                case 'grid'      , fun = fun4; Bs = B(4); param = handles.T3add{indT3,4};
            end
            
            switch func2str(fun)
                case 'psda_BMT2017'    , integrator = 1;
                case 'psda_BT2007_1'   , integrator = 1;
                case 'psda_BT2007_2'   , integrator = 2;
                case 'psda_AM1988'     , integrator = 2;
                case 'psda_RathjeR'    , integrator = 3;
                case 'psda_RathjeF'    , integrator = 4;
                case 'psda_BT2007_pce' , integrator = 5;
            end
            
            IMslope         = methods(Bs).Safactor*Ts.*(methods(Bs).Safactor>=0)+methods(Bs).Safactor.*(methods(Bs).Safactor<0);
            [~,period_ptr]  = intersect(optL.IM,IMslope);
            im              = optL.im(:,period_ptr);
            site            = cell2mat(handles.h.p(site_ptr,2:4));
            Vs30            = handles.h.Vs30(site_ptr);
            gmpetype        = handles.model(indT1).source(source_ptr).gmpe.type;
            hd              = nan(1,Nd);
            if integrator==1 && strcmp(gmpetype,'regular')
                deagg       = handles.deagg (site_ptr,:,period_ptr,source_ptr,indT1);  deagg   = permute(deagg ,[2,1]);
                hd          = integrate_Mw_dependent(fun,d,ky,Ts,im,deagg,condTerm);
            end
            
            if integrator==2 && strcmp(gmpetype,'regular')
                lambda      = handles.lambda(site_ptr,:,period_ptr,source_ptr,indT1);  lambda  = permute(lambda,[2,1]);
                hd          = integrate_Mw_independent(fun,d,ky,Ts,im,lambda);
            end
            
            if integrator==3 && strcmp(gmpetype,'regular')
                model          = handles.model(indT1);
                model.source   = model.source(source_ptr);
                [~,~,MRD]      = runhazardV1(im,IMslope,site,Vs30,optL,model,1,param.rho);
                hd             = fun(d,ky,im,MRD);
            end
            
            if integrator==4 && strcmp(gmpetype,'regular')
                model          = handles.model(indT1);
                model.source   = model.source(source_ptr);
                [~,~,MRD]      = runhazardV1(im,IMslope,site,Vs30,optL,model,1,param.rho);
                [Tm,~,dP]      = trlognpdf_psda([param.Tm_mean,param.Tm_cov,param.Tm_Nsta]);
                hd             = zeros(size(d));
                for ix = 1:length(Tm)
                    [im2,MRDkvkm]  = computeMRDkmkv(Ts, Tm(ix), im, MRD,param.rhok);
                    dhd            = fun(d, Ts, ky, im2,MRDkvkm);
                    hd             = hd+dP(ix)*dhd;
                end
            end
            
            if integrator==5 && strcmp(gmpetype,'pce')
                model      = handles.model(indT1);
                source     = model.source(source_ptr);
                ellipsoid  = opt0.ellipsoid;
                Ts_param   = handles.Ts_param;
                ky_param   = handles.ky_param;
                RandType   = 'Stable';
                [PCE,~,~,~,Cz]  = runPCE(source,IMslope,site,Vs30,im,ellipsoid,param.realSa,RandType);
                hd         = fun(d, Ts_param, ky_param, im, PCE,param,Cz,'pce');
                
            end
            
            %             if integrator==5 && strcmp(gmpetype,'pce')
            %                 % Jorge's Validation Example
            %                 Ts_param   = [1.0 0.3 100];
            %                 ky_param   = [0.1 0.3 100];
            %                 load MACEDO2 Sa_vector Hazard_MC_samples PC_Coefficients_Sum_Hazards
            %                 im   = Sa_vector(:);
            %                 PCE  = Hazard_MC_samples';
            %                 Cz   = PC_Coefficients_Sum_Hazards;
            %                 hd   = fun(d, Ts_param, ky_param, im, PCE,param,Cz,'pce');
            %             end
            
            handles.lambdaD(site_ptr,:,source_ptr,branch_ptr) = hd;
        end
        fprintf(spat,site_ptr,branch_ptr,Nbranches,toc(ti))
    end
end

fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

handles.butt2.Value=1;