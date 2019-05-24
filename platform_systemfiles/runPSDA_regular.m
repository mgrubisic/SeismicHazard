function[handles]=runPSDA_regular(handles)
delete(findall(handles.ax1,'type','line'));
drawnow
d        = handles.d;
Nd       = length(handles.d);

opt0     = handles.opt;
optL     = handles.opt;
[optL.IM,imethod]  = PSDA_im_list(handles,'model',handles.T2,handles.T3(:,2:5));
optL.im  = retrieve_im(opt0.im,opt0.IM,optL.IM);
handles.site_selection = 1:size(handles.h.p,1);
Nsites = length(handles.site_selection);
Nsources=0;
for i=1:length(handles.model)
    Nsources=max(Nsources,length(handles.model(i).source));
end

% Pre-run Hazard Analyses
% -------------------------------------------------------------------------
handles.opt = optL;
if any(ismember(imethod,[1 2]))
    [handles.lambda,handles.deagg]=runlogictree2(handles);
end
handles.opt=opt0;
% -------------------------------------------------------------------------

Nbranches       = size(handles.table.Data,1);
handles.lambdaD = zeros(Nsites,Nd,Nsources,Nbranches);

fprintf('\n');
spat  = 'Site %-17g | Branch %-3g of %-49g Runtime:  %-4.3f s\n';
t0 = tic;
fprintf('                               SLOPE DISPLACEMENT HAZARD \n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');

SMLIB = handles.sys.SMLIB;
id={handles.sys.SMLIB.id};
for site_ptr=1:Nsites
    for branch_ptr=1:Nbranches
        ti=tic;
        indT1    = handles.IJK(branch_ptr,1); % pointer to scenario and Tm value
        indT2    = handles.IJK(branch_ptr,2); % pointer to Ky and Ts values
        indT3    = handles.IJK(branch_ptr,3); % pointer to analyses models
        
        Ts       = handles.T2{indT2,2};
        ky       = handles.T2{indT2,3};
        B        = zeros(4,1);
        
        [~,B(1)] = intersect(id,handles.T3{indT3,2}); fun1 = SMLIB(B(1)).func; % interface
        [~,B(2)] = intersect(id,handles.T3{indT3,3}); fun2 = SMLIB(B(2)).func; % slab
        [~,B(3)] = intersect(id,handles.T3{indT3,4}); fun3 = SMLIB(B(3)).func; % crustal
        [~,B(4)] = intersect(id,handles.T3{indT3,5}); fun4 = SMLIB(B(4)).func; % grid, others
        
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
                case 'interface' , fun = fun1; Bs=B(1);
                case 'intraslab' , fun = fun2; Bs=B(2);
                case 'crustal'   , fun = fun3; Bs=B(3);
                case 'grid'      , fun = fun4; Bs=B(4);
            end
            
            param      = SMLIB(Bs).param;
            integrator = SMLIB(Bs).integrator;
            Safactor   = SMLIB(Bs).Safactor;
            IMslope    = Safactor*Ts.*(Safactor>=0)+Safactor.*(Safactor<0);
            
            [~,period_ptr]  = intersect(optL.IM,IMslope);
            im              = optL.im(:,period_ptr);
            
            if integrator==1 % magnitude dependent models
                deagg    = handles.deagg (site_ptr,:,period_ptr,source_ptr,indT1);
                deagg    = permute(deagg ,[2,1]);
                condTerm = handles.ConditionalIntegration.Checked;
                hd       = integrate_Mw_dependent(fun,d,ky,Ts,im,deagg,condTerm);
            end
            
            if integrator==2 % magnitude independent models
                lambda  = handles.lambda(site_ptr,:,period_ptr,source_ptr,indT1);  
                lambda  = permute(lambda,[2,1]);
                hd      = integrate_Mw_independent(fun,d,ky,Ts,im,lambda);
            end
            
            if integrator==3 % Ellen´s rigid slope
                site          = cell2mat(handles.h.p(site_ptr,2:4));
                Vs30          = handles.h.Vs30(site_ptr);
                model         = handles.model(indT1);
                model.source  = model.source(source_ptr);
                [~,~,MRD]     = runhazardV1(im,IMslope,site,Vs30,optL,model,1,param.rho);
                hd            = fun(d,ky,im,MRD);
            end
            
            if integrator==4 % Ellen´s flexible slope
                site          = cell2mat(handles.h.p(site_ptr,2:4));
                Vs30          = handles.h.Vs30(site_ptr);
                model         = handles.model(indT1);
                model.source  = model.source(source_ptr);
                [~,~,MRD]     = runhazardV1(im,IMslope,site,Vs30,optL,model,1,param.rho);
                [Tm,~,dP]     = trlognpdf_psda([param.Tm_mean,param.Tm_cov,param.Tm_Nsta]);
                hd            = zeros(size(d));
                for ix = 1:length(Tm)
                    [im2,MRDkvkm]  = computeMRDkmkv(Ts, Tm(ix), im, MRD,param.rhok);
                    dhd            = fun(d, Ts, ky, im2,MRDkvkm);
                    hd             = hd+dP(ix)*dhd;
                end
            end
            
           
            handles.lambdaD(site_ptr,:,source_ptr,branch_ptr) = hd;
        end
        fprintf(spat,site_ptr,branch_ptr,Nbranches,toc(ti))
    end
end

fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

handles.butt2.Value=1;
