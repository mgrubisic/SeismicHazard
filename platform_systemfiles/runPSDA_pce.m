function[handles]=runPSDA_pce(handles)
delete(findall(handles.ax1,'type','line'));
drawnow
d          = handles.d;
realSa     = handles.realSa;
realD      = handles.realD;
Nscenarios = handles.realSa*handles.realD;
Nd         = length(handles.d);
opt0       = handles.opt;
optL       = handles.opt;
T2         = handles.tablePCE.Data(:,1:2);
T3         = handles.tablePCE.Data(:,[4,5,6,7]);
optL.IM    = PSDA_im_list(handles,'modelpce',T2,T3);
optL.im    = retrieve_im(opt0.im,opt0.IM,optL.IM);
handles.site_selection = 1:size(handles.h.p,1);
Nsites = length(handles.site_selection);
Nsources=0;
for i=1:length(handles.modelpce)
    Nsources=max(Nsources,length(handles.modelpce(i).source));
end

Nmodels           = size(handles.tablePCE.Data,1);
handles.lambdaPCE = zeros(Nscenarios,Nsites,Nd,Nsources,Nmodels);

fprintf('\n');
spat  = 'Site %-17g | Branch %-3g of %-49g Runtime:  %-4.3f s\n';
t0 = tic;
fprintf('                               SLOPE DISPLACEMENT HAZARD \n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');

SMLIB = handles.sys.SMLIB;
id    = {handles.sys.SMLIB.id};
TPCE  = handles.tablePCE.Data;

ellipsoid = opt0.ellipsoid;

for site_ptr=1:Nsites
    for model_ptr=1:Nmodels
        site  = cell2mat(handles.h.p(site_ptr,2:4));
        xyz   = gps2xyz(site,ellipsoid);
        Vs30  = handles.h.Vs30(site_ptr);
        ti=tic;
        
        Ts_param = str2double(regexp(TPCE{model_ptr,2},'\, ','split'));
        ky_param = str2double(regexp(TPCE{model_ptr,3},'\, ','split'));
        B        = zeros(4,1);
        
        [~,B(1)] = intersect(id,TPCE{model_ptr,4}); fun1 = SMLIB(B(1)).func; % interface
        [~,B(2)] = intersect(id,TPCE{model_ptr,5}); fun2 = SMLIB(B(2)).func; % slab
        [~,B(3)] = intersect(id,TPCE{model_ptr,6}); fun3 = SMLIB(B(3)).func; % crustal
        [~,B(4)] = intersect(id,TPCE{model_ptr,7}); fun4 = SMLIB(B(4)).func; % grid, others
        
        
        [~,Bmodel]=intersect({handles.modelpce.id},handles.tablePCE.Data{model_ptr,1});
        modelpce = handles.modelpce(Bmodel);
        ind=selectsource(opt0.MaxDistance,xyz,modelpce.source,ellipsoid);
        ind = find(ind);
        
        % run sources
        for source_ptr=ind
            source    = modelpce.source(source_ptr);
            
            switch source.mechanism
                case 'interface' , fun = fun1; Bs=B(1);
                case 'intraslab' , fun = fun2; Bs=B(2);
                case 'crustal'   , fun = fun3; Bs=B(3);
                case 'grid'      , fun = fun4; Bs=B(4);
            end
            
            param      = SMLIB(Bs).param;
            integrator = SMLIB(Bs).integrator;
            Safactor   = SMLIB(Bs).Safactor;
            Ts         = Ts_param(1);
            IMslope    = Safactor*Ts.*(Safactor>=0)+Safactor.*(Safactor<0);
            
            [~,period_ptr]  = intersect(optL.IM,IMslope);
            im              = optL.im(:,period_ptr);
            
            if integrator==5
                sampling   = param.sampling;
                [PCE,~,~,~,Cz]  = runPCE(source,IMslope,site,Vs30,im,ellipsoid,realSa,sampling);
                lambda  = fun(d, Ts_param, ky_param, im, PCE,param,Cz,realSa,realD); %psda_BT2007_pce
            end

            handles.lambdaPCE(:,site_ptr,:,source_ptr,model_ptr) = lambda;
        end
        fprintf(spat,site_ptr,model_ptr,Nmodels,toc(ti))
    end
end

fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

