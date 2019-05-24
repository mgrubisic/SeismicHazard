function[deagg]=runhazard3(im,IM,site,Vs30,opt,model,Nsource,site_selection)

ellipsoid = opt.ellipsoid;
xyz       = gps2xyz(site,ellipsoid);
Nsite     = size(xyz,1);
NIM       = length(IM);
Nim       = size(im,1);
deagg     = cell(Nsite,Nim,NIM,Nsource);
epsilon   = -4.1:0.2:4.1;
emin      = epsilon (1);
emax      = epsilon (end);
deps      = mean(diff(epsilon));

ind  = zeros(Nsite,length(model.source));
for i=site_selection
    ind(i,:)=selectsource(opt.MaxDistance,xyz(i,:),model.source,ellipsoid);
end

for k=site_selection
    ind_k      = ind(k,:);
    sptr       = find(ind_k);
    xyzk       = xyz(k,:);
    Vs30k      = Vs30(k);
    source     = model.source(ind_k==1);
    Nsource_k  = length(source);
    for i=1:Nsource_k
        ind_s = sptr(i);
        source_i = mGMPEVs30(source(i),Vs30k);
        deagg(k,:,:,ind_s)=runsourceDeagg(source_i,xyzk,IM,im,ellipsoid,emin,emax,deps);
    end
end


return

function[deagg]=runsourceDeagg(source,r0,T,im,ellipsoid,emin,emax,deps)

mscl = source.mscl;
gmpe = source.gmpe;

%% MAGNITUDE RATE OF EARTHQUAKES
Nper        = length(T);
Nim         = size(im,1);
rateM       = mscl.dPm;
NMmin       = mscl.msparam.NMmin;

%% ASSEMBLE GMPE PARAMERSER
gmpefun  = gmpe.handle;
imunits  = gmpe.un;
Rmetric  = gmpe.Rmetric;

switch source.type
    case 'point'
        [param,rate_MR] = mGMPEassemblePoint(r0,source,Rmetric,ellipsoid);
    case 'line'
        [param,rate_MR] = mGMPEassembleLine(r0,source,Rmetric,ellipsoid);
    case 'area'
        switch source.mechanism
            case {'interface','intraslab','grid'}
                [param,rate_MR] = mGMPEassembleArea(r0,source,Rmetric,ellipsoid);
            case 'crustal'
                [param,rate_MR] = mGMPEassembleCrustal(r0,source,rateM,Rmetric,ellipsoid);
        end
end

%% HAZARD INTEGRAL
sigma     = gmpe.usp.sigma;
std_exp   = 1;
sig_overw = 1;

if ~isempty(sigma)
    switch sigma{1}
        case 'overwrite', std_exp = 0; sig_overw = sigma{2};
        case 'truncate' , emin = min([emin,sigma{2},emin+sigma{2}]); emax = sigma{2};
    end
end

deagg    = cell(Nim,Nper);
ep1      = emin:deps:(emax-deps);
ep2      = (emin+deps):deps:emax;
epsilon  = 1/2*(ep1+ep2);
Neps     = length(epsilon);

ep1(1)   = -inf;
ep2(end) = inf;

rate_e   = normcdf(ep2)-normcdf(ep1);

switch gmpe.type
    case 'regular'
    NMR      = size(param{1},1);
    Mag      = repmat(param{1},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{2},Neps,1);%creating copies for the distance deagregation
    case 'udm'
    NMR      = size(param{5},1);
    Mag      = repmat(param{5},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{6},Neps,1);%creating copies for the distance deagregation
    case 'cond'
    NMR      = size(param{5},1);
    Mag      = repmat(param{5},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{6},Neps,1);%creating copies for the distance deagregation        
end

epsilon  = repmat(epsilon',NMR,1); epsilon  = epsilon(:);
rate_e   = repmat(rate_e',NMR,1);  rate_e   = rate_e(:);
rate_MR  = repmat(rate_MR,Neps,1);
rate     = rate_MR.*rate_e;



for j=1:Nper
    [mu,std] = gmpefun(T(j),param{:});
    std      = std.^std_exp*sig_overw;
    mu       = repmat(mu,Neps,1);
    std      = repmat(std,Neps,1);
    lnSa     = mu+std.*epsilon;
    imj = im(:,j); 
    
    for i=1:Nim
        im_i = log(imunits*imj(i));
        ccdf = (lnSa>=im_i);
        deagg{i,j}  = [Mag,Rup,epsilon,NMmin*ccdf.*rate];%that will provide the rates for all possible combinations of M, D, and epsilon
    end
end

return

