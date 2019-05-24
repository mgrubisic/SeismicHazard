function[MCSk]=runMCS(source,site,IM,im,ellipsoid,Nsim)


r0      = gps2xyz(site,ellipsoid);
mscl    = source.mscl;
gmpe    = source.gmpe;
Nim     = size(im,1);
rateM   = mscl.dPm;
NMmin   = mscl.msparam.NMmin;

%% ASSEMBLE GMPE PARAMERTER
gmpefun  = gmpe.handle;
imunits  = gmpe.un;
Rmetric  = gmpe.Rmetric;

switch source.type
    case 'point'
        [param,rate] = mGMPEassemblePoint(r0,source,Rmetric,ellipsoid);
    case 'line'
        [param,rate] = mGMPEassembleLine(r0,source,Rmetric,ellipsoid);
    case 'area'
        switch source.mechanism
            case {'interface','intraslab','grid'}
                [param,rate] = mGMPEassembleArea(r0,source,Rmetric,ellipsoid);
            case 'crustal',    [param,rate] = mGMPEassembleCrustal(r0,source,rateM,Rmetric,ellipsoid);
        end
end

%% HAZARD INTEGRAL
Nscen = length(param{1});
zrnd  = normrnd(0,1,1, Nsim);

[lnY,sigma]  = gmpefun(IM,param{:});
mu  = mean(lnY,1);
smu = std (lnY,0,1);

Hazard_F1_MC = zeros(1, Nim, Nreal);
for i = 1:Nscen
    for j = 1:Nreal
        Hazard_F1_MC(1, :, j) = Hazard_F1_MC(1, :, j) + rate(i) * (1 - normcdf((log(imunits*im') - (mu(i) + smu(i) * zrnd(j)))/sigma(i)));
    end
end

MCSk = zeros(Nim, Nsim);
for k = 1:Nsim
    MCSk(:, k) = NMmin*sum(Hazard_F1_MC(:, :, k), 1);
end




