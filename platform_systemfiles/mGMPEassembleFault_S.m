function [param,rate,iR]=mGMPEassembleFault_S(r0,source,rateM,Rmetric,ellipsoid) %#ok<STOUT>

geom = source.geom;
mscl = source.mscl;
gmpe = source.gmpe;
rupt = source.rupt;

%% RUPTURE AREA AND SCENARIOS
M         = mscl.M;
rupArea   = rupRelation(M,0,geom.Area,rupt.RA{1},rupt.RA{2});
aratio    = rupt.aratio;
rupLength = sqrt(rupArea*aratio);
rupWidth  = sqrt(rupArea/aratio);

nM    = 1;
if Rmetric(1),  rrup  = cell(1,nM);end
if Rmetric(2),  rhyp  = cell(1,nM);end
if Rmetric(3),  rjb   = cell(1,nM);end
if Rmetric(4),  repi  = cell(1,nM);end
if Rmetric(5),  rseis = cell(1,nM);end
if Rmetric(6),  rx    = cell(1,nM);end
if Rmetric(7),  ry0   = cell(1,nM);end
if Rmetric(8),  zhyp  = cell(1,nM);end
if Rmetric(9),  ztor  = cell(1,nM);end
if Rmetric(10), zbor  = cell(1,nM);end
if Rmetric(11), zbot  = cell(1,nM);end

p     = geom.p;
pmean = geom.pmean;
rot   = geom.rot;

for i=1:nM
    RA     = rupRelation(M,0,Inf,rupt.RA{1},rupt.RA{2});
    aratio = rupt.aratio;
    RL     = sqrt(RA*aratio);
    RW     = sqrt(RA/aratio);
    rf     = geom.hypm;
    
    if Rmetric(1),  rrup {i} = dist_rrup4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(2),  rhyp {i} = dist_rhyp4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(3),  rjb  {i} = dist_rjb4   (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(4),  repi {i} = dist_repi4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(5),  rseis{i} = dist_rseis4 (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(6),  rx   {i} = dist_rx4    (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(7),  ry0  {i} = dist_ry04   (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(8),  zhyp {i} = dist_zhyp4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(9),  ztor {i} = dist_ztor4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(10), zbor {i} = dist_zbor4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(11), zbot {i} = dist_zbot4  (r0,rf,RW,RL,geom,ellipsoid);end
    
end

if Rmetric(1)==1
    Mw    = cell(size(rrup));
    rate  = cell(size(rrup));
    for i=1:nM
        nri =size(rrup{i},1);
        Mw{i}=M(i)*ones(nri,1);
        rate{i}=1/nri*ones(nri,1)*rateM(i);
    end
elseif Rmetric(3)==1 && Rmetric(1)==0
    Mw    = cell(size(rjb));
    rate  = cell(size(rjb));
    for i=1:nM
        nri =size(rjb{i},1);
        Mw{i}=M(i)*ones(nri,1);
        rate{i}=1/nri*ones(nri,1)*rateM(i);
    end
end

M     = vertcat(Mw{:});
rate  = vertcat(rate{:});
if Rmetric(1),  rrup  = vertcat(rrup {:});end
if Rmetric(2),  rhyp  = vertcat(rhyp {:});end
if Rmetric(3),  rjb   = vertcat(rjb  {:});end
if Rmetric(4),  repi  = vertcat(repi {:});end  %#ok<*NASGU> % not used by any GMPE listed below
if Rmetric(5),  rseis = vertcat(rseis{:});end  % appears in Cambell, but is that GMPE is not implemented
if Rmetric(6),  rx    = vertcat(rx   {:});end
if Rmetric(7),  ry0   = vertcat(ry0  {:});end
if Rmetric(8),  zhyp  = vertcat(zhyp {:});end
if Rmetric(9),  ztor  = vertcat(ztor {:});end
if Rmetric(10), zbor  = vertcat(zbor {:});end  % used in NGA West2, but not necesary ig W is given
if Rmetric(11), zbot  = vertcat(zbot {:});end  % used in NGA West2, but not necesary ig W is given
usp   = gmpe.usp;

%% GMM PARAMETERS

switch func2str(gmpe.handle)
    % subduction sources
    case 'Arteta2018'
        param   = {M,rhyp,usp.media,usp.region};
    case 'Youngs1997'
        param   = {M,rrup,zhyp,usp.mechanism,usp.media};
    case 'AtkinsonBoore2003'
        param   = {M,rrup,zhyp,usp.mechanism,usp.media,usp.region};
    case 'Zhao2006'
        param   = {M,rrup,zhyp,usp.mechanism,usp.Vs30};
    case 'Mcverry2006'
        param   = {M,rrup,zhyp,usp.mechanism,usp.media,usp.rvol};
    case 'BCHydro2012'
        switch usp.mechanism
            case 'interface'
                param     = {M,rrup,zhyp,usp.mechanism,usp.region,usp.DeltaC1,usp.Vs30};
            case 'intraslab'
                param     = {M,rhyp,zhyp,usp.mechanism,usp.region,usp.DeltaC1,usp.Vs30};
        end
    case 'MontalvaBastias2017'
        switch usp.mechanism
            case 'interface'
                param   = {M,rrup,zhyp,usp.mechanism,usp.region,usp.Vs30};
            case 'intraslab'
                param   = {M,rhyp,zhyp,usp.mechanism,usp.region,usp.Vs30};
        end
    case 'SiberRisk2019'
        switch usp.mechanism
            case 'interface'
                param   = {M,rrup,zhyp,usp.mechanism,usp.Vs30};
            case 'intraslab'
                param   = {M,rhyp,zhyp,usp.mechanism,usp.Vs30};
        end
    case 'Idini2016'
        switch usp.mechanism
            case 'interface'
                R=rrup; R(M<7.7)=rhyp(M<7.7);
                param   = {M,R,zhyp,usp.mechanism,usp.spectrum,usp.Vs30};
            case 'intraslab'
                param   = {M,rhyp,zhyp,usp.mechanism,usp.spectrum,usp.Vs30};
        end
    case 'ContrerasBoroschek2012'
        param     = {M,rrup,zhyp,usp.media};
    case 'Jaimes2006'
        param   = {M,rrup};
    case 'Jaimes2015'
        param   = {M,rrup};
    case 'Jaimes2016'
        param   = {M,rrup};
    case 'GarciaJaimes2017'
        param   = {M,rrup};
    case 'Sadigh1997'
        param   = {M,rrup,usp.mechanism,usp.media};
    case 'AS1997h'
        usp.location = 'footwall';
        param   = {M,rrup,usp.media,usp.mechanism,usp.location,usp.sig};
    case 'Idriss2008_nga'
        param     = {M,rrup,usp.mechanism,usp.Vs30};
    case 'ChiouYoungs2008_nga'
        dip     = geom.dip;
        z1      = exp(28.5-3.82/8*log(usp.Vs30^8+378.8^8));
        param   = {M, rrup, rjb, rx, ztor, dip, z1, usp.mechanism, usp.event, usp.Vs30, usp.Vs30type};
    case 'BooreAtkinson_2008_nga'
        param   = {M,rjb,usp.mechanism,usp.Vs30};
    case 'CampbellBozorgnia_2008_nga'
        dip     = geom.dip;
        param   = {M,rrup, rjb, ztor, dip, usp.mechanism, usp.Vs30, usp.Z25, usp.sigmatype};
    case 'AbrahamsonSilva2008_nga'
        dip     = geom.dip;
        dipvec  = [1 -1]*gps2xyz(source.vertices(1:2,:),ellipsoid);
        W       = norm(dipvec);
        Z10     = Z10_default_AS08_NGA(usp.Vs30);
        param   = {M, rrup, rjb, rx, ztor, dip,W, Z10, usp.Vs30,usp.mechanism,usp.event,usp.Vs30type};
    case 'I_2014_nga'
        param     = {M,rrup,usp.mechanism,usp.Vs30};
    case 'CY_2014_nga'
        dip     = geom.dip;
        param   = {M, rrup, rjb, rx, ztor, dip, usp.mechanism, usp.Z10, usp.Vs30, usp.Vs30type, usp.region};
    case 'BSSA_2014_nga'
        param   = {M, rjb, usp.mechanism, usp.region, usp.BasinDepth, usp.Vs30};
    case 'CB_2014_nga'
        dip     = geom.dip;
        dipvec  = [1 -1]*gps2xyz(source.vertices(1:2,:),ellipsoid);
        W       = norm(dipvec);
        param   = {M, rrup, rjb, rx, W, ztor, zhyp, dip, usp.mechanism, usp.HW, usp.Vs30, usp.Z25, usp.region};
    case 'ASK_2014_nga'
        dip     = geom.dip;
        dipvec  = [1 -1]*gps2xyz(source.vertices(1:2,:),ellipsoid);
        W       = norm(dipvec);
        param={M, rrup, rjb, rx, ry0, ztor, dip, W, usp.mechanism, usp.event, usp.Z10, usp.Vs30, usp.Vs30type, usp.region};
    case 'udm'
        var  = gmpe.var;
        txt  = regexp(var.syntax,'\(','split');
        args = regexp(txt{2}(1:end-1),'\,','split');
        args = strtrim(args);
        args(1)  = [];
        param    = cell(1,4+length(args));
        param{1} = str2func(strtrim(txt{1}));
        param{2} = var.vector;
        param{3} = var.units;
        param{4} = var.residuals;
        for cont=1:length(args)
            f = var.(args{cont});
            if strcmpi(f.tag,'magnitude')
                param{4+cont}=M;
            end
            
            if strcmpi(f.tag,'distance')
                fval = find(f.value);
                switch fval
                    case 1 , param{4+cont}=rrup;
                    case 2 , param{4+cont}=rhyp;
                    case 3 , param{4+cont}=rjb;
                    case 4 , param{4+cont}=repi;
                    case 5 , param{4+cont}=rseis;
                    case 6 , param{4+cont}=rx;
                    case 7 , param{4+cont}=ry0;
                    case 8 , param{4+cont}=zhyp;
                    case 9 , param{4+cont}=rztor;
                    case 10, param{4+cont}=rzbor;
                    case 11, param{4+cont}=rzbot;
                end
            end
            
            if strcmpi(f.tag,'param')
                switch f.type
                    case 'string'
                        param{4+cont}=gmpe.usp.(args{cont});
                    case 'double'
                        if isnumeric(gmpe.usp.(args{cont}))
                            param{4+cont}=gmpe.usp.(args{cont});
                        else
                            param{4+cont}=str2double(gmpe.usp.(args{cont}));
                        end
                end
            end
        end
end