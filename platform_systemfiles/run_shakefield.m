function[handles]=run_shakefield(handles)

shakefield = buildshake(handles.model,handles.scenarios,handles.opt);
opt        = handles.opt; 
Ns         = length(shakefield);
sites      = cell2mat(handles.h.p(:,2:4));
Vs30       = handles.h.Vs30;
r0         = gps2xyz(sites,handles.opt.ellipsoid);

param.opp       = 0;
param.residual  = 'tau';
param.direction = 'horizontal';
param.isVs30clustered = 1;

for i=1:Ns
    param.mechanism = shakefield.mechanism;
    shakefield(i)   = runscenario(shakefield(i),r0,Vs30,opt,param);
end
handles.shakefield=shakefield;

function shakefield = runscenario(shakefield,r0,Vs30,opt,param)
IM1         = opt.IM1;
IM2         = opt.IM2;
NumSim      = opt.NumSim;
ellipsoid   = opt.ellipsoid;
IMs         = [IM1;IM2];
NIM         = length(IMs);
Nsites      = size(r0,1);
mulogIM     = zeros(Nsites,NIM);
[tau,sig,sigma] = deal(zeros(1,NIM));

for i=1:Nsites
    for j=1:NIM
        [mulogIM(i,j),sigma(j),tau(j),sig(j),M]=run_gmpe(shakefield,r0(i,:),IMs(j),Vs30(i),ellipsoid);
    end
end

% Computes Spatial Correlation Structures
h    = computeh(r0);
Lii  = zeros(Nsites,Nsites,NIM);
for i=1:NIM
    rho  = opt.Spatial(IMs(i), h,param);
    Cii  = (tau(i)^2+sig(i)^2*rho);
    Lii(:,:,i)  = chol(Cii,'lower');
end

% Computes Spectral Correlation Structures
rhoIM   = zeros(NIM,NIM);
param.M = M;
for i=1:NIM
    for j=1:NIM
        rhoIM(i,j)=opt.Spectral(IMs(i),IMs(j),param);
    end
end

% Builds Extended Covariance Matrix
C = zeros(Nsites*NIM,Nsites*NIM);
for i=1:NIM
    II = (1:Nsites)+Nsites*(i-1);
    for j=1:NIM
        JJ       = (1:Nsites)+Nsites*(j-1);
        C(II,JJ) = rhoIM(i,j)*Lii(:,:,i)*Lii(:,:,j)';
    end
end

shakefield.mulogIM  = mulogIM(:);
shakefield.L        = chol(C,'lower');
shakefield.Z        = normrnd(0,1,[Nsites*NIM,NumSim]);

function[mu,sigma,tau,sig,M,rrup]=run_gmpe(scenario,r0,IM,Vs30,ellipsoid)

mscl     = scenario.mscl;
gmpe     = scenario.gmpe;
rateM    = mscl.dPm;
gmpefun  = gmpe.handle;
imunits  = gmpe.un;
Rmetric  = gmpe.Rmetric;
scenario = mGMPEVs30(scenario,Vs30);
switch scenario.type
    case 'point',      param = mGMPEassemblePoint(r0,scenario,Rmetric,ellipsoid);
    case 'line',       param = mGMPEassembleLine(r0,scenario,Rmetric,ellipsoid);
    case 'area',       param = mGMPEassembleArea(r0,scenario,Rmetric,ellipsoid);  
    case 'background', param = mGMPEassembleArea(r0,scenario,Rmetric,ellipsoid);
    case 'fault',      param = mGMPEassembleFault_S(r0,scenario,rateM,Rmetric,ellipsoid);
end

[mu,sigma,tau,sig] = gmpefun(IM,param{:});
mu  = mu-log(imunits); % verified! such that exp(mu) is in units of g's

% Compute Magnitude and Rrup (this is necesary since not all GMPEs use rrup
M       = scenario.mscl.M;
geom    = scenario.geom;
rf      = geom.hypm;
switch scenario.type
    case {'point','line'}
        rrup    = dist_rrup(r0,rf,0,[]);
    case {'area','background'}
        rupt    = scenario.rupt;
        rupArea = rupRelation(M,0,geom.Area,rupt.RA{1},rupt.RA{2});
        rrup    = dist_rrup(r0,rf,rupArea,geom.normal);
    case 'fault'
        rupt    = scenario.rupt;
        RA     = rupRelation(M,0,Inf,rupt.RA{1},rupt.RA{2});
        aratio = rupt.aratio;
        RL     = sqrt(RA*aratio);
        RW     = sqrt(RA/aratio);
        rrup = dist_rrup4  (r0,rf,RW,RL,geom,ellipsoid);
end
