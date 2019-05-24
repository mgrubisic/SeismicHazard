function[PCE,MCS,time1,time2,Cz]=runPCE(source,IM,site,Vs30,im,ellipsoid,realSa,RandType)

source  = mGMPEVs30(source,Vs30);
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
            case 'crustal'
                [param,rate] = mGMPEassembleCrustal(r0,source,rateM,Rmetric,ellipsoid);
        end
end

%% HAZARD INTEGRAL

sqpi  = sqrt(pi); % squared root of pi
Nscen = length(param{1});

switch RandType
    case 'random'
        zrnd  = normrnd(0,1,1, realSa);
    case 'stable'
        tol   = 1e-2;
        if realSa<=5
            x2=norminv(1-tol,0,1);
            zrnd=gaussquad(realSa,[-x2 x2]);
        else
            zrnd = norminv(linspace(tol,1-tol,realSa),0,1); % mod by GC 15/11/2018
        end
end

Psi   = [H(0,zrnd);H(1,zrnd);H(2,zrnd);H(3,zrnd);H(4,zrnd)];
Cz    = zeros(Nscen, Nim, 5);
lnzi  = log(imunits*im);
[lnY,sigma]  = gmpefun(IM,param{:});
mu  = mean(lnY,1);
smu = std (lnY,0,1);

% compute PCE
tic
for k = 1:Nscen
    muk   = mu(k);
    smuk  = smu(k);   % Assume a standard deviation for the median
    sk    = sigma(k); % Actual Standard Deviation from GMPE (ASK_nga2014)
    stot  = sqrt(sk^2+smuk^2);
    
    for  i = 1:Nim
        lnz  = lnzi(i);
        a    = -smuk^2/(2*sk^2) - 1/2;
        b    = (lnz - muk) * smuk/(sk^2);
        c    = - (lnz - muk).^2/(2*sk^2);
        d    = smuk/(2*sk*sqpi)*exp(c - b^2/(4*a));
        xhat = (lnz - muk)/stot;
        Cz(k, i, 1) = 1/1  * NMmin * 0.5*(1-erf(xhat/sqrt(2)));
        Cz(k, i, 2) = 1/1  * NMmin * d * 1/((-a)^0.5);
        Cz(k, i, 3) = 1/2  * NMmin * d * b/(2*(-a)^1.5);
        Cz(k, i, 4) = 1/6  * NMmin * d * ((-2*a*(1 + 2*a) + b^2)/(4*(-a)^2.5));
        Cz(k, i, 5) = 1/24 * NMmin * d * (-b*(6*a*(1 + 2*a)-b^2)/(8*(-a)^3.5));
    end
end
Cz  = permute(sum(bsxfun(@times,rate,Cz),1),[2 3 1]);
PCE = Cz*Psi;
PCE(PCE<0)=nan;
time1 = toc;

% compute MCS
tic
MCS = zeros(Nim, realSa);
MU = bsxfun(@plus,mu,zrnd'*smu)';
for i = 1:Nim
    xhat = bsxfun(@times,lnzi(i)-MU,1./sigma);
    ccdf = 0.5*(1-erf(xhat/sqrt(2)));
    MCS(i,:) = NMmin* rate'*ccdf;
end
time2 = toc;

% EOF
