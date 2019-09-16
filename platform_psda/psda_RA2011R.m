function vout = psda_RA2011R(ky, ~,pgv  , pga,varargin)

% Rathje, E. M., & Antonakos, G. (2011). A unified model for predicting 
% earthquake-induced sliding displacements of rigid and flexible slopes.
% Engineering Geology, 122(1-2), 51-60.

%% el GMM entrega PGV en m/s, aca lo estoy pasando a cm/s
% pgv must be input in cm/s

% Coefficients of 2GM displacement model
y1  = -1.560;
y2  = -4.583;
y3  = -20.841;
y4  =  44.754;
y5  = -30.504;
y6  = -0.640;
y7  =  1.550;
[PGVz,PGAz]=meshgrid(pgv,pga);

kj  = ky./PGAz;
lnd = y1+y2*kj+y3*kj.^2+y4*kj.^3+y5*kj.^4+y6*log(PGAz)+y7*log(PGVz);
sig = 0.41+0.52*kj;

if nargin==4
    vout=[lnd,sig];
    return
end

d     = varargin{1};
dist  = varargin{2};
Nd    = length(d);

switch dist
    case 'convolute'
        vout  = zeros(1,Nd);
        MRD   = varargin{3};
        for k = 1:Nd
            xhat    = (log(d(k))-lnd)./sig;
            ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
            dlambda = ccdf.*MRD.*(ky<=PGAz);
            vout(k)   = sum(dlambda(:));
        end
    case 'pdf'
        vout = lognpdf(d,lnd,sig);
    case 'cdf'
        vout = logncdf(d,lnd,sig);
    case 'ccdf'
        vout = 1-logncdf(d,lnd,sig);
end