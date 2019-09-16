function vout = psda_RA2011F(ky,Ts,kvmax,kmax,varargin)

% Rathje, E. M., & Antonakos, G. (2011). A unified model for predicting 
% earthquake-induced sliding displacements of rigid and flexible slopes.
% Engineering Geology, 122(1-2), 51-60.

%% Coefficients of 2GM displacement model
y1  = -1.56;
y2  = -4.58;
y3  = -20.84;
y4  =  44.75;
y5  = -30.5;
y6  = -0.64;
y7  =  1.55;

% kvmax must be input in cm/s
[kvmax,kmax] = meshgrid(kvmax(:),kmax(:));

kj  = ky./kmax;
lnd = y1+y2*kj+y3*kj.^2+y4*kj.^3+y5*kj.^4+y6*log(kmax)+y7*log(kvmax)+0.71*(Ts>0.5) + 1.42*Ts*(Ts<=0.5);
sig = 0.4 + 0.284 * kj;

if nargin==4
   vout=[lnd,sig];
   return
end

d     = varargin{1};
dist  = varargin{2};
Nd    = length(d);

switch dist
    case 'convolute'
        vout    = zeros(1,length(d));
        MRDkvkm = varargin{3};
        for k = 1:Nd
            xhat    = (log(d(k))-lnd)./sig;
            ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
            dlambda = ccdf.*MRDkvkm.*(kmax>=ky);
            vout(k) = sum(dlambda(:));
        end
    case 'pdf'
        vout = lognpdf(d,lnd,sig);
    case 'cdf'
        vout = logncdf(d,lnd,sig);
    case 'ccdf'
        vout = 1-logncdf(d,lnd,sig);
end
