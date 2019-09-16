function[vout]=psda_RS09M(ky,~,im,M,varargin)

% Rathje, E. M., & Saygili, G. (2009). Probabilistic assessment of 
% earthquake-induced sliding displacements of natural slopes. Bulletin 
% of the New Zealand Society for Earthquake Engineering, 42(1), 18.

%%
im    = max(im,ky+eps);
aratio = ky./im;
%------------------------------------------------------------
%mu    = displacement in cm;
e1 =  4.89;
e2 = -4.85;
e3 =-19.64;
e4 = 42.49;
e5 =-29.06;
e6 =  0.72;
e7 =  0.89;

lnEDP=e1+e2*(aratio)+e3*(aratio).^2+e4*(aratio).^3+e5*(aratio).^4+e6*log(im)+e7*(M-6); 
std=0.732+0.789*(aratio)-0.539*(aratio).^2;
sig    = std.*ones(size(lnEDP));

if nargin==4
   vout=[lnEDP,sig];
   return
end

y    = varargin{1};
dist = varargin{2};
if strcmp(dist,'pdf')
    vout = lognpdf(y,lnEDP,sig);
elseif strcmp(dist,'cdf')
    vout = logncdf(y,lnEDP,sig);
elseif strcmp(dist,'ccdf')
    vout = 1-logncdf(y,lnEDP,sig);
end
