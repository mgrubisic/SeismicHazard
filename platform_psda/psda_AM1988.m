function[vout]=psda_AM1988(ky,~,im,varargin)

% Ambraseys, N. N., & Menu, J. M. (1988). Earthquake-induced ground 
% displacements. Earthquake engineering & structural dynamics, 
% 16(7), 985-1006.
%
%
%------------------------------------------------------------
im    = max(im,ky+eps);
aratio = ky./im;
%------------------------------------------------------------
%mu    = displacement in cm;
logEDP = 0.9+log10(((1-aratio).^2.53).*(aratio.^(-1.09)));
logsig = 0.30;
lnEDP  = logEDP *log(10);
sig    = logsig*log(10)*ones(size(lnEDP));

if nargin==3
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
