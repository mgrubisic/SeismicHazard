function[varargout]=psda_RS09M(ky,~,im,M,varargin)

% ky = yield acceleration
% im = PGA, intensity measure
%
% Reference:
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
e1=4.89;
e2=-4.85;
e3=-19.64;
e4=42.49;
e5=-29.06;
e6=0.72;
e7=0.89;

lnEDP=e1+e2*(aratio)+e3*(aratio).^2+e4*(aratio).^3+e5*(aratio).^4+e6*log(im)+e7*(M-6); 
std=0.732+0.789*(aratio)-0.539*(aratio).^2;
sig    = std.*ones(size(lnEDP));

if nargin==4
    % empty varargin returns medium values and the standard deviation
    varargout{1}=lnEDP;
    varargout{2}=sig;
    return
end

y    = varargin{1};
dist = varargin{2};
if strcmp(dist,'pdf')
    varargout{1} = lognpdf(y,lnEDP,sig);
elseif strcmp(dist,'cdf')
    varargout{1} = logncdf(y,lnEDP,sig);
elseif strcmp(dist,'ccdf')
    varargout{1} = 1-logncdf(y,lnEDP,sig);
end
