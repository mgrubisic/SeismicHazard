function[varargout]=psda_J07M(ky,~,im,M,varargin)
                          

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
logEDP=(-2.710+log10(((1-aratio).^2.335).*((aratio).^-1.478))+0.424*M);
lnEDP=logEDP *log(10); 
logsig = 0.454;
sig    = logsig*log(10)*ones(size(lnEDP));

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
