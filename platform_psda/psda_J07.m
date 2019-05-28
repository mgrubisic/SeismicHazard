function[varargout]=psda_J07(ky,~,im,varargin)
                          

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
%im    = max(im,ky+eps);
%aratio = ky./im;
%------------------------------------------------------------
%mu    = displacement in cm;
logEDP=(2.401*log10(im)-3.481*log10(ky)-3.23);
lnEDP=logEDP *log(10); 
logsig = 0.656;
sig    = logsig*log(10)*ones(size(lnEDP));

if nargin==3
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
