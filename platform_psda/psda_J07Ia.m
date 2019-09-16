function[vout]=psda_J07Ia(ky,~,Ia,varargin)
                          
% Jibson, R. W. (2007). Regression models for estimating coseismic 
% landslide displacement. Engineering geology, 91(2-4), 209-218.

%%
% Ia  in m/s
% 
logEDP =2.401*log10(Ia)-3.481*log10(ky)-3.230;
lnEDP  =logEDP *log(10); 
logsig = 0.656;
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
