function[vout]=psda_J07M(ky,~,im,M,varargin)
                          
% Jibson, R. W. (2007). Regression models for estimating coseismic 
% landslide displacement. Engineering geology, 91(2-4), 209-218.

%%
im    = max(im,ky+eps);
aratio = ky./im;
%------------------------------------------------------------
%mu    = displacement in cm;
logEDP=(-2.710+log10(((1-aratio).^2.335).*((aratio).^-1.478))+0.424*M);
lnEDP=logEDP *log(10); 
logsig = 0.454;
sig    = logsig*log(10)*ones(size(lnEDP));

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
