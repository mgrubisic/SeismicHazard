function[vout]=psda_BT2007(ky,Ts,im,varargin)

% Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating 
% earthquake-induced deviatoric slope displacements. Journal of Geotechnical
% and Geoenvironmental Engineering, 133(4), 381-392.

%%
if Ts<0.05
    lnEDP = -0.22-2.83*log(ky)-0.333*(log(ky))^2+0.566*log(ky)*log(im)+3.04*log(im)-0.244*(log(im)).^2+1.5*Ts;
else
    lnEDP = -1.10-2.83*log(ky)-0.333*(log(ky))^2+0.566*log(ky)*log(im)+3.04*log(im)-0.244*(log(im)).^2+1.5*Ts;
end

sig   = 0.67;

if nargin==3
   vout=[lnEDP,sig];
   return
end

y0    = 1;
Pzero = 1-normcdf(-1.76-3.22*log(ky)-0.484*Ts*log(ky)+3.52*log(im),0,1);    
y     = max(varargin{1},1);
dist  = varargin{2};
if strcmp(dist,'pdf')
    if y > y0
        vout = (1-Pzero).*lognpdf(y,lnEDP,sig);
    elseif y==y0
        vout = Pzero/y0;
    end

elseif strcmp(dist,'cdf')
    if y >= y0
        vout = Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig));
    elseif y==y0
        vout = Pzero;
    end    
elseif strcmp(dist,'ccdf')
    if y >= y0
        vout = 1- (Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig)));
    elseif y==y0
        vout = 1-Pzero;
    end    
end


