function[vout]=psda_BT2007M(ky,Ts,im,M,varargin)

% Model with Magnitude Dependence
%Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating 
%earthquake-induced deviatoric slope displacements. Journal of Geotechnical
% and Geoenvironmental Engineering, 133(4), 381-392.

% ky = yield coefficient
% Ts = fundamental period
% im = Sa(1.5Ts) spetral acceleration at the degraded system period
if Ts<0.05
    lnEDP = -0.22-2.83*log(ky)-0.333*(log(ky))^2+0.566*log(ky)*log(im)+3.04*log(im)-0.244*(log(im)).^2+1.5*Ts+0.278*(M-7);
else
    lnEDP = -1.10-2.83*log(ky)-0.333*(log(ky))^2+0.566*log(ky)*log(im)+3.04*log(im)-0.244*(log(im)).^2+1.5*Ts+0.278*(M-7);
end

sig   = 0.66;
Pzero = 1-normcdf(-1.76-3.22*log(ky)-0.484*Ts*log(ky)+3.52*log(im),0,1);    
y     = max(varargin{1},1);
dist  = varargin{2};
if strcmp(dist,'pdf')
    if y > 1
        vout = (1-Pzero).*lognpdf(y,lnEDP,sig);
    elseif y==1
        vout = Pzero;
    end

elseif strcmp(dist,'cdf')
    if y >= 1
        vout = Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig));
    elseif y==1
        vout = Pzero;
    end    
elseif strcmp(dist,'ccdf')
    if y >= 1
        vout = 1- (Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig)));
    elseif y==1
        vout = 1-Pzero;
    end    
end


