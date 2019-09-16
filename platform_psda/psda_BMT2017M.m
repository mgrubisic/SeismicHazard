function[vout]=psda_BMT2017M(ky,Ts,im,M,varargin)

% Bray, J. D., Macedo, J., & Travasarou, T. (2017). Simplified procedure 
% for estimating seismic slope displacements for subduction zone 
% earthquakes. Journal of Geotechnical and Geoenvironmental Engineering, 
% 144(3), 04017124.

%% FROM JON'S SPREADSHEET
if Ts<0.7
    Pzero = 1-normcdf(-2.640-3.200*log(ky)-0.170*log(ky)^2-0.490*Ts*log(ky)+2.094*Ts+2.908*log(im),0,1);
else
    Pzero = 1-normcdf(-3.531-4.783*log(ky)-0.342*log(ky)^2-0.300*Ts*log(ky)-0.672*Ts+2.658*log(im),0,1);
end
if Ts<0.1
    lnEDP = -5.864-3.353*log(ky)-0.390*(log(ky))^2+0.538*log(ky)*log(im)+3.060*log(im)-0.225*(log(im)).^2-9.421*Ts+0.55*M-0*Ts^2;
else
    lnEDP = -6.896-3.353*log(ky)-0.390*(log(ky))^2+0.538*log(ky)*log(im)+3.060*log(im)-0.225*(log(im)).^2+3.801*Ts+0.55*M-0.803*Ts^2;
end

sig   = 0.73;

if nargin==4
   vout=[lnEDP,sig];
   return
end

y0    = 0.5;
y     = max(varargin{1},y0);
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


