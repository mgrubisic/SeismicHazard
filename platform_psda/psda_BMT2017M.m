function[vout]=psda_BMT2017M(ky,Ts,im,M,varargin)

% Bray, J. D., Macedo, J. L., & Travasarou, T. (2017). Simplified 
% Procedure for Estimating Seismic Slope Deviatoric Displacements 
% in Subduction Zones. Journal Geotechnical and Geoenvironmental 
% Engineering (ASCE).

% ky = yield coefficient
% Ts = fundamental period
% im = Sa(1.5Ts) spetral acceleration at the degraded system period


%% FROM peru PRESENTATION
% if Ts<0.7
%     Pzero = 1-normcdf(-2.75-3.30*log(ky)-0.18*log(ky)^2-0.56*Ts*log(ky)+1.94*Ts+2.95*log(im),0,1);
% else
%     Pzero = 1-normcdf(-3.77-5.17*log(ky)-0.40*log(ky)^2-0.43*Ts*log(ky)-1.03*Ts+2.91*log(im),0,1);
% end
% 
% if Ts<0.05
%     lnEDP = -6.37-3.05*log(ky)-0.33*(log(ky))^2+0.45*log(ky)*log(im)+2.61*log(im)-0.23*(log(im)).^2+1.41*Ts+0.64*M;
% else
%     lnEDP = -6.97-3.05*log(ky)-0.33*(log(ky))^2+0.45*log(ky)*log(im)+2.61*log(im)-0.23*(log(im)).^2+1.41*Ts+0.64*M;
% end
% sig   = 0.73;

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

%%

y     = max(varargin{1},0.5);
dist  = varargin{2};
if strcmp(dist,'pdf')
    if y > 0.5
        vout = (1-Pzero).*lognpdf(y,lnEDP,sig);
    elseif y==0.5
        vout = Pzero;
    end

elseif strcmp(dist,'cdf')
    if y >= 0.5
        vout = Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig));
    elseif y==0.5
        vout = Pzero;
    end    
elseif strcmp(dist,'ccdf')
    if y >= 0.5
        vout = 1- (Pzero + (1-Pzero).*(logncdf(y,lnEDP,sig)));
    elseif y==0.5
        vout = 1-Pzero;
    end    
end


