function[lny,sigma,tau,sig]=Jaimes2006(To,M,rrup)

% Jaimes, M.A., Reinoso, E. y Ordaz, M. (2006). Comparison of methods to predict response spectra
% at instrumented sites given the magnitude and distance of an earthquake, Journal of Earthquake Engineering,
% 10(6), 887-902
%
% To        = spectral period
% M         = moment magnitude
% rrup      = closest distance to fault rupture
% h         = focal depth (km)
% mechanism ='interface'
% media     = 'rock' for specifically Mexico City

if  To<0 || To> 6
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    sig   = nan(size(M));
    IM    = IM2str(To);
    h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    uiwait(h);
    return
end

To      = max(To,0.01); %PGA is associated to To=0.01; 
period  = [0.01 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 3 4 5 6];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6); % Identify the period

if T_lo==T_hi
    [lny,sigma] = gmpe(index,M,rrup);
else
    [lny_lo,sigma_lo] = gmpe(index,  M,rrup);
    [lny_hi,sigma_hi] = gmpe(index+1,M,rrup);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma    = [sigma_lo,sigma_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    sigma      = interp1(x,Y_sigma,log(To))';
end

tau=0*sigma;
sig=sigma;

function [lny,sigma]=gmpe(index,M,rrup)
DATA = [5.8402 1.0316  -0.50  -0.0061 0.32
    5.9134 0.9830 -0.50 -0.0059 0.30
    6.5114 0.8385 -0.50 -0.0065 0.29
    6.6971 0.9241 -0.50 -0.0069 0.30
    6.3551 0.9185 -0.50 -0.0058 0.36
    6.5482 0.9275 -0.50 -0.0060 0.33
    6.8881 1.1769 -0.50 -0.0073 0.35
    6.3906 1.2098 -0.50 -0.0059 0.41
    6.2638 1.2022 -0.50 -0.0057 0.47
    5.9080 1.3671 -0.50 -0.0056 0.49
    5.3173 1.4247 -0.50 -0.0055 0.58
    3.7525 1.5224 -0.50 -0.0032 0.57
    3.6685 1.5340 -0.50 -0.0045 0.68
    3.6034 1.6933 -0.50 -0.0060 0.81];
C     = DATA(index,:);
lny   = C(1)+ C(2)*(M-6) + C(3)*log(rrup) + C(4)*rrup;
sigma = C(5)*ones(size(M));

