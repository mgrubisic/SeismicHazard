function[lny,sigma,tau,sig]=Jaimes2016(To,M,rrup)

%Jaimes, M.A., Lermo, J. y García-Soto, A. (2016). Ground-Motion Prediction Model from Local Earthquakes
%of the Mexico Basin at the Hill Zone of Mexico City, Bulletin of the Seismological Society of America,
%106(6), 2532-2544
%
% To        = spectral period
% M         = moment magnitude
% rrup      = closest distance to fault rupture
% h         = focal depth (km)
% mechanism ='intraslab'
% media     = 'rock' for specifically Mexico City

if  To<0 || To> 10
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
period  = [0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 10];
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
DATA = [
    -5.8904 3.1992  -1.300  -0.0315 1.05
    -5.9225 3.2973 -1.300  -0.0355 1.13
    -5.3886 3.3471 -1.300 -0.0489 1.23
    -4.8389 3.1510 -1.300 -0.0392 1.07
    -4.9727 2.8882 -1.300 -0.0092 0.90
    -5.6202 2.9602 -1.300 -0.0057 0.78
    -6.2649 3.0069 -1.300 -0.0036 0.70
    -7.5119 3.1248 -1.300 -0.0024 0.60
    -8.2369 3.2187 -1.300 -0.0056 0.59
    -10.0150 3.5280 -1.300 -0.0004 0.51
    -10.7760 3.5746 -1.300  0.0014 0.53
    -11.4350 3.5005 -1.300  0.0039 0.55
    -11.9930 3.5312 -1.300 -0.0017 0.56
    -12.4270 3.5710 -1.300 -0.0060 0.59
    -12.9800 3.4714 -1.300 -0.0122 0.68];

C       = DATA(index,:);
lny     = C(1)+ C(2)*M + C(3)*log(rrup) + C(4)*rrup;
sigma   = C(5)*ones(size(M));