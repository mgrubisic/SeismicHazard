function[lny,sigma,tau,sig]=GarciaJaimes2017(To,M,rrup)

% García-Soto & Jaimes (2017). Ground-Motion Prediction Model for Vertical
% Response Spectra from Mexican Interplate Earthquakes. Bulletin of the
% Society of America, 107(2), 887-900
%
% To        = spectral period
% M         = moment magnitude
% rrup      = closest distance to fault rupture
% h         = focal depth (km)
% mechanism ='interface'
% media     = 'rock'

if  To<0 || To> 3
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    sig   = nan(size(M));
    %IM    = IM2str(To);
    %h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    %uiwait(h);
    return
end

To      = max(To,0.01); %PGA is associated to To=0.01; 
period  = [0.01 0.06 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 3];
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

DATA = [0.5115 0.9092  -0.500  -0.0119  0.78
    1.8411 0.8082 -0.500  -0.0139 0.87
    1.8188 0.8454 -0.500 -0.0137 0.83
    0.6324 0.9828 -0.500 -0.0112 0.76
    -0.3584 1.0822 -0.500 -0.0094 0.73
    -1.1473 1.1551 -0.500 -0.0081 0.73
    -1.8420 1.2211 -0.500 -0.0069 0.73
    -2.9901 1.3359 -0.500 -0.0057 0.72
    -3.8785 1.3874 -0.500 -0.0046 0.69
    -5.2292 1.4958 -0.500 -0.0040 0.70
    -6.5005 1.6227 -0.500 -0.0039 0.64
    -7.9826 1.7396 -0.500 -0.0041 0.65];

C       = DATA(index,:);
lny     = C(1)+ C(2)*M + C(3)*log(rrup) + C(4)*rrup;
sigma   = C(5)*ones(size(M));