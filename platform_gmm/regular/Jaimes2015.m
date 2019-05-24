function[lny,sigma,tau,sig]=Jaimes2015(To,M,rrup)
% 
% Jaimes M.A., Ramirez-Gaytán A. & Reinoso E (2015).Ground-Motion Prediction Model From Intermediate-Depth
% Intraslab Earthquakes at the Hill and Lake-Bed Zones of Mexico City. Journal of Earthquake Engineering,
% 19(8), 1260-1278, 2015
%
% To        = spectral period
% M         = moment magnitude
% rrup      = closest distance to fault rupture
% h         = focal depth (km)
% mechanism ='intraslab'
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

To      = max(To,0.001); %PGA is associated to To=0.001; 
period  = [0.001 0.1 0.2 0.4 0.5 0.6 0.8 1 1.4 2 3 4 5 6];
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

DATA = [-1.7918 1.6132  -1.00  -0.0058  0.60
    -1.2779 1.5902 -1.00  -0.0061 0.55
    -1.0515 1.6707 -1.00  -0.0064 0.60
    -1.1719 1.6122 -1.00  -0.0053 0.65
    -0.9565 1.5647 -1.00  -0.0047 0.66
    -1.6592 1.6529 -1.00  -0.0043 0.68
    -2.1755 1.6601 -1.00  -0.0030 0.69
    -2.8579 1.7282 -1.00  -0.0026 0.71
    -4.2636 1.9042 -1.00  -0.0030 0.70
    -5.8253 2.0204 -1.00  -0.0012 0.73
    -7.9145 2.2219 -1.00  -0.0006 0.69
    -8.6208 2.2856 -1.00  -0.0010 0.68
    -9.5327 2.3711 -1.00  -0.0009 0.66
    -9.7489 2.3438 -1.00  -0.0005 0.60];

C       = DATA(index,:);
lny     = C(1)+ C(2)*M + C(3)*log(rrup) + C(4)*rrup;
sigma   = C(5)*ones(size(M));