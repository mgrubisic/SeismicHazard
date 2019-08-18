function [lny,sigma,tau,phi] = Macedo2019(T,M,func,mechanism,un,varargin)

% M         = Moment magnitude
% func      = function handle to Sa model
% mechanism = 1 for interface earthquakes
%             2 for intraslab earthquakes
% un        = units of IM
% varargin  = input parameters that follow T in the matlab function that
%              runs the Sa GMPE

if T~=-5
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    %IM    = IM2str(T);
    %h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    %uiwait(h);
    return
end

switch mechanism
    case 1 %interface
        % south america
        c1 = 0.95;
        c2 =-0.36;
        c3 = 0.53;
        c4 = 1.54;
        c5 = 0.17;
        
        s1 = 0.130;
        s2 = 2.370;
        s3 = 0.030;
        s4 = 0.270;
        
        phi = 0.30;
        tau = 0.19;
      
    case 2 %intraslab
        c1 =-0.75;
        c2 =-0.24;
        c3 = 0.66;
        c4 = 1.58;
        c5 = 0.14;
        
        s1 = 0.190;
        s2 = 2.500;
        s3 = 0.020;
        s4 = 0.230;        
        
        phi = 0.30;
        tau = 0.14;        
end

[lnPGA,sPGA]  = func(0,varargin{:}); % evaluates GMM for PGA
[lnSa1,sSa1]  = func(1,varargin{:}); % evaluates GMM for Sa(T=1)

% correct for units
lnPGA  = lnPGA - log(un);
lnSa1  = lnSa1 - log(un);

% extracts Vs30
Vs30 = varargin{6};

% evaluates conditional GMM
lny   = c1+c2*log(Vs30)+c3*M+c4*lnPGA+c5*lnSa1;
sigma = sqrt(s1+s2*sPGA.^2+s3*sSa1.^2+s4*sPGA.*sSa1);
