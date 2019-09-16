function [lnY,sigma] = PCE_bchydro(To,M,Rrup,Vs30)

if  To<0 || To> 3
    lnY   = nan(size(M));
    sigma = nan(size(M));
    return
end

load PCE_bchydro_coefs.mat COEF
Nreal      = 2000;
To         = max(To,0.01); %PGA is associated to To=0.01;
sigma      = 0.74*ones(size(M)); % this is period independent
period     = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];
T_lo       = max(period(period<=To));
T_hi       = min(period(period>=To));
index      = find(abs((period - T_lo)) < 1e-6); % Identify the period
X1         = ones(size(M));
X2         = M;
X3         = (10-M).^2;
X4         = log(Rrup+10*exp(0.4*(M-6)));
X5         = M.*log(Rrup+10*exp(0.4*(M-6)));
X6         = Rrup;
X7         = log(Vs30*ones(size(M))/760);
X          = [X1,X2,X3,X4,X5,X6,X7];

if T_lo==T_hi
    mean_theta = COEF(:,1,index);
    cov_ini    = COEF(:,2:end,index);
    c          = mvnrnd(mean_theta,cov_ini,Nreal);
    lnY        = c*X';
else
    x1    = log(T_lo);
    x2    = log(T_hi);
    x0    = log(To);
    S     = [1-(x0-x1)/(x2-x1),(x0-x1)/(x2-x1)];     
    mu1   = COEF(:,1,index);
    cov1  = COEF(:,2:end,index);
    mu2   = COEF(:,1,index+1);
    cov2  = COEF(:,2:end,index+1);
    mu    = [mu1,mu2]*S';
    cov12 = S(1)^2*cov1+S(2)^2*cov2;
    c     = mvnrnd(mu,cov12,Nreal);    
    lnY   = c*X';
end
