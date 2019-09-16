function[vout]=psda_RS09V(ky,~,pgv,pga,varargin)

% Rathje, E. M., & Saygili, G. (2009). Probabilistic assessment of 
% earthquake-induced sliding displacements of natural slopes. Bulletin 
% of the New Zealand Society for Earthquake Engineering, 42(1), 18.

%%
% pgv must be input in cm/s
aratio = ky./pga;
%------------------------------------------------------------
e1 = -1.56;
e2 = -4.58;
e3 =-20.84;
e4 = 44.75;
e5 =-30.50;
e6 = -0.64;
e7 =  1.55;

lnEDP=e1+e2*(aratio)+e3*(aratio).^2+e4*(aratio).^3+e5*(aratio).^4+e6*log(pga)+e7*log(pgv); 
std=0.405+0.524*aratio;
sig    = std.*ones(size(lnEDP));

if nargin==4
   vout=[lnEDP,sig];
   return
end

y    = varargin{1};
dist = varargin{2};
if strcmp(dist,'pdf')
    vout = lognpdf(y,lnEDP,sig);
elseif strcmp(dist,'cdf')
    vout = logncdf(y,lnEDP,sig);
elseif strcmp(dist,'ccdf')
    vout = 1-logncdf(y,lnEDP,sig);
end
