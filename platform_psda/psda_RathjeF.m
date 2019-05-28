function hd = psda_RathjeF(d,Ts,ky,im,MRDkvkm)

% Coefficients of 2GM displacement model
y1  = -1.56;
y2  = -4.58;
y3  = -20.84;
y4  =  44.75;
y5  = -30.5;
y6  = -0.64;
y7  =  1.55;
Nd  = length(d);
[kvmax,kmax] = meshgrid(im(:,1),im(:,2));

kj  = ky./kmax;
mu  = y1+y2*kj+y3*kj.^2+y4*kj.^3+y5*kj.^4+y6*log(kmax)+y7*log(kvmax)+0.71*(Ts>0.5) + 1.42*Ts*(Ts<=0.5);
sig = 0.4 + 0.284 * kj;
hd  = zeros(1,length(d));
for k = 1:Nd
    xhat    = (log(d(k))-mu)./sig;
    ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
    dlambda = ccdf.*MRDkvkm.*(kmax>=ky);
    hd(k)   = sum(dlambda(:));
end