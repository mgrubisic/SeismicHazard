function hd = psda_RathjeR(d,ky,im,MRD)

im(:,1) = im(:,1)*100; % el GMM entrega PGV en m/s, aca lo estoy pasando a cm/s

% Coefficients of 2GM displacement model
y1  = -1.560;
y2  = -4.583;
y3  = -20.841;
y4  =  44.754;
y5  = -30.504;
y6  = -0.640;
y7  =  1.550;
Nd  = length(d);
[PGVz,PGAz]=meshgrid(im(:,1),im(:,2));


kj  = ky./PGAz;
mu  = y1+y2*kj+y3*kj.^2+y4*kj.^3+y5*kj.^4+y6*log(PGAz)+y7*log(PGVz);
sig = 0.41+0.52*kj;
hd  = zeros(1,Nd);
for k = 1:Nd
    xhat    = (log(d(k))-mu)./sig;
    ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
    dlambda = ccdf.*MRD.*(ky<=PGAz);
    hd(k)   = sum(dlambda(:));
end