function[vout]=psda_BM2019M(ky,Ts,im,M,varargin)
% Bray, J. D., Macedo, J. L.(2019). Procedure for Estimating Shear-Induced Seismic
% Slope Displacement for Shallow Crustal Earthquakes  (ASCE).

% ky = yield coefficient
% Ts = fundamental period
% im = Sa(1.3Ts) spetral acceleration at the degraded system period
%coefficientes for Ts<0.7 (PD=0)
a1=-2.457;
a2=-2.979;
a3=-0.124;
a4=-0.707;
a5=1.688;
a6=2.759;
%coefficientes for Ts>0.7 (PD=0)
b1=-3.531;
b2=-4.783;
b3=-0.342;
b4=-0.300;
b5=-0.672;
b6=2.658;
% model standard deviations
sig=[0.736];
%coeficientes model Ts<0.1 (Ignore the second and third columns, which
%corresponds to variations of the model and are not used)
%%%
c1=[-4.684     ,-3.1384                ,-4.9202+1.335];
c2=[-2.482      ,0.471                  ,0.4529];
c3=[-0.244      ,-1.199                 ,-1.527];
c4=[0.344       ,3.1563-12.53           ,3.1391-12.515];
c5=[2.649       ,0.6038                 ,0.6024];
c6=[-0.090      ,0                      ,0];
c7=[-9.471      ,0.373                  ,0.056];
c8=[0           ,1                      ,-0.051];
c9=[0.603       ,1                      ,1];

%coeficientes model Ts>=0.1(Ignore the second and third columns, which
%corresponds to variations of the model and are not used)
%%%
d1=[-5.981      ,-4.4778                ,-4.9202];
d2=[-2.482      ,0.471                  ,0.4529];
d3=[-0.244      ,-1.199                 ,-1.527];
d4=[0.344       ,3.1563                 ,3.1391];
d5=[2.649       ,0.6038                 ,0.6024];
d6=[-0.090      ,-0.9058                ,-0.8971];
d7=[3.223       ,0.373                  ,0.056];
d8=[-0.945     ,1                      ,-0.051];
d9=[0.603       ,1                      ,1];

if Ts<=0.1
    P =1-normcdf(a1 +a2*log(ky)+a3*log(ky)^2 +a4*Ts*log(ky)+ a5*Ts+ a6*log(im));%stata2
    Pzero=P;
    lnEDP=(c1(1)+c2(1)*log(ky)+c3(1)*(log(ky))^2+c4(1)*log(ky)*log(im)+c5(1)*log(im)+c6(1)*(log(im)).^2+c7(1)*Ts+c8(1)*(Ts)^2+c9(1)*(M));
end

if Ts>0.1
    if Ts<=0.7
        P=1-normcdf(a1 +a2*log(ky)+a3*log(ky)^2 +a4*Ts*log(ky)+ a5*Ts+ a6*log(im));%stata2
    end
    
    if Ts>0.7
        P=1-normcdf(b1 +b2*log(ky)+b3 *log(ky)^2 +b4*Ts*log(ky)+b5*Ts+ b6*log(im));%stata2
    end
    if Ts==0.7
        aa=1-normcdf(a1 +a2*log(ky)+a3*log(ky)^2 +a4*Ts*log(ky)+ a5*Ts+ a6*log(im));%stata2
        bb=1-normcdf(b1 +b2*log(ky)+b3 *log(ky)^2 +b4*Ts*log(ky)+b5*Ts+ b6*log(im));%stata2
        P=0.5*(aa+bb);
    end
    Pzero=P;
    
    lnEDP=exp(d1(1)+d2(1)*log(ky)+d3(1)*(log(ky))^2+d4(1)*log(ky)*log(im)+d5(1)*log(im)+d6(1)*(log(im)).^2+d7(1)*Ts+d8(1)*(Ts)^2+d9(1)*(M));
    
    
end
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