function[lamdaD]=psda_BMT2017_cdmM(d, Ts_param, ky_param, im, M,dPm,Cz,realSa,realD) %psda_BT2007_pce

% Bray, J. D., Macedo, J., & Travasarou, T. (2017). Simplified procedure 
% for estimating seismic slope displacements for subduction zone 
% earthquakes. Journal of Geotechnical and Geoenvironmental Engineering, 
% 144(3), 04017124.

%%
Nd       = length(d);   % N° of displacement values of lambdaD
Nim      = length(im);  % N° of acceleration values of lambdaSa

if length(Ts_param)==2, Ts_param = [Ts_param,100]; end; NTs = Ts_param(3);
if length(ky_param)==2, ky_param = [ky_param,100]; end; Nky = ky_param(3);

[ky,Ts]  = meshgrid(trlognpdf_psda(ky_param),trlognpdf_psda(Ts_param));
ky       = ky(:)';
Ts       = Ts(:)';

xrnd = normrnd(0,1,1, realSa*realD)';
zrnd = normrnd(0,1,1, realSa*realD)';

Nmag = length(M);
lnd  = zeros(Nim,Nky*NTs,Nmag);
Tlow = Ts<0.1;
for m=1:Nmag
    for j=1:Nim
        lnd(j, Tlow,m) = -5.864-3.353*log(ky( Tlow))-0.390*(log(ky( Tlow)))^2+0.538*log(ky( Tlow))*log(im(j))+3.060*log(im(j))-0.225*(log(im(j))).^2-9.421*Ts( Tlow)+0.55*M(m);
        lnd(j,~Tlow,m) = -6.896-3.353*log(ky(~Tlow))-0.390*(log(ky(~Tlow)))^2+0.538*log(ky(~Tlow))*log(im(j))+3.060*log(im(j))-0.225*(log(im(j))).^2+3.801*Ts(~Tlow)+0.55*M(m)-0.803*Ts(~Tlow)^2;
    end
end
mean_d  = mean(lnd, 2);
std_d   = std(lnd,[], 2);
sigmaD  = 0.67;

PC_term_0_GM = Cz(1:Nim, 1)';
PC_term_1_GM = Cz(1:Nim, 2)';
PC_term_2_GM = Cz(1:Nim, 3)';
PC_term_3_GM = Cz(1:Nim, 4)';
PC_term_4_GM = Cz(1:Nim, 5)';

Hazard_PC_samples = zeros(realSa*realD, Nim);
for i = 1:Nim
    Hazard_PC_samples(:, i) = ...
        PC_term_0_GM(i) + ...
        PC_term_1_GM(i) * zrnd + ...
        PC_term_2_GM(i) * (zrnd.^2-1) + ...
        PC_term_3_GM(i) * (zrnd.^3 - 3*zrnd) + ...
        PC_term_4_GM(i) * (zrnd.^4 - 6*zrnd.^2 + 3);
end

dlambdaPC        = diff(-prctile(Hazard_PC_samples,50),1,2);
dlambdaPC(end+1) = dlambdaPC(end);

PC_term_0_array = zeros(Nim,Nd);
PC_term_1_array = zeros(Nim,Nd);
PC_term_2_array = zeros(Nim,Nd);
PC_term_3_array = zeros(Nim,Nd);
PC_term_4_array = zeros(Nim,Nd);

PC_term_0_summed = zeros(1, Nd);
PC_term_1_summed = zeros(1, Nd);
PC_term_2_summed = zeros(1, Nd);
PC_term_3_summed = zeros(1, Nd);
PC_term_4_summed = zeros(1, Nd);

logd = log(d(:));
for m=1:length(M)
    for i = 1:Nim
        A_s = - std_d(i,1,m)^2/(2*sigmaD^2) - 1/2;
        B_s = (logd - mean_d(i,1,m))*std_d(i,1,m)/(sigmaD^2);
        C_s = -(logd - mean_d(i,1,m)).^2 * 1/(2*sigmaD^2);
        PC_term_0_array(i,:) = 1/1 * (1 - normcdf((logd - mean_d(i,1,m))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambdaPC(i);
        PC_term_1_array(i,:) = 1/1 * (std_d(i,1,m)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(i);
        PC_term_2_array(i,:) = 1/2 * (std_d(i,1,m)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(i);
        PC_term_3_array(i,:) = 1/6 * (std_d(i,1,m)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(i);
        PC_term_4_array(i,:) = 1/24* (std_d(i,1,m)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(i);
    end
    
    PC_term_0_summed = PC_term_0_summed + sum(PC_term_0_array ,1) * dPm(m); % sums over the rows
    PC_term_1_summed = PC_term_1_summed + sum(PC_term_1_array ,1) * dPm(m);
    PC_term_2_summed = PC_term_2_summed + sum(PC_term_2_array ,1) * dPm(m);
    PC_term_3_summed = PC_term_3_summed + sum(PC_term_3_array ,1) * dPm(m);
    PC_term_4_summed = PC_term_4_summed + sum(PC_term_4_array ,1) * dPm(m);
end

lamdaD = zeros(realSa*realD,Nd);

for i = 1:Nd
    lamdaD(:, i) =...
        PC_term_0_summed(i) * (xrnd.^0) + ...
        PC_term_1_summed(i) * (xrnd.^1) + ...
        PC_term_2_summed(i) * (xrnd.^2-1) + ...
        PC_term_3_summed(i) * (xrnd.^3-3*xrnd) + ...
        PC_term_4_summed(i) * (xrnd.^4-6*xrnd.^2+3);
end

