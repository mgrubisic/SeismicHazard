function [lny, sigma,tau,sig] = ASK_2014_nga(To,M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, mechanism, event, Z10, Vs30, Vs30type, reg)

% Norman A. Abrahamson, Walter J. Silva, and Ronnie Kamai (2014) Summary 
% of the ASK14 Ground Motion Relation for Active Crustal Regions. 
% Earthquake Spectra: August 2014, Vol. 30, No. 3, pp. 1025-1055.

if  and(To<0 || To> 10,To~=-1)
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    sig   = nan(size(M));
    IM    = IM2str(To);
    h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    uiwait(h);
    return
end

if To>=0
    To      = max(To,0.001); %PGA is associated to To=0.01;
end
period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0.001	-1];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6); % Identify the period

if T_lo==T_hi
    [lny,sigma,tau] = gmpe(index,M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, mechanism, event, Z10, Vs30, Vs30type, reg);
    sig              = sqrt(sigma.^2-tau.^2);
else
    [lny_lo,sigma_lo,tau_lo] = gmpe(index,  M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, mechanism, event, Z10, Vs30, Vs30type, reg);
    [lny_hi,sigma_hi,tau_hi] = gmpe(index+1,M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, mechanism, event, Z10, Vs30, Vs30type, reg);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma    = [sigma_lo,sigma_hi]';
    Y_tau      = [tau_lo,tau_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    sigma      = interp1(x,Y_sigma,log(To))';
    tau        = interp1(x,Y_tau,log(To))';
    sig        = sqrt(sigma.^2-tau.^2);
end

function[lny,sigma,tau]=gmpe(index,M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, mechanism, event, Z10, Vs30, Vs30type, reg)

switch mechanism
    case 'strike-slip',     frv = 0; fnm=0;
    case 'normal',          frv = 0; fnm=1;
    case 'normal-oblique',  frv = 0; fnm=0;
    case 'reverse',         frv = 1; fnm=0;
    case 'reverse-oblique', frv = 1; fnm=0;
    case 'thrust',          frv = 1; fnm=0;
end

switch reg
    case 'global',     region=0;
    case 'california', region=1;
    case 'japan',      region=2;
    case 'china',      region=3;
    case 'italy',      region=4;
    case 'turkey',     region=5;
end

switch Vs30type
    case 'measured', FVS30=1; % thereis an error in bakers_ASK_nga, see line
    case 'inferred', FVS30=0;
end

switch event
    case 'mainshock' , fas=1;
    case 'aftershock', fas=0;
end

if Ztor == 999
    if frv == 1
        Ztor = max(2.704 - 1.226 .* max(M-5.849,0),0).^2;
    else
        Ztor = max(2.673 - 1.136 .* max(M-4.970,0),0).^2;
    end
end

if W == 999
    W = min(18./sin(deg2rad(delta)),10.^(-1.75+0.45.*M));
end

if ischar(Z10) && strcmp(Z10,'unk')
    if region == 2
        Z10 = exp(-5.23 ./ 2 .* log((Vs30 ^ 2 + 412 ^ 2) ./ (1360 ^ 2 + 412 ^ 2))) ./ 1000;
    else %'non-japanese
        Z10 = exp((-7.67 ./ 4) .* log((Vs30 ^ 4 + 610 ^ 4) ./ (1360 ^ 4 + 610 ^ 4))) ./ 1000;
    end
end

[lny,sigma,tau] = ASK_2014_sub_1(index,M, Rrup, Rjb, Rx, Ry0, Ztor, delta, W, frv, fnm, fas, Z10, Vs30, FVS30, region);

function [ln_Sa, sigma,tau]=ASK_2014_sub_1 (index,M, R_RUP, R_JB, Rx, Ry0, Ztor, delta, W, F_RV, F_NM, F_AS, Z10, Vs30, FVS30, region)

HW = (Rx>=0);

% Coefficients
period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0	-1];

Vlin =[	660.0000	680.0000	770.0000	915.0000	960.0000	910.0000	740.0000	590.0000	495.0000	430.0000	360.0000	340.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	660.0000	330.0000];
b	=[-1.4700	-1.4590	-1.3900	-1.2190	-1.1520	-1.2300	-1.5870	-2.0120	-2.4110	-2.7570	-3.2780	-3.5990	-3.8000	-3.5000	-2.4000	-1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	-1.4700	-2.0200];
n	=[1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000];
M1	=[6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.8200	6.9200	7.0000	7.0600	7.1450	7.2500	6.7500	6.7500];
c	=[2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2400.0000];
c4	=[4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000];
a1	=[0.5870	0.5980	0.6020	0.7070	0.9730	1.1690	1.4420	1.6370	1.7010	1.7120	1.6620	1.5710	1.2990	1.0430	0.6650	0.3290	-0.0600	-0.2990	-0.5620	-0.8750	-1.3030	-1.9280	0.5870	5.9750];
a2	=[-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7650	-0.7110	-0.6340	-0.5290	-0.7900	-0.9190];
a3	=[0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750];
a4	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000];
a5	=[-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100];
a6	=[2.1541	2.1461	2.1566	2.0845	2.0285	2.0408	2.1208	2.2241	2.3124	2.3383	2.4688	2.5586	2.6821	2.7630	2.8355	2.8973	2.9061	2.8888	2.8984	2.8955	2.8700	2.8431	2.1541	2.3657];
a7	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a8	=[-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0220	-0.0300	-0.0380	-0.0450	-0.0550	-0.0650	-0.0950	-0.1100	-0.1240	-0.1380	-0.1720	-0.1970	-0.2180	-0.2350	-0.2550	-0.2850	-0.0150	-0.0940];
a10	=[1.7350	1.7180	1.6150	1.3580	1.2580	1.3100	1.6600	2.2200	2.7700	3.2500	3.9900	4.4500	4.7500	4.3000	2.6000	0.5500	-0.9500	-0.9500	-0.9300	-0.9100	-0.8700	-0.8000	1.7350	2.3600];
a11	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a12	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.2000	-0.2000	-0.2000	-0.1000	-0.1000];
a13	=[0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.5800	0.5600	0.5300	0.5000	0.4200	0.3500	0.2000	0.0000	0.0000	0.0000	0.0000	0.0000	0.6000	0.2500];
a14	=[-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.2400	-0.1900	-0.1100	-0.0400	0.0700	0.1500	0.2700	0.3500	0.4600	0.5400	0.6100	0.6500	0.7200	0.8000	-0.3000	0.2200];
a15	=[1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.0300	0.9200	0.8400	0.6800	0.5700	0.4200	0.3100	0.1600	0.0500	-0.0400	-0.1100	-0.1900	-0.3000	1.1000	0.3000];
a17	=[-0.0072	-0.0073	-0.0075	-0.0080	-0.0089	-0.0095	-0.0095	-0.0086	-0.0074	-0.0064	-0.0043	-0.0032	-0.0025	-0.0025	-0.0022	-0.0019	-0.0015	-0.0010	-0.0010	-0.0010	-0.0010	-0.0010	-0.0072	-0.0005];
a43	=[0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1400	0.1700	0.2200	0.2600	0.3400	0.4100	0.5100	0.5500	0.4900	0.4200	0.1000	0.2800];
a44	=[0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0700	0.1000	0.1400	0.1700	0.2100	0.2500	0.3000	0.3200	0.3200	0.3200	0.2750	0.2200	0.0500	0.1500];
a45 =[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0300	0.0600	0.1000	0.1400	0.1700	0.2000	0.2200	0.2300	0.2300	0.2200	0.2000	0.1700	0.1400	0.0000	0.0900];
a46	=[-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0300	0.0000	0.0300	0.0600	0.0900	0.1300	0.1400	0.1600	0.1600	0.1600	0.1400	0.1300	0.1000	0.0900	0.0800	-0.0500	0.0700];
a25	=[-0.0015	-0.0015	-0.0016	-0.0020	-0.0027	-0.0033	-0.0035	-0.0033	-0.0029	-0.0027	-0.0023	-0.0020	-0.0010	-0.0005	-0.0004	-0.0002	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	-0.0015	-0.0001];
a28	=[0.0025	0.0024	0.0023	0.0027	0.0032	0.0036	0.0033	0.0027	0.0024	0.0020	0.0010	0.0008	0.0007	0.0007	0.0006	0.0003	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0025	0.0005];
a29	=[-0.0034	-0.0033	-0.0034	-0.0033	-0.0029	-0.0025	-0.0025	-0.0031	-0.0036	-0.0039	-0.0048	-0.0050	-0.0041	-0.0032	-0.0020	-0.0017	-0.0020	-0.0020	-0.0020	-0.0020	-0.0020	-0.0020	-0.0034	-0.0037];
a31	=[-0.1503	-0.1479	-0.1447	-0.1326	-0.1353	-0.1128	0.0383	0.0775	0.0741	0.2548	0.2136	0.1542	0.0787	0.0476	-0.0163	-0.1203	-0.2719	-0.2958	-0.2718	-0.2517	-0.1400	-0.0216	-0.1503	-0.1462];
a36	=[0.2650	0.2550	0.2490	0.2020	0.1260	0.0220	-0.1360	-0.0780	0.0370	-0.0910	0.1290	0.3100	0.5050	0.3580	0.1310	0.1230	0.1090	0.1350	0.1890	0.2150	0.1500	0.0920	0.2650	0.3770];
a37	=[0.3370	0.3280	0.3200	0.2890	0.2750	0.2560	0.1620	0.2240	0.2480	0.2030	0.2320	0.2520	0.2080	0.2080	0.1080	0.0680	-0.0230	0.0280	0.0310	0.0240	-0.0700	-0.1590	0.3370	0.2120];
a38	=[0.1880	0.1840	0.1800	0.1670	0.1730	0.1890	0.1080	0.1150	0.1220	0.0960	0.1230	0.1340	0.1290	0.1520	0.1180	0.1190	0.0930	0.0840	0.0580	0.0650	0.0000	-0.0500	0.1880	0.1570];
a39	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a40	=[0.0880	0.0880	0.0930	0.1330	0.1860	0.1600	0.0680	0.0480	0.0550	0.0730	0.1430	0.1600	0.1580	0.1450	0.1310	0.0830	0.0700	0.1010	0.0950	0.1330	0.1510	0.1240	0.0880	0.0950];
a41	=[-0.1960	-0.1940	-0.1750	-0.0900	0.0900	0.0060	-0.1560	-0.2740	-0.2480	-0.2030	-0.1540	-0.1590	-0.1410	-0.1440	-0.1260	-0.0750	-0.0210	0.0720	0.2050	0.2850	0.3290	0.3010	-0.1960	-0.0380];
a42	=[0.0440	0.0610	0.1620	0.4510	0.5060	0.3350	-0.0840	-0.1780	-0.1870	-0.1590	-0.0230	-0.0290	0.0610	0.0620	0.0370	-0.1430	-0.0280	-0.0970	0.0150	0.1040	0.2990	0.2430	0.0440	0.0650];
s1	=[0.7540	0.7600	0.7810	0.8100	0.8100	0.8100	0.8010	0.7890	0.7700	0.7400	0.6990	0.6760	0.6310	0.6090	0.5780	0.5550	0.5480	0.5270	0.5050	0.4770	0.4570	0.4290	0.7540	0.6620];
s2	=[0.5200	0.5200	0.5200	0.5300	0.5400	0.5500	0.5600	0.5650	0.5700	0.5800	0.5900	0.6000	0.6150	0.6300	0.6400	0.6500	0.6400	0.6300	0.6300	0.6300	0.6300	0.6300	0.5200	0.5100];
s3	=[0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.3800];
s4	=[0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3800];
s1_m =[0.7410	0.7470	0.7690	0.7980	0.7980	0.7950	0.7730	0.7530	0.7290	0.6930	0.6440	0.6160	0.5660	0.5410	0.5060	0.4800	0.4720	0.4470	0.4250	0.3950	0.3780	0.3590	0.7410	0.6600];
s2_m =[0.5010	0.5010	0.5010	0.5120	0.5220	0.5270	0.5190	0.5140	0.5130	0.5190	0.5240	0.5320	0.5480	0.5650	0.5760	0.5870	0.5760	0.5650	0.5680	0.5710	0.5750	0.5850	0.5010	0.5100];
s5_JP=[0.5400	0.5400	0.5500	0.5600	0.5700	0.5700	0.5800	0.5900	0.6100	0.6300	0.6600	0.6900	0.7300	0.7700	0.8000	0.8000	0.8000	0.7600	0.7200	0.7000	0.6700	0.6400	0.5400	0.5800];
s6_JP=[0.6300	0.6300	0.6300	0.6500	0.6900	0.7000	0.7000	0.7000	0.7000	0.7000	0.7000	0.7000	0.6900	0.6800	0.6600	0.6200	0.5500	0.5200	0.5000	0.5000	0.5000	0.5000	0.6300	0.5300];

To     = period (index);
Vlin   = Vlin   (index);
b      = b      (index);
n      = n      (index);
M1     = M1     (index);
c      = c      (index);
c4     = c4     (index);
a1     = a1     (index);
a2     = a2     (index);
a3     = a3     (index);
a4     = a4     (index);
a5     = a5     (index);
a6     = a6     (index);
a7     = a7     (index);
a8     = a8     (index);
a10    = a10    (index);
a11    = a11    (index);
a12    = a12    (index);
a13    = a13    (index);
a14    = a14    (index);
a15    = a15    (index);
a17    = a17    (index);
a43    = a43    (index);
a44    = a44    (index);
a45    = a45    (index);
a46    = a46    (index);
a25    = a25    (index);
a28    = a28    (index);
a29    = a29    (index);
a31    = a31    (index);
a36    = a36    (index);
a37    = a37    (index);
a38    = a38    (index);
a39    = a39    (index);
a40    = a40    (index);
a41    = a41    (index);
a42    = a42    (index);
s1     = s1     (index);
s2     = s2     (index);
s3     = s3     (index);
s4     = s4     (index);
s1_m   = s1_m   (index);
s2_m   = s2_m   (index);
s5_JP  = s5_JP  (index);
s6_JP  = s6_JP  (index);

M2 = 5;
CRjb = 999.9;

% Term f1 - Basic form
c4m1 = c4.*ones(length(M),1);
c4m2 =c4-(c4-1).*(5-M);
c4m =ones(length(M),1);
c4m(M>4)=c4m2(M>4);
c4m(M>5)=c4m1(M>5);

R = sqrt(R_RUP.^2+ c4m.^2);

f1a = a1 + a5.*(M - M1) + a8.*(8.5 - M).^2 + (a2 + a3.*(M - M1)).*log(R) + a17.*R_RUP;
f1b = a1 + a4.*(M - M1) + a8.*(8.5 - M).^2 + (a2 + a3.*(M - M1)).*log(R) + a17.*R_RUP;
f1 = a1 + a4.*(M2 - M1) + a8.*(8.5 - M2).^2 + a6.*(M - M2) + a7.*(M - M2).^2 + (a2 + a3.*(M2 - M1)).*log(R) + a17.*R_RUP;
f1(M>=M2)=f1b(M>=M2);
f1(M>M1)=f1a(M>M1);

% term f4 - Hanging wall model
R1 = W .* cos(deg2rad(delta));
R2 = 3 .* R1;
Ry1 = Rx .* tan(deg2rad(20));
h1 = 0.25;
h2 = 1.5;
h3 = -0.75;

if delta > 30
    T1 = (90- delta)./45;
else
    T1 = 60./45;
end

a2hw = 0.2;

T2a = 1 + a2hw .* (M - 6.5);
T2b = 1 + a2hw .* (M - 6.5) - (1 - a2hw) .* (M - 6.5).^2;
T2 = zeros(length(M),1);
T2(M>5.5)=T2b(M>5.5);
T2(M>6.5)=T2a(M>6.5);

T3a = h1 + h2.*(Rx./R1) + h3.*(Rx./R1).^2;
T3b = 1 - (Rx - R1)./(R2 - R1);
T3 = zeros(length(M),1);
T3(Rx<R2)=T3b(Rx<R2);
T3(Rx<=R1)=T3a(Rx<=R1);

T4a = 1 - Ztor.^2./100;
T4 = zeros(length(M),1);
T4(Ztor<10)=T4a(Ztor<10);

T5 = R_JB*0;

for i=1:length(M)
    if Ry0(i) == 999 || Ry0(i) == 0
        if R_JB(i) == 0
            T5(i,1) = 1;
        elseif R_JB(i) <30
            T5(i,1) = 1 - R_JB(i)/30;
        end
    else
        if Ry0(i) - Ry1(i) <= 0
            T5(i,1) = 1;
        elseif Ry0(i) - Ry1(i) < 5
            T5(i,1) = 1- (Ry0(i)-Ry1(i))/5;
        end
    end
end

f4 = a13 * T1 * T2 .* T3 .* T4 .* T5 .*HW;

% Term f6 - Depth to top rupture model
f6a = a15.*Ztor./20;
f6 = a15.*ones(length(M),1);
f6(Ztor<20)=f6a(Ztor<20);

% Term: f7 and f8 - Style of Faulting
f7a = a11.*ones(length(M),1);
f7b = a11.*(M - 4);
f7 = zeros(length(M),1);
f7(M>4)=f7b(M>4);
f7(M>5)=f7a(M>5);

f8a = a12.*ones(length(M),1);
f8b = a12.*(M - 4);
f8 = zeros(length(M),1);
f8(M>4)=f8b(M>4);
f8(M>5)=f8a(M>5);

if To <= 0.5
    V1 = 1500;
elseif To < 3
    V1 = exp(-0.35.*log(To./0.5)+log(1500));
else
    V1 = 800;
end

if Vs30 < V1
    Vs30s = Vs30;
else
    Vs30s = V1;
end

if 1180 >= V1
    Vs30star1180 = V1;
else
    Vs30star1180 = 1180;
end

% term  Regional:
Ftw = (region == 6);
Fcn = (region == 3);
Fjp = (region == 2);

% Japan
if Vs30 < 150
    y1 = a36;
    y2 = a36;
    x1 = 50;
    x2 = 150;
elseif Vs30 < 250
    y1 = a36;
    y2 = a37;
    x1 = 150;
    x2 = 250;
elseif Vs30 < 350
    y1 = a37;
    y2 = a38;
    x1 = 250;
    x2 = 350;
elseif Vs30 < 450
    y1 = a38;
    y2 = a39;
    x1 = 350;
    x2 = 450;
elseif Vs30 < 600
    y1 = a39;
    y2 = a40;
    x1 = 450;
    x2 = 600;
elseif Vs30 < 850
    y1 = a40;
    y2 = a41;
    x1 = 600;
    x2 = 850;
elseif Vs30 < 1150
    y1 = a41;
    y2 = a42;
    x1 = 850;
    x2 = 1150;
else
    y1 = a42;
    y2 = a42;
    x1 = 1150;
    x2 = 3000;
end
f13Vs30 = y1 + (y2 - y1) ./ (x2 - x1) .* (Vs30 - x1);

% Taiwan
f12Vs30 = a31 .* log(Vs30s./Vlin);
f12Vs30_1180 = a31 .* log(Vs30star1180./Vlin);

Regional = Ftw .* (f12Vs30 + a25 .* R_RUP) + Fcn .* (a28 .* R_RUP) + Fjp .* (f13Vs30 + a29 .* R_RUP);
Regional_1180 = Ftw .* (f12Vs30_1180 + a25 .* R_RUP) + Fcn .* (a28 .* R_RUP) + Fjp .* (f13Vs30 + a29 .* R_RUP);

% Term f5 - site response model
%Sa 1180
f5_1180 = (a10 + b .* n) .* log(Vs30star1180 ./ Vlin);

Sa1180 = exp(f1 + f6 + F_RV.*f7 + F_NM .* f8 +  HW .* f4 + f5_1180 + Regional_1180);

if Vs30 >= Vlin
    f5 = (a10+ b.*n).*log(Vs30s./Vlin);
else
    f5 = a10.*log(Vs30s./Vlin) - b.*log(Sa1180 + c) + b.*log(Sa1180 + c.*(Vs30s./Vlin).^n);
end

% Term f10 - soil depth model
if region ~= 2 % california
    Z1ref = 1/1000.* exp(-7.67./4 .* log((Vs30.^4 + 610^4)./(1360^4 + 610^4)));
else  % Japan
    Z1ref =  1/1000.* exp(-5.23/2 .* log((Vs30.^2 + 412^2)./(1360^2 + 412^2)));
end

if Vs30 <= 150
    y1z = a43;
    y2z = a43;
    x1z = 50;
    x2z = 150;
elseif Vs30 <= 250
    y1z = a43;
    y2z = a44;
    x1z = 150;
    x2z = 250;
elseif Vs30 <= 400
    y1z = a44;
    y2z = a45;
    x1z = 250;
    x2z = 400;
elseif Vs30 <= 700
    y1z = a45;
    y2z = a46;
    x1z = 400;
    x2z = 700;
else
    y1z = a46;
    y2z = a46;
    x1z = 700;
    x2z = 1000;
end

%'f10 term goes to zero at 1180 m/s (reference)
f10 = (y1z + (Vs30 - x1z) .* (y2z - y1z) ./ (x2z - x1z)) .* log((Z10 + 0.01) ./ (Z1ref + 0.01));

% term f11 - Aftershock scaling
if CRjb <= 5
    f11 = a14;
elseif CRjb < 15
    f11 = a14.*(1-(CRjb-5)./10);
else
    f11 = 0;
end

if F_AS == 0
    f11 = 0;
end

% Sa
ln_Sa = f1+ f6 + F_RV*f7 + F_NM * f8 +  HW .* f4 + F_AS * f11 + f5 + f10 + Regional;

% Standard deviation
if FVS30 == 1 % measured
    s1 = s1_m;
    s2 = s2_m;
end

phi_ALa = s1.*ones(length(M),1);
phi_ALb = s1+(s2-s1)./2.*(M-4);
phi_AL = s2.*ones(length(M),1);
phi_AL(M<=6)= phi_ALb(M<=6);
phi_AL(M<4)= phi_ALa(M<4);

tau_ALa = s3.*ones(length(M),1);
tau_ALb = s3 + (s4-s3)./2.*(M-5);
tau_AL = s4.*ones(length(M),1);
tau_AL(M<=7)= tau_ALb(M<=7);
tau_AL(M<5)= tau_ALa(M<5);

tau_B = tau_AL;

if Fjp == 1
    phi_ALa = s5_JP.*ones(length(M),1);
    phi_ALb = s5_JP + (s6_JP - s5_JP)./50 .* (R_RUP - 30);
    phi_AL = s6_JP.*ones(length(M),1);
    phi_AL(R_RUP<=80)=phi_ALb(R_RUP<=80);
    phi_AL(R_RUP<30)=phi_ALa(R_RUP<30);
end

phi_amp = 0.4;

phi_B = sqrt(phi_AL.^2 - phi_amp.^2);

if Vs30 >= Vlin
    dln = 0;
else
    dln = -b.*Sa1180./(Sa1180 + c) + b.*Sa1180./(Sa1180 + c.*(Vs30./Vlin).^n);
end

phi = sqrt(phi_B.^2 .* (1 + dln).^2 + phi_amp.^2);

tau = tau_B.*(1+ dln);

sigma = sqrt(phi.^2+ tau.^2);


