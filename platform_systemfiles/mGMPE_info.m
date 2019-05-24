function [un,IM,Rmetric,Residuals]=mGMPE_info(model)
%
% [un,IM,Rmetric,Residuals]=mGMPE_info(model)
%
% Function that retuns information about GMPEs
%
% un = units of GMPE.
% For Sa or PGA set 'un' as the following:
%      un=1       for GMPE that provide IM in dimensionless units
%      un=9.8066  for GMPE that provide IM in m/s2
%      un=980.66  for GMPE that provide IM in cm/s2
%
% IM = Intensity Measure Id.
%    = To  for Sa       (Pseudoacceleration at period To)
%    =  0  for PGA      (Peak Ground Acceleration)
%    = -1  for PGV      (Peak Ground Velocity)
%    = -2  for MMI      (Modified Mercaly Intensity)
%    = -3  for Duration (Duration)
%    = -4  for H/V      (Ratio H/V)
%    = -5  for Ia       (Arias Intensity)
%    = -10 for PGD      (Permanent Ground Displacement)
%
% Rmetric is a vector of the distance metrics used by the GMPE (1 or 0)
% if Rmetric(i)=1 means that the i-th metric will be computed, and 0 will
% not be computed. The metris available are:
%
% Rmetric    = [Rrup Rhyp Rjb Repi Zseis Rx Ry0 zhyp ztor zbor zbot]
% 1.-  rrup  = closest diatance from site to rupture plane
% 2.-  rhyp  = distance from the site to the hypocenter
% 3.-  rjb   = Joyner-Boore distance, i.e., closest distance from site to
%              surface projection of rupture area
% 4.-  repi  = epicentral distance
% 5.-  rseis = Shortest distance between the recording site and the presumed
%              zone of seismogenic rupture on the fault
% 6.-  rx    = Horizontal distance from top of rupture measured perpendicular
%              to fault strike
% 7.-  ry0   = Horizontal distance off the end of the rupture measured
%              parallel to strike
% 8.-  zhyp  = depth of hypocenter
% 9.-  ztor  = Depth to top of coseismic rupture (km)
% 10.- zbor  = Depth to the bottom of the rupture plane (km)
% 11.- zbot  = The depth to the bottom of the seismogenic crust (km)

%
% Residuals = Probability distribution of IM rsiduals, e.g., lognormal,
%             beta, etc.
switch lower(model)
    case 'youngs1997',                 un=1;      IM = [0 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3];
    case 'atkinsonboore2003',          un=980.66; IM = [0 0.04 0.1 0.2 0.4 1 2 1/0.33];
    case 'zhao2006',                   un=980.66; IM = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 2 2.5 3 4 5];
    case 'mcverry2006',                un=1;      IM = [0 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0];
    case 'bchydro2012',                un=1;      IM = [0 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'arteta2018',                 un=1;      IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'montalvabastias2017',        un=1;      IM = [0 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'siberrisk2019',              un=1;      IM = [-1 0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];    
    case 'idini2016',                  un=1;      IM = [0 0.01 0.02 0.03 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'contrerasboroschek2012',     un=1;      IM = [0 0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 2];
    case 'jaimes2006',                 un=980.66; IM = [0.01 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 3 4 5 6];
    case 'jaimes2015',                 un=980.66; IM = [0 0.1 0.2 0.4 0.5 0.6 0.8 1 1.4 2 3 4 5 6];
    case 'jaimes2016',                 un=980.66; IM = [0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 10];
    case 'garciajaimes2017',           un=980.66; IM = [0.01 0.06 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 3];
    case 'sadigh1997',                 un=1;      IM = [0.00 0.075 0.10 0.20 0.30 0.40 0.50 0.75 1.00 1.50 2.00 3.00 4.00];
    case 'idriss2008_nga',             un=1;      IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'chiouyoungs2008_nga',        un=1;      IM = [-1 0 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0];
    case 'booreatkinson_2008_nga',     un=1;      IM = [0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'campbellbozorgnia_2008_nga', un=1;      IM = [-10 -1 0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'abrahamsonsilva2008_nga',    un=1;      IM = [-1 0 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'as1997h',                    un=1;      IM = [0.01 0.02 0.03 0.04 0.05 0.06 0.075 0.09 0.1 0.12 0.15 0.17 0.2 0.24 0.3 0.36 0.4 0.46 0.5 0.6 0.75 0.85 1 1.5 2 3 4 5];
    case 'campbell1997h',              un=1;      IM = [0.05 0.075 0.10 0.15 0.20 0.30 0.50 0.75 1.00 1.50 2.00 3.00 4.00];
    case 'i_2014_nga',                 un=1;      IM = [0 0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10];
    case 'cy_2014_nga',                un=1;      IM = [-1 0 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.12 0.15 0.17 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'bssa_2014_nga',              un=1;      IM = [-1 0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'cb_2014_nga',                un=1;      IM = [-1 0 0.010 0.020 0.030 0.050 0.075 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10];
    case 'ask_2014_nga',               un=1;      IM = [-1 0 0.01 0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10];
    case 'condsa',                     un=1;      IM = 1;
    case 'condpgv',                    un=1;      IM = -1;
    case 'macedo2019',                 un=1;      IM = -5;
    case 'pce_ngawest2',               un=1;      IM = 0.2; % discuss with Jorge
    case 'pce_bchydro',                un=1;      IM = 0.2; % discuss with Jorge
    case 'udm',                        un=[];     IM = [];
end

if nargout>=3
    switch lower(model)
        case 'youngs1997',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'atkinsonboore2003',          Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'zhao2006',                   Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'mcverry2006',                Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'bchydro2012',                Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'arteta2018',                 Rmetric=[0 1 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'montalvabastias2017',        Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'siberrisk2019',              Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'idini2016',                  Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'contrerasboroschek2012',     Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2006',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2015',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2016',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'garciajaimes2017',           Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'sadigh1997',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'idriss2008_nga',             Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'chiouyoungs2008_nga',        Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'booreatkinson_2008_nga',     Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'campbellbozorgnia_2008_nga', Rmetric=[1 0 1 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'abrahamsonsilva2008_nga',    Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'as1997h',                    Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'campbell1997h',              Rmetric=[0 0 0 0 1 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'i_2014_nga',                 Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'cy_2014_nga',                Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'bssa_2014_nga',              Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'cb_2014_nga',                Rmetric=[1 0 1 0 0 1 0 1 1 0 0];  Residuals = 'lognormal';
        case 'ask_2014_nga',               Rmetric=[1 0 1 0 0 1 1 0 1 0 0];  Residuals = 'lognormal';
        case 'condsa',                     Rmetric=[0 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'condpgv',                    Rmetric=[0 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'macedo2019',                 Rmetric=[0 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'pce_ngawest2',               Rmetric=[1 0 1 0 0 1 1 0 1 0 0];  Residuals = 'lognormal';
        case 'pce_bchydro',                Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'udm',                        Rmetric=[];                       Residuals = '';
    end
end

