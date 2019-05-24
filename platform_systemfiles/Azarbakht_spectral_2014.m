function [rho] = Azarbakht_spectral_2014(T1, T2,param)

% Azarbakht, A., Mousavi, M., Nourizadeh, M. and Shahri, M., 2014. 
% Dependence of correlations between spectral accelerations at multiple 
% periods on magnitude and distance. Earthquake Engineering & Structural 
% Dynamics, 43(8), 1193-1204.
%
% ---- help for CP_spectral_2018 ---
% Developed by G.Candia, March 2018. Santiago.
% National Research Center for Integrated Natural Disaster Management
% CONICYT/FONDAP/15110017. FONDECYT N� 1170836 �SIBER-RISK: SImulation
% Based Earthquake Risk and Resilience of Interdependent Systems and NetworKs

T = [0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1 1.5 2 2.5 3 4 5];
M = param.M;
R = param.R;

if and(M <  6.33,R <  19.47), V=Bin1; end
if and(M <  6.33,R >= 19.47), V=Bin2; end
if and(M >= 6.33,R <  19.47), V=Bin3; end
if and(M >= 6.33,R >= 19.47), V=Bin4; end

T_min = min(T1, T2); 
T_max = max(T1, T2); 
[X,Y] = meshgrid(T);
rho   = interp2(X,Y,V,T_min,T_max,'spline');
rho(T_min==T_max)=1;
rho   = min(rho,1);
rho   = max(rho,min(V(:)));
rho(T_min<0.05) = nan;
rho(T_max>5)    = nan;

function V=Bin1()
V=[1 0.93 0.88 0.83 0.71 0.63 0.60 0.55 0.45 0.45 0.42 0.41 0.41 0.47 0.49
0.93 1 0.91 0.83 0.68 0.59 0.53 0.51 0.39 0.38 0.34 0.33 0.33 0.38 0.40
0.88 0.91 1 0.89 0.75 0.62 0.57 0.55 0.46 0.45 0.44 0.41 0.41 0.43 0.44
0.83 0.83 0.89 1 0.83 0.71 0.68 0.65 0.57 0.53 0.51 0.48 0.48 0.51 0.52
0.71 0.68 0.75 0.83 1 0.87 0.77 0.69 0.65 0.59 0.58 0.55 0.57 0.57 0.57
0.63 0.59 0.62 0.71 0.87 1 0.90 0.77 0.75 0.67 0.61 0.59 0.60 0.63 0.65
0.60 0.53 0.57 0.68 0.77 0.90 1 0.84 0.78 0.64 0.60 0.58 0.59 0.63 0.63
0.55 0.51 0.55 0.65 0.69 0.77 0.84 1 0.91 0.82 0.77 0.74 0.74 0.73 0.74
0.45 0.39 0.46 0.57 0.65 0.75 0.78 0.91 1 0.88 0.82 0.78 0.81 0.80 0.79
0.45 0.38 0.45 0.53 0.59 0.67 0.64 0.82 0.88 1 0.93 0.88 0.87 0.85 0.83
0.42 0.34 0.44 0.51 0.58 0.61 0.60 0.77 0.82 0.93 1 0.95 0.92 0.88 0.84
0.41 0.33 0.41 0.48 0.55 0.59 0.58 0.74 0.78 0.88 0.95 1 0.97 0.93 0.90
0.41 0.33 0.41 0.48 0.57 0.60 0.59 0.74 0.81 0.87 0.92 0.97 1 0.96 0.93
0.47 0.38 0.43 0.51 0.57 0.63 0.63 0.73 0.80 0.85 0.88 0.93 0.96 1 0.98
0.49 0.40 0.44 0.52 0.57 0.65 0.63 0.74 0.79 0.83 0.84 0.90 0.93 0.98 1];

function V=Bin2()
V=[  1 0.89 0.83 0.74 0.58 0.47 0.39 0.27 0.19 0.07 0.04 0.04 0.08 0.15 0.22
  0.89 1 0.84 0.73 0.53 0.42 0.32 0.23 0.15 0.04 0.02 0.02 0.05 0.12 0.18 
 0.83 0.84 1 0.87 0.65 0.54 0.46 0.34 0.26 0.16 0.13 0.14 0.17 0.23 0.29   
 0.74 0.73 0.87 1 0.76 0.60 0.52 0.38 0.30 0.19 0.16 0.17 0.19 0.25 0.30    
 0.58 0.53 0.65 0.76 1 0.85 0.74 0.57 0.48 0.39 0.36 0.35 0.36 0.40 0.45    
 0.47 0.42 0.54 0.60 0.85 1 0.87 0.70 0.59 0.50 0.44 0.42 0.42 0.45 0.48    
 0.39 0.32 0.46 0.52 0.74 0.87 1 0.84 0.73 0.63 0.56 0.53 0.53 0.55 0.57    
 0.27 0.23 0.34 0.38 0.57 0.70 0.84 1 0.88 0.76 0.68 0.65 0.65 0.65 0.65   
 0.19 0.15 0.26 0.30 0.48 0.59 0.73 0.88 1 0.86 0.78 0.75 0.76 0.72 0.71      
 0.07 0.04 0.16 0.19 0.39 0.50 0.63 0.76 0.86 1 0.91 0.85 0.82 0.76 0.73    
 0.04 0.02 0.13 0.16 0.36 0.44 0.56 0.68 0.78 0.91 1 0.95 0.90 0.84 0.80      
 0.04 0.02 0.14 0.17 0.35 0.42 0.53 0.65 0.75 0.85 0.95 1 0.96 0.90 0.87    
 0.08 0.05 0.17 0.19 0.36 0.42 0.53 0.65 0.76 0.82 0.90 0.96 1 0.94 0.91      
 0.15 0.12 0.23 0.25 0.40 0.45 0.55 0.65 0.72 0.76 0.84 0.90 0.94 1 0.97      
 0.22 0.18 0.29 0.30 0.45 0.48 0.57 0.65 0.71 0.73 0.80 0.87 0.91 0.97 1 ];

function V=Bin3()
V=[ 1 0.92 0.87 0.83 0.76 0.72 0.62 0.58 0.46 0.25 0.16 0.10 0.04 -0.08 -0.10   
 0.92 1 0.92 0.85 0.73 0.66 0.52 0.46 0.33 0.15 0.06 0.02 -0.02 -0.07 -0.09
 0.87 0.92 1 0.93 0.79 0.69 0.53 0.46 0.35 0.18 0.07 0.03 -0.01 -0.05 -0.08  
 0.83 0.85 0.93 1 0.85 0.74 0.61 0.50 0.41 0.23 0.14 0.09 0.04 -0.01 -0.03    
 0.76 0.73 0.79 0.85 1 0.87 0.73 0.62 0.54 0.37 0.25 0.19 0.16 0.06 0.05      
 0.72 0.66 0.69 0.74 0.87 1 0.86 0.71 0.60 0.46 0.31 0.23 0.17 0.02 0.03      
 0.62 0.52 0.53 0.61 0.73 0.86 1 0.79 0.68 0.53 0.40 0.33 0.27 0.10 0.10      
 0.58 0.46 0.46 0.50 0.62 0.71 0.79 1 0.84 0.67 0.57 0.48 0.42 0.24 0.22     
 0.46 0.33 0.35 0.41 0.54 0.60 0.68 0.84 1 0.83 0.72 0.61 0.58 0.40 0.36        
 0.25 0.15 0.18 0.23 0.37 0.46 0.53 0.67 0.83 1 0.85 0.73 0.69 0.53 0.46      
 0.16 0.06 0.07 0.14 0.25 0.31 0.40 0.57 0.72 0.85 1 0.92 0.84 0.72 0.66        
 0.10 0.02 0.03 0.09 0.19 0.23 0.33 0.48 0.61 0.73 0.92 1 0.94 0.84 0.80      
 0.04 -0.02 -0.01 0.04 0.16 0.17 0.27 0.42 0.58 0.69 0.84 0.94 1 0.91 0.86      
-0.08 -0.07 -0.05 -0.01 0.06 0.02 0.10 0.24 0.40 0.53 0.72 0.84 0.91 1 0.94    
-0.10 -0.09 -0.08 -0.03 0.05 0.03 0.10 0.22 0.36 0.46 0.66 0.80 0.86 0.94 1];

function V=Bin4()
V=[1 0.93 0.89 0.86 0.82 0.79 0.75 0.64 0.54 0.44 0.33 0.28 0.26 0.21 0.15
0.93 1 0.92 0.86 0.77 0.71 0.65 0.53 0.42 0.32 0.24 0.20 0.17 0.15 0.11 
0.89 0.92 1 0.92 0.81 0.73 0.65 0.52 0.38 0.29 0.21 0.17 0.15 0.14 0.10   
0.86 0.86 0.92 1 0.86 0.78 0.69 0.54 0.41 0.31 0.22 0.18 0.17 0.15 0.11    
0.82 0.77 0.81 0.86 1 0.88 0.81 0.67 0.53 0.41 0.31 0.27 0.25 0.23 0.18    
0.79 0.71 0.73 0.78 0.88 1 0.90 0.74 0.61 0.49 0.38 0.32 0.31 0.27 0.22    
0.75 0.65 0.65 0.69 0.81 0.90 1 0.84 0.71 0.58 0.47 0.41 0.40 0.35 0.29    
0.64 0.53 0.52 0.54 0.67 0.74 0.84 1 0.86 0.72 0.61 0.55 0.51 0.44 0.38   
0.54 0.42 0.38 0.41 0.53 0.61 0.71 0.86 1 0.83 0.72 0.65 0.61 0.52 0.43      
0.44 0.32 0.29 0.31 0.41 0.49 0.58 0.72 0.83 1 0.88 0.79 0.74 0.66 0.57    
0.33 0.24 0.21 0.22 0.31 0.38 0.47 0.61 0.72 0.88 1 0.91 0.86 0.76 0.69      
0.28 0.20 0.17 0.18 0.27 0.32 0.41 0.55 0.65 0.79 0.91 1 0.93 0.84 0.76    
0.26 0.17 0.15 0.17 0.25 0.31 0.40 0.51 0.61 0.74 0.86 0.93 1 0.90 0.81      
0.21 0.15 0.14 0.15 0.23 0.27 0.35 0.44 0.52 0.66 0.76 0.84 0.90 1 0.93      
0.15 0.11 0.10 0.11 0.18 0.22 0.29 0.38 0.43 0.57 0.69 0.76 0.81 0.93 1];