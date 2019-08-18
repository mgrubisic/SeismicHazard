function [lnY,sigma] = PCE_bchydro(To,M,Rrup,h,mechanism,region,DeltaC1,Vs30)

% computes sigma from ASK2014

if  To<0 || To> 10
    lnY   = nan(size(M));
    sigma = nan(size(M));
    %IM    = IM2str(To);
    %h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    %uiwait(h);
    return
end

index = 1; % set to PGA temporarily
[Vlin,b,t1,t2,t6,t7,t8,t10,t11,t12,t13,t14,t15,t16,~,~,~,n,c,t3,t4,t5,t9,C4,C1]=getcoefs(index);
[~,sigma]=BCHydro2012(To,M,Rrup,h,mechanism,region,DeltaC1,Vs30);


% DeltaC1
if strcmp(mechanism,'interface')
    logP = log([eps;0.3;0.5;1;2;3;10]);
    DC1data = [0	0.2     0.4
        0	0.2     0.4
        -0.1	0.1     0.3
        -0.2	0       0.2
        -0.3    -0.1	0.1
        -0.4	-0.2	0
        -0.4	-0.2	0];
    switch DeltaC1
        case 'lower',  DC1 = interp1(logP,DC1data(:,1),log(max(To,eps)),'pchip');
        case 'central',DC1 = interp1(logP,DC1data(:,2),log(max(To,eps)),'pchip');
        case 'upper',  DC1 = interp1(logP,DC1data(:,3),log(max(To,eps)),'pchip');
        case 'none',   DC1 = 0;
    end
end

if strcmp(mechanism,'intraslab')
    switch DeltaC1
        case 'lower',  DC1 = -0.5;
        case 'central',DC1 = -0.3;
        case 'upper',  DC1 = -0.1;
        case 'none',   DC1 = 0;
    end
end

% FEVENT
switch mechanism
    case 'interface'
        Fevent = 0;
    case 'intraslab'
        Fevent = 1;
end

% FFABA
switch region
    case {'forearc','unkown'}, FFABA = 0;
    case 'backarc',            FFABA = 1;
end

% Vs*
if Vs30>1000
    Vsa = 1000;
else
    Vsa = Vs30;
end

%% RANDOMIZED coefficients
mean_t  = [t1,t2,t3,t7,t8,t11,t12,t13,t15,t16];
cov_t   = diag(abs([t1,t2,t3,t7,t8,t11,t12,t13,t15,t16]))/500; % ignoring cross-correlations
Nreal   = 2000;
t       = mvnrnd(mean_t,cov_t,Nreal);
Nscen = length(M);
lnY   = zeros(Nreal,Nscen);

switch mechanism
    case 'interface'
        for k=1:Nscen
            Mk = M(k);
            Rk = Rrup(k);
            
            t1  = t(:,1);
            t2  = t(:,2);
            t3  = t(:,3);
            t7  = t(:,4);
            t8  = t(:,5);
            t12 = t(:,7);
            t13 = t(:,8);
            t15 = t(:,9);
            t16 = t(:,10);
            
            % magnitude term
            if Mk<=(C1+DC1)
                fmag = t4*(Mk-(C1+DC1))+t13*(10-Mk)^2;
            else
                fmag = t5*(Mk-(C1+DC1))+t13*(10-Mk)^2;
            end
            
            % site term
            fsite = t12*log(Vsa/Vlin)+b*n*log(Vsa/Vlin);
            
            % fore arc back arc
            if Fevent ==1
                fFABA = (t7+t8*log(max(Rk,85)/40))*FFABA;
            else
                fFABA = (t15+t16*log(max(Rk,100)/40))*FFABA;
            end
            
            lnY(:,k) = t1+t4*DC1+(t2+t3*(Mk-7.8)).*log(Rk+C4*exp(t9*(Mk-6)))+t6*Rk+fmag+fFABA+fsite;
        end
        
    case 'intraslab'
        for k=1:Nscen
            Mk = M(k);
            hk = h(k);
            Rk = Rrup(k);
            
            t1  = t(:,1);
            t2  = t(:,2);
            t3  = t(:,3);
            t7  = t(:,4);
            t8  = t(:,5);
            t11 = t(:,6);
            t12 = t(:,7);
            t13 = t(:,8);
            t15 = t(:,9);
            t16 = t(:,10);
            
            % magnitude term
            if Mk<=(C1+DC1)
                fmag = t4*(Mk-(C1+DC1))+t13*(10-Mk).^2;
            else
                fmag = t5*(Mk-(C1+DC1))+t13*(10-Mk).^2;
            end
            
            % depth
            fdepth = t11*(min(hk,120)-60)*Fevent;
            
            % site term
            fsite = t12*log(Vsa/Vlin)+b*n*log(Vsa/Vlin);
            
            % fore arc back arc
            if Fevent ==1
                fFABA = (t7+t8*log(max(Rk,85)/40))*FFABA;
            else
                fFABA = (t15+t16*log(max(Rk,100)/40))*FFABA;
            end
            
            lnY(:,k) = t1+t4*DC1+(t2+t14*Fevent+t3*(Mk-7.8)).*log(Rk+C4*exp(t9*(Mk-6)))+t6*Rk+t10*Fevent+fmag+fdepth+fFABA+fsite;
        end
end


function[Vlin,b,t1,t2,t6,t7,t8,t10,t11,t12,t13,t14,t15,t16,phi,tau,sigma,n,c,t3,t4,t5,t9,C4,C1]=getcoefs(index)

COEFF = [865.1	-1.186	4.2203	-1.35	-0.0012	1.0988	-1.42	3.12	0.013	0.98	-0.0135	-0.4	0.9996	-1	0.6	0.43	0.74
    865.1	-1.186	4.2203	-1.35	-0.0012	1.0988	-1.42	3.12	0.013	0.98	-0.0135	-0.4	0.9996	-1	0.6	0.43	0.74
    1053.5	-1.346	4.5371	-1.4	-0.0012	1.2536	-1.65	3.37	0.013	1.288	-0.0138	-0.4	1.103	-1.18	0.6	0.43	0.74
    1085.7	-1.471	5.0733	-1.45	-0.0012	1.4175	-1.8	3.37	0.013	1.483	-0.0142	-0.4	1.2732	-1.36	0.6	0.43	0.74
    1032.5	-1.624	5.2892	-1.45	-0.0012	1.3997	-1.8	3.33	0.013	1.613	-0.0145	-0.4	1.3042	-1.36	0.6	0.43	0.74
    877.6	-1.931	5.4563	-1.45	-0.0014	1.3582	-1.69	3.25	0.013	1.882	-0.0153	-0.4	1.26	-1.3	0.6	0.43	0.74
    748.2	-2.188	5.2684	-1.4	-0.0018	1.1648	-1.49	3.03	0.0129	2.076	-0.0162	-0.35	1.223	-1.25	0.6	0.43	0.74
    654.3	-2.381	5.0594	-1.35	-0.0023	0.994	-1.3	2.8	0.0129	2.248	-0.0172	-0.31	1.16	-1.17	0.6	0.43	0.74
    587.1	-2.518	4.7945	-1.28	-0.0027	0.8821	-1.18	2.59	0.0128	2.348	-0.0183	-0.28	1.05	-1.06	0.6	0.43	0.74
    503	-2.657	4.4644	-1.18	-0.0035	0.7046	-0.98	2.2	0.0127	2.427	-0.0206	-0.23	0.8	-0.78	0.6	0.43	0.74
    456.6	-2.669	4.0181	-1.08	-0.0044	0.5799	-0.82	1.92	0.0125	2.399	-0.0231	-0.19	0.662	-0.62	0.6	0.43	0.74
    430.3	-2.599	3.6055	-0.99	-0.005	0.5021	-0.7	1.7	0.0124	2.273	-0.0256	-0.16	0.58	-0.5	0.6	0.43	0.74
    410.5	-2.401	3.2174	-0.91	-0.0058	0.3687	-0.54	1.42	0.012	1.993	-0.0296	-0.12	0.48	-0.34	0.6	0.43	0.74
    400	-1.955	2.7981	-0.85	-0.0062	0.1746	-0.34	1.1	0.0114	1.47	-0.0363	-0.07	0.33	-0.14	0.6	0.43	0.74
    400	-1.025	2.0123	-0.77	-0.0064	-0.082	-0.05	0.7	0.01	0.408	-0.0493	0	0.31	0	0.6	0.43	0.74
    400	-0.299	1.4128	-0.71	-0.0064	-0.2821	0.12	0.7	0.0085	-0.401	-0.061	0	0.3	0	0.6	0.43	0.74
    400	0	0.9976	-0.67	-0.0064	-0.4108	0.25	0.7	0.0069	-0.723	-0.0711	0	0.3	0	0.6	0.43	0.74
    400	0	0.6443	-0.64	-0.0064	-0.4466	0.3	0.7	0.0054	-0.673	-0.0798	0	0.3	0	0.6	0.43	0.74
    400	0	0.0657	-0.58	-0.0064	-0.4344	0.3	0.7	0.0027	-0.627	-0.0935	0	0.3	0	0.6	0.43	0.74
    400	0	-0.4624	-0.54	-0.0064	-0.4368	0.3	0.7	0.0005	-0.596	-0.098	0	0.3	0	0.6	0.43	0.74
    400	0	-0.9809	-0.5	-0.0064	-0.4586	0.3	0.7	-0.0013	-0.566	-0.098	0	0.3	0	0.6	0.43	0.74
    400	0	-1.6017	-0.46	-0.0064	-0.4433	0.3	0.7	-0.0033	-0.528	-0.098	0	0.3	0	0.6	0.43	0.74
    400	0	-2.2937	-0.4	-0.0064	-0.4828	0.3	0.7	-0.006	-0.504	-0.098	0	0.3	0	0.6	0.43	0.74 ];

Coeff = COEFF(index,:);
Vlin=Coeff(1);
b=Coeff(2);
t1=Coeff(3);
t2=Coeff(4);
t6=Coeff(5);
t7=Coeff(6);
t8=Coeff(7);
t10=Coeff(8);
t11=Coeff(9);
t12=Coeff(10);
t13=Coeff(11);
t14=Coeff(12);
t15=Coeff(13);
t16=Coeff(14);
phi=Coeff(15);
tau=Coeff(16);
sigma=Coeff(17);

n  = 1.18;
c  = 1.88;
t3 = 0.1;
t4 = 0.9;
t5 = 0;
t9 = 0.4;
C4 = 10;
C1 = 7.8;
