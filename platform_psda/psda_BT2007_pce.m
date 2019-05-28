function[lambdaD]=psda_BT2007_pce(d, Ts_param, ky_param, im, HAZ,param,Cz,realSa,realD)

Nd       = length(d);     % N° of displacement values of lambdaD
Nim      = length(im);    % N° of acceleration values of lambdaSa

if length(Ts_param)==2, Ts_param = [Ts_param,100]; end; NTs = Ts_param(3);
if length(ky_param)==2, ky_param = [ky_param,100]; end; Nky = ky_param(3);

[ky,Ts]  = meshgrid(trlognpdf_psda(ky_param),trlognpdf_psda(Ts_param));
ky       = ky(:)';
Ts       = Ts(:)';
Percent  = 50;
switch param.sampling
    case 'random'
        xrnd = normrnd(0,1,1, realD);
        zrnd = normrnd(0,1,1, realSa);
    case 'stable'
        tol   = 1e-2;
        xrnd  = norminv(linspace(tol,1-tol,realD)  ,0,1); % mod by GC 15/11/2018
        zrnd  = norminv(linspace(tol,1-tol,realSa)',0,1); % mod by GC 15/11/2018
    case 'custom'
        zrnd  = param.zrnd;
        xrnd  = param.xrnd;
end

switch sprintf('%s&%s',param.avgHazard,param.integration)
    case 'N&mcs0'
        % assumes normal distributoin on lnd given realizations of ky and Ts
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        sigmaD  = 0.67;
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        lambdaD = zeros(realD * realSa, length(d));
        dHAZ    = diff(HAZ,1); dHAZ(end+1,:)=dHAZ(end,:);
        
        for i = 1:length(d)
%             fprintf('MCS, displacement %g of %g\n',i,length(d))
            for j = 1:realD
                xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                IND  = (1:realSa)+realSa*(j-1);
                lambdaD(IND,i) = -dHAZ'*Pd;   % rectangular rule
            end
        end
        
    case 'N&mcs'
        % uses true distributoin on lnd given realizations of ky and Ts
        mu   = zeros(Nim,realD);
        tol  = 1e-2;
        y    = linspace(tol,1-tol,realD);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd  = zeros(1,Nky*NTs);
            lnd( Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
            [f,x]=ecdf(lnd);
            mu(j,:) = interp1(f,x,y,'linear');
        end
        sigmaD  = 0.67;
        lambdaD = zeros(realD * realSa, Nd);
        
        for i = 1:Nd
            fprintf('MCS, displacement %g of %g\n',i,Nd)
            cont = 1;
            for k = 1:realD
                xhat = (log(d(i)) - mu(:,k))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                for j = 1:realSa
                    lambdaD(cont,i) = -trapz(HAZ(:,k),Pd);
                    cont = cont + 1;
                end
            end
        end
        
    case 'N&pce'
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        sigmaD  = 0.67;
        
        PC_term_0_GM = Cz(1:Nim, 1)';
        PC_term_1_GM = Cz(1:Nim, 2)';
        PC_term_2_GM = Cz(1:Nim, 3)';
        PC_term_3_GM = Cz(1:Nim, 4)';
        PC_term_4_GM = Cz(1:Nim, 5)';
        
        Hazard_PC_samples = zeros(realSa, Nim);
        
        for i = 1:Nim
            Hazard_PC_samples(:, i) = PC_term_0_GM(i) + ...
                PC_term_1_GM(i) * zrnd + ...
                PC_term_2_GM(i) * (zrnd.^2-1) + ...
                PC_term_3_GM(i) * (zrnd.^3 - 3*zrnd) + ...
                PC_term_4_GM(i) * (zrnd.^4 - 6*zrnd.^2 + 3);
        end
        
        dlambdaPC = zeros(5, Nim);
        dlambdaPC(1, 1:end-1) = diff(-PC_term_0_GM(1, :),1,2);
        dlambdaPC(2, 1:end-1) = diff(-PC_term_1_GM(1, :),1,2);
        dlambdaPC(3, 1:end-1) = diff(-PC_term_2_GM(1, :),1,2);
        dlambdaPC(4, 1:end-1) = diff(-PC_term_3_GM(1, :),1,2);
        dlambdaPC(5, 1:end-1) = diff(-PC_term_4_GM(1, :),1,2);
        dlambdaPC(:,end) = dlambdaPC(:,end);
        
        PC_term_0_array = zeros(Nd,5,Nim);
        PC_term_1_array = zeros(Nd,5,Nim);
        PC_term_2_array = zeros(Nd,5,Nim);
        PC_term_3_array = zeros(Nd,5,Nim);
        PC_term_4_array = zeros(Nd,5,Nim);
        logd = log(d(:));
        for k = 1:5
            for i = 1:Nim
                A_s = - std_d(i)^2/(2*sigmaD^2) - 1/2;
                B_s = (logd - mean_d(i))*std_d(i)/(sigmaD^2);
                C_s = -(logd - mean_d(i)).^2 * 1/(2*sigmaD^2);
                PC_term_0_array(:,k,i) = 1/1 * (1 - normcdf((logd - mean_d(i))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambdaPC(k, i);
                PC_term_1_array(:,k,i) = 1/1 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_2_array(:,k,i) = 1/2 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_3_array(:,k,i) = 1/6 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_4_array(:,k,i) = 1/24 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
            end
            
        end % End loop over PC terms Hazard!!
        PC0 = sum(PC_term_0_array,3);
        PC1 = sum(PC_term_1_array,3);
        PC2 = sum(PC_term_2_array,3);
        PC3 = sum(PC_term_3_array,3);
        PC4 = sum(PC_term_4_array,3);
        
        HD      = zeros(5,realD);
        HSa     = zeros(5,realSa);
        for i=1:5
            HD (i,:)= H(i-1,xrnd);
            HSa(i,:)= H(i-1,zrnd');
        end
        
        lambdaD = zeros(realD*realSa, Nd);
        for i = 1:Nd
            cont = 1;
            for j = 1:realD
                HD1=HD(1,j);
                HD2=HD(2,j);
                HD3=HD(3,j);
                HD4=HD(4,j);
                HD5=HD(5,j);
                for k = 1:realSa
                    for l = 1:5 % Loop over PC terms hazard
                        Hlk = HSa(l,k);
                        lambdaD(cont, i) =lambdaD(cont, i)+Hlk*(PC0(i, l)*HD1+PC1(i,l)*HD2+PC2(i,l)*HD3+PC3(i,l)*HD4+ PC4(i,l)*HD5);
                    end
                    cont = cont + 1;
                end
            end
        end
        
    case 'Y&mcs0'
        % assumes normal distributoin on lnd given realizations of ky and Ts
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        sigmaD  = 0.67;
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        lambdaD = zeros(realD, length(d));
        HAZP   = prctile(HAZ,Percent,2);
        for i = 1:length(d)
            cont    = 1;
            for j = 1:realD
                xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                lambdaD(cont,i) = -trapz(HAZP,Pd);
                cont = cont + 1;
            end
        end
        
    case 'Y&mcs'
        % uses true distributoin on lnd given realizations of ky and Ts
        mu   = zeros(Nim,realD);
        tol  = 1e-2;
        y    = linspace(tol,1-tol,realD);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd  = zeros(1,Nky*NTs);
            lnd( Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
            [f,x]=ecdf(lnd);
            mu(j,:) = interp1(f,x,y,'linear');
        end
        sigmaD  = 0.67;
        lambdaD = zeros(realD, Nd);
        
        for i = 1:Nd
            fprintf('MCS, displacement %g of %g\n',i,Nd)
            cont = 1;
            for k = 1:realD
                xhat = (log(d(i)) - mu(:,k))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                %for j = 1:realSa
                lambdaD(cont,i) = -trapz(HAZP(:,1),Pd);%-trapz(HAZP(:,j),Pd)
                cont = cont + 1;
                %end
            end
        end
        
    case 'Y&pce'
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        sigmaD  = 0.67;
        
        PC_term_0_GM = Cz(1:Nim, 1)';
        PC_term_1_GM = Cz(1:Nim, 2)';
        PC_term_2_GM = Cz(1:Nim, 3)';
        PC_term_3_GM = Cz(1:Nim, 4)';
        PC_term_4_GM = Cz(1:Nim, 5)';
        
        Hazard_PC_samples = zeros(realSa, Nim);
        
        for i = 1:Nim
            Hazard_PC_samples(:, i) = PC_term_0_GM(i) + ...
                PC_term_1_GM(i) * zrnd + ...
                PC_term_2_GM(i) * (zrnd.^2-1) + ...
                PC_term_3_GM(i) * (zrnd.^3 - 3*zrnd) + ...
                PC_term_4_GM(i) * (zrnd.^4 - 6*zrnd.^2 + 3);
        end
        
        Hazard_PC_samplesP  = prctile(Hazard_PC_samples, Percent, 1);
        dlambdaPCP          = diff(Hazard_PC_samplesP,1);
        dlambdaPCP(1,end+1) = dlambdaPCP(1,end);
        
        PC_term_0_array = zeros(Nim, Nd);
        PC_term_1_array = zeros(Nim, Nd);
        PC_term_2_array = zeros(Nim, Nd);
        PC_term_3_array = zeros(Nim, Nd);
        PC_term_4_array = zeros(Nim, Nd);
        logd = log(d(:));
        
        for i = 1:Nim
            A_s = - std_d(i)^2/(2*sigmaD^2) - 1/2;
            B_s = (logd - mean_d(i))*std_d(i)/(sigmaD^2);
            C_s = -(logd - mean_d(i)).^2 * 1/(2*sigmaD^2);
            PC_term_0_array(i,:) = 1/1 * (1 - normcdf((logd - mean_d(i))/(sqrt(sigmaD^2 + std_d(i)^2))))*-dlambdaPCP(i);
            PC_term_1_array(i,:) = 1/1 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*-dlambdaPCP(i);
            PC_term_2_array(i,:) = 1/2 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*-dlambdaPCP(i);
            PC_term_3_array(i,:) = 1/6 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*-dlambdaPCP(i);
            PC_term_4_array(i,:) = 1/24 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*-dlambdaPCP(i);
        end
        
        PC0 = sum(PC_term_0_array,1);
        PC1 = sum(PC_term_1_array,1);
        PC2 = sum(PC_term_2_array,1);
        PC3 = sum(PC_term_3_array,1);
        PC4 = sum(PC_term_4_array,1);
        
        HD = zeros(5,realD);
        for i=1:5
            HD (i,:)= H(i-1,xrnd);
        end
        
        lambdaD = zeros(realD, Nd);
        for i = 1:Nd
            lambdaD(:, i)= PC0(i)*HD(1,:)'+PC1(i)*HD(2,:)'+PC2(i)*HD(3,:)'+PC3(i)*HD(4,:)'+PC4(i)*HD(5,:)';
        end
        
        % trucazo para poder combinarlo con el caso 'N'. Jorge, do you
        % agree??
        lambdaD = repmat(lambdaD,realSa,1);
end
