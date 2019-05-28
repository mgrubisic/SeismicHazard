function[hd]=psda_BM2019_NF(d,ky,im,MRD,mode)
%mode=2 for D100GM50 (i.e. the 100th percentil of D over all rotations)
%mode=4 for D50GM50  (i.e. the 50th percentil of D over all rotations)
% There are two additional sub functions in this function
% ALLBMP18 sets the coefficients for the diffent pulse models and BMP18 calculates
% the PD0, ln median and ln std for the pulse models. BMP 18 is used within
% ALLBMP18.

% Bray, J. D., Macedo, J. L.(2019). Procedure for Estimating Shear-Induced Seismic 
% Slope Displacement for Shallow Crustal Earthquakes  (ASCE).

% ky = yield coefficient
% Ts = fundamental period
% im = Sa(1.3Ts) spetral acceleration at the degraded system period
 %pulses
 Nd  = length(d);
 [PGVz,PGAz]=meshgrid(im(:,1),im(:,2));
  hd  = zeros(1,Nd);
 
 [pd0, mu,sig ]=ALLBMP18(Ts,ky,M,PGAz,PGVz,mode);%D100GM50, n=2
 
 for k = 1:Nd
    xhat    = (log(d(k))-mu)./sig;
    ccdf_a    = 0.5*(1-erf(xhat/sqrt(2)));
    cdf_a     = 1-ccdf_a;
    ccdf     = 1- (pd0+(1-pd0).*cdf_a);
    dlambda = ccdf.*MRD;%.*(ky<=PGAz);
    hd(k)   = sum(dlambda(:));
end
 


end

function [P, Lnd, stdv ]=ALLBMP18(ts_05,ky_05,Mag,saused_05,PGV,n)

%use n=1 for estimating DispD100 in terms of the IM 100
%use n=2 for estimating DispD100 in terms of the IM 50
%use n=3 for estimating DispD50  in terms of the IM 100
%use n=4 for estimating DispD50  in terms of the IM 50

%coefficientes for Ts<0.7 (PD=0)

A= [ 14.130,     10.787,    14.930,     14.930;
     9.223,      8.787,     10.383,     10.383;
    -1.928,      -1.660,    -1.971,     -1.971;
    -2.774,      -3.150,    -3.763,     -3.763;
    -7.67,       -7.560,    -8.812,     -8.812];
 %A =[A,A,A,A];
 
 B= [15.121,     12.771,    14.671,     14.671;
     10.023,     9.979,     10.489,     10.489;
    -2.443,      -2.286     -2.222,     -2.222;
     5.183,      4.965,     4.759,      4.759;
    -4.696,      -4.817,    -5.549,     -5.549];
 %B =[B,B,B,B];
 %FOR TS<=0.1 S
 C= [-7.415+0.229,      -6.462+0.227,       -8.615+0.224,       -7.718+0.221;
     -2.758,            -2.632,             -3.025,             -2.931;
     -0.280,            -0.278,             -0.313,             -0.319;
      0.529,             0.527,              0.571,              0.584;
      2.080,             1.978,              2.321,              2.261;
     -0.240,            -0.233,             -0.243,             -0.241;
      0.976-3.834,       1.069-3.813,        0.928-3.789,        1.031-3.762;
      0.000,             0.000,              0.000,              0.000;
      1.476,             1.547,              1.354,              1.458;
      0.000,            -0.06,               0.1314,             0.05];
 %FOR TS>0.1 S 
  D= [-7.415,           -6.462,             -8.615,             -7.718;
     -2.758,            -2.632,             -3.025,             -2.931;
     -0.280,            -0.278,             -0.313,             -0.319;
      0.529,             0.527,              0.571,              0.584;
      2.080,             1.978,              2.321,              2.261;
     -0.240,            -0.233,             -0.243,             -0.241;
      0.976,             1.069,              0.928,              1.031;
     -0.467,            -0.498,             -0.443,             -0.480;
      1.476,             1.547,              1.354,              1.458;
      0.000,            -0.06,               0.1314,             0.05]; 

STDV=[0.59,0.59,0.56,0.54];
 
                temp1=C(1,1);
                temp2=C(9,1);
                temp3=D(1,1);
                temp4=D(9,1);        

                temp5=C(1,2);
                temp6=C(9,2);
                temp7=D(1,2);
                temp8=D(9,2);        

                temp9=C(1,3);
                temp10=C(9,3);
                temp11=D(1,3);
                temp12=D(9,3);

                temp13=C(1,4);
                temp14=C(9,4);
                temp15=D(1,4);
                temp16=D(9,4);            
%%
nmodels=4;

%for n=1:1:nmodels
%     countcol4=0;
%     if n==1||n==3;
%        saused_05=saused100;
%        PGV=PGV100;
%     end
%     
%     if n==2||n==4;
%        saused_05=saused50;
%        PGV=PGV50;
%     end
           if PGV >=150
                C(1,1)=temp1+5.452;
                C(9, 1)=temp2-0.991;
                D(1,1)=temp3+5.452;
                D(9,1)=temp4-0.991;        

                C(1,2)=temp5+8.715;
                C(9,2)=temp6-1.644;
                D(1,2)=temp7+8.715;
                D(9,2)=temp8-1.644;        

                C(1,3)=temp9+6.776;
                C(9,3)=temp10-1.257;
                D(1,3)=temp11+6.776;
                D(9,3)=temp12-1.257;

                C(1,4)=temp13+7.349;
                C(9,4)=temp14-1.433;
                D(1,4)=temp15+7.349;
                D(9,4)=temp16-1.433;
            else
                C(1,1)=temp1;
                C(9,1)=temp2;
                D(1,1)=temp3;
                D(9,1)=temp4;        

                C(1,2)=temp5;
                C(9,2)=temp6;
                D(1,2)=temp7;
                D(9,2)=temp8;        

                C(1,3)=temp9;
                C(9,3)=temp10;
                D(1,3)=temp11;
                D(9,3)=temp12;

                C(1,4)=temp13;
                C(9,4)=temp14;
                D(1,4)=temp15;
                D(9,4)=temp16;                
            end
        
        
                a=A(:,n);
                b=B(:,n);
                c=C(:,n);
                d=D(:,n);
                stdm=STDV(n);
%%
            if n==1
                [P, Lnd, stdv ]=BMP18(ts_05,ky_05,Mag,saused_05,PGV,a,b,c,d,stdm);
            end
            
            if n==2
                [P, Lnd, stdv ]=BMP18(ts_05,ky_05,Mag,saused_05,PGV,a,b,c,d,stdm);
            end

            if n==3
                [P, Lnd, stdv ]=BMP18(ts_05,ky_05,Mag,saused_05,PGV,a,b,c,d,stdm);
            end

            if n==4
                [P, Lnd, stdv]=BMP18(ts_05,ky_05,Mag,saused_05,PGV,a,b,c,d,stdm);
            end
    
    
    
%end
end
%%
function [P, Lnd, stdv]=BMP18(ts,ky,Mag,sa,PGV,a,b,c,d,stdm)
      stdv=stdm;
      a1=a(1);
      a2=a(2);
      a3=a(3);
      a4=a(4);
      a5=a(5);
%       a6=a(6);
%       a7=a(7);
%       a8=a(8);
%       a9=a(9);
      
      b1=b(1);
      b2=b(2);
      b3=b(3);
      b4=b(4);
      b5=b(5);
%       b6=b(6);
%       b7=b(7);
%       b8=b(8);
%       b9=b(9);
            
      c1=c(1);
      c2=c(2);
      c3=c(3);
      c4=c(4);
      c5=c(5);
      c6=c(6);
      c7=c(7);
      c8=c(8);
      c9=c(9);
      c10=c(10);
      
      d1=d(1);
      d2=d(2);
      d3=d(3);
      d4=d(4);
      d5=d(5);
      d6=d(6);
      d7=d(7);
      d8=d(8);
      d9=d(9);
      d10=d(10);
      
      
      if ts<=0.1

      Disp4a=exp(c1+c2*log(ky)+c3*(log(ky))^2+c4*log(ky)*log(sa)+c5*log(sa)+c6*(log(sa))^2+c7*ts+c8*(ts)^2+c9*log(PGV)+c10*(Mag));
      P=1/(1+exp(-(a1+a2*log(ky)+a3*log(PGV)+a4*ts+a5*log(sa))));
      Lnd=log(Disp4a);
%             if P<0.5  
%                 Disp4=exp(log(Disp4a)+stdm*norminv(1-0.5/(1-P),0,1));
%             else
%                 Disp4=0.5;
%             end
%             
%             if P<0.84
%                 Disp4plus=exp(log(Disp4a)+stdm*norminv(1-0.16/(1-P),0,1));
%             else
%                 Disp4plus=0.5;
%             end
%             
%             if P<0.16
%                 Disp4minus=exp(log(Disp4a)+stdm*norminv(1-0.84/(1-P),0,1));
%             else
%                 Disp4minus=0.5;
%             end

      end
      
      if ts>0.1
            Disp4a=exp(d1+d2*log(ky)+d3*(log(ky))^2+d4*log(ky)*log(sa)+d5*log(sa)+d6*(log(sa))^2+d7*ts+d8*(ts)^2+d9*log(PGV)+d10*(Mag));
            Lnd=log(Disp4a);
        if ts<=0.7
           P=1/(1+exp(-(a1+a2*log(ky)+a3*log(PGV)+a4*ts+a5*log(sa))));
        end
        
        if ts>0.7
           P=1/(1+exp(-(b1+b2*log(ky)+b3*log(PGV)+b4*ts+b5*log(sa))));
        end
        
        if ts==0.7
          aa=1/(1+exp(-(a1+a2*log(ky)+a3*log(PGV)+a4*ts+a5*log(sa))));
          bb=1/(1+exp(-(b1+b2*log(ky)+b3*log(PGV)+b4*ts+b5*log(sa))));
          P=0.5*(aa+bb);
        end
        
%             if P<0.5  
%                 Disp4=exp(log(Disp4a)+stdm*norminv(1-0.5/(1-P),0,1));
%             else
%                 Disp4=0.5;
%             end
%             
%             if P<0.84
%                 Disp4plus=exp(log(Disp4a)+stdm*norminv(1-0.16/(1-P),0,1));
%             else
%                 Disp4plus=0.5;
%             end
%             
%             if P<0.16
%                 Disp4minus=exp(log(Disp4a)+stdm*norminv(1-0.84/(1-P),0,1));
%             else
%                 Disp4minus=0.5;
%             end
      end
      
end