function[M,dPm]=getMdPm(deagg)

M    = unique(deagg(:,1));
dPm  = zeros(size(M));
rate = deagg(:,3);
for i=1:length(M)
    dPm(i) =sum(rate(deagg(:,1)==M(i)));
end
dPm = dPm/sum(dPm);