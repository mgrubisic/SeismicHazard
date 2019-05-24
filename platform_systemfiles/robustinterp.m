function [y]=robustinterp(X,Y,x,style)

ind = isnan(X+Y);
X(ind)=[];
Y(ind)=[];

[Xu,ind]= unique(X,'stable');
Yu      = Y(ind);

switch style
    case 'loglog'
        y = exp(interp1(log(Xu),log(Yu),log(x)));
    case 'linear'
        y = interp1(Xu,Yu,x);
end
