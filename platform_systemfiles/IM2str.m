function[str]=IM2str(T)

if isnumeric(T)
    str = cell(1,length(T));
    for i=1:length(T)
        val = T(i);
        if val==-10
            str{i} = 'PGD';
        elseif val==-1
            str{i} = 'PGV';
        elseif val==0
            str{i} = 'PGA';
        elseif val==-3
            str{i} = 'Duration (s)';
        elseif val==-5
            str{i} = 'Ia';            
        elseif val>0
            str{i} = ['Sa(T=',num2str(val),')'];
        else
            str{i} = ['T=',num2str(val)];
        end
    end
end

if iscell(T)
    str = T;
end