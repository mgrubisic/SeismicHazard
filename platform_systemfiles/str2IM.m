function[IM]=str2IM(str)

NIM = length(str);
IM  = zeros(NIM,1);
for i=1:NIM
    switch str{i}
        case 'PGD', IM(i)=-10;
        case 'PGV', IM(i)=-1;
        case 'PGA', IM(i)=0;
        case 'Ia',  IM(i)=-5;
        otherwise
            stri = lower(str{i});
            stri = strrep(stri,' ','');
            stri = strrep(stri,'sa','');
            stri = strrep(stri,'t','');
            stri = strrep(stri,'s','');
            stri = strrep(stri,'o','');
            stri = strrep(stri,'=','');
            stri = strrep(stri,'(' ,'');
            stri = strrep(stri,')' ,'');
            IM(i)=str2double(stri);
    end
end