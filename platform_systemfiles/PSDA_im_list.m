function[IM,imethod]=PSDA_im_list(handles,fld,T2,T3)

methods  = pshatoolbox_methods(5);
s=vertcat(handles.(fld).source);s=s(:);
mechs = unique({s.mechanism});
[~,B] = intersect({'interface','intraslab','crustal','grid'},mechs);
func = {methods.str}';

Smodels = T3(:,B);
Smodels = unique(Smodels(:));
usedM   = cell(size(Smodels));

for i=1:length(usedM)
    [~,B]=intersect({handles.sys.SMLIB.id},Smodels{i});
    usedM{i}=handles.sys.SMLIB(B).str;
end

[~,b]=intersect(func,usedM);
IMfactor = zeros(0,1);
for i=1:length(b)
    IMfactor = [IMfactor;methods(b(i)).Safactor(:)];
end
IMfactor = unique(IMfactor);

Tnat = zeros(size(T2,1),1);
for i=1:length(Tnat)
    switch ischar(T2{i,2})
        case 0,Tnat(i)=T2{i,2};
        case 1,dd = regexp(T2{i,2},'\, ','split'); Tnat(i)=str2double(dd{1});
    end
end

Tnat     = unique(Tnat);
IM  = [];
for i=1:length(IMfactor)
    if IMfactor(i)<0
        IM = [IM;IMfactor(i)]; %#ok<*AGROW>
    else
        IM = [IM;IMfactor(i)*Tnat]; %#ok<*AGROW>
    end
end

imethod = zeros(1,length(usedM));
for i=1:length(imethod)
    [~,B]=intersect(func,usedM{i});
    imethod(i)=methods(B).integrator;
end

imethod = unique(imethod);

