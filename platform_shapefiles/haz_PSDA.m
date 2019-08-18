function[haz]=haz_PSDA(handles)
handles.site_selection = 1:size(handles.h.p,1);
opt = handles.opt;

% finds what integration methods are declared in the logic tree
T3       = handles.T3(:,2:5);
[imethod,usedM]  = getmethod(handles,T3);
haz=struct('imstandard',[],'IMstandard',[],'lambda',[],'deagg',[],'imvector',[],'IMvector',[],'corrlist',[],'MRD',[]);

if any(ismember(imethod,[1 2]))
    handles.opt = opt_update(handles,usedM,opt,1);
    haz.imstandard = handles.opt.im;
    haz.IMstandard = handles.opt.IM;
    [haz.lambda,haz.deagg]=runlogictree2(handles);
end

if any(ismember(imethod,[3 4])) % Ellen´s rigid and flexible slopes
    SMLIB=handles.sys.SMLIB;
    handles.opt = opt_update(handles,usedM,opt,3);

    corrlist =[];
    T3models = unique(T3(:));
    for i=1:length(T3models)
        [~,B]=intersect({SMLIB.id},T3models{i});
        switch SMLIB(B).str
            case 'psda_RathjeF'
                corrlist =[corrlist;SMLIB(B).param.rho]; %#ok<*AGROW>
            case 'psda_RathjeR'
                corrlist =[corrlist;SMLIB(B).param.rho];
        end
    end
    sourcelist = cell(0,1);
    [~,B1]=intersect({SMLIB.id},handles.T3(:,2)); D1=intersect({SMLIB(B1).str},{'psda_RathjeF','psda_RathjeR'});
    [~,B2]=intersect({SMLIB.id},handles.T3(:,3)); D2=intersect({SMLIB(B2).str},{'psda_RathjeF','psda_RathjeR'});
    [~,B3]=intersect({SMLIB.id},handles.T3(:,4)); D3=intersect({SMLIB(B3).str},{'psda_RathjeF','psda_RathjeR'});
    [~,B4]=intersect({SMLIB.id},handles.T3(:,5)); D4=intersect({SMLIB(B4).str},{'psda_RathjeF','psda_RathjeR'});
    if ~isempty(D1), sourcelist=[sourcelist,'interface'];end
    if ~isempty(D2), sourcelist=[sourcelist,'instraslab'];end
    if ~isempty(D3), sourcelist=[sourcelist,'crustal','shallowcrustal'];end
    if ~isempty(D4), sourcelist=[sourcelist,'grid'];end
    
    haz.imvector = handles.opt.im;
    haz.IMvector = handles.opt.IM;
    haz.corrlist=corrlist;
    for ii=1:length(corrlist)
        [haz.MRD(:,:,:,:,:,ii)]=runlogictree2V(handles,corrlist(1),sourcelist);
    end
end

function[imethod,usedM]=getmethod(handles,T3)

methods  = pshatoolbox_methods(5);
s=vertcat(handles.model.source);s=s(:);
mechs = unique({s.mechanism});
mechs = strrep(mechs,'shallowcrustal','crustal');
[~,B] = intersect({'interface','intraslab','crustal','grid'},mechs);
func = {methods.str}';

Smodels = T3(:,B);
Smodels = unique(Smodels(:));
usedM   = cell(size(Smodels));

for i=1:length(usedM)
    [~,B]=intersect({handles.sys.SMLIB.id},Smodels{i});
    usedM{i}=handles.sys.SMLIB(B).str;
end

imethod = zeros(1,length(usedM));
for i=1:length(imethod)
    [~,B]=intersect(func,usedM{i});
    imethod(i)=methods(B).integrator;
end

imethod = unique(imethod);

function[opt]=opt_update(handles,usedM,opt,mtype)
T2       = handles.T2;
methods  = pshatoolbox_methods(5);

switch mtype
    case 1, B = ismember(usedM,{'psda_BMT2017M','psda_BT2007M','psda_J07M','psda_RS09M','psda_BM2019M','psda_BT2007','psda_AM1988','psda_J07'}); usedM = usedM(B);
    case 3, B = ismember(usedM,{'psda_RathjeF','psda_RathjeR'}); usedM = usedM(B);
end

func = {methods.str}';
[~,b]=intersect(func,usedM);
IMfactor = zeros(0,1);
for i=1:length(b)
    IMfactor = [IMfactor;methods(b(i)).Safactor(:)];
end
IMfactor = unique(IMfactor);
Tnat     = unique(cell2mat(T2(:,2)));
IM  = [];
for i=1:length(IMfactor)
    if IMfactor(i)<=0
        IM = [IM;IMfactor(i)]; %#ok<*AGROW>
    else
        IM = [IM;IMfactor(i)*Tnat]; %#ok<*AGROW>
    end
end

im0=opt.im;
[r,c]= size(im0);
NIM  = length(IM);
im = zeros(r,NIM);

if c==1
    im = repmat(im0,1,length(IM));
    if any(IM==-1)
        im(:,IM==-1)=logsp(0.1,60,r)';
    end
    opt.im=im;
    opt.IM=IM;
    return;
end

for i=1:NIM
    if IM(i)<0
        ind = IM0==IM(i);
        im(:,i)=im0(:,ind);
    else
        disc = abs(IM0-IM(i));
        [~,ind] = min(disc);
        im(:,i)=im0(:,ind);
    end
end
opt.im=im;
opt.IM=IM;


