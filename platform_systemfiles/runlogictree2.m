function[lambda,deagg]=runlogictree2(handles)

opt            = handles.opt;
model          = handles.model;
site           = cell2mat(handles.h.p(:,2:4));
Vs30           = handles.h.Vs30;
weights        = handles.sys.WEIGHT(:,4);
sitelist       = handles.site_selection;

%% variable initialization
IM        = opt.IM;
im        = opt.im;
Nsites    = size(site,1);
Nim       = size(im,1);
NIM       = length(IM);
Nbranch   = length(model);

Nsource   = 0;
for i=1:Nbranch
    Nsource = max(Nsource,length(model(i).source));
end

%% do not run analysis if ind is empty
lambda  = nan(Nsites,Nim,NIM,Nsource,Nbranch);
deagg   = cell(size(lambda));

%% run logic tree
fprintf('\n');
spat  = '%-20s   |   %-30s   |   %-20s     Runtime:  %-4.3f s\n';
t0 = tic;
fprintf('                               SEISMIC HAZARD + (M,R) DEAGGREGATION \n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');
for i=1:Nbranch
    ti=tic;
    ID1=model(i).id1; if length(ID1)>27,ID1=[ID1(1:27),'...'];end
    ID2=model(i).id2; if length(ID2)>27,ID2=[ID2(1:27),'...'];end
    ID3=model(i).id3; if length(ID3)>27,ID3=[ID3(1:27),'...'];end
    if weights(i)~=0
        deagg(:,:,:,:,i)= runhazard2(im,IM,site,Vs30,opt,model(i),Nsource,sitelist);
    end
    fprintf(spat,ID1,ID2,ID3,toc(ti));
end

for i=1:numel(deagg)
    if ~isempty(deagg{i})
        lambda(i)=sum(deagg{i}(:,3));
    end
end
fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

