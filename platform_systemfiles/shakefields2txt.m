function[]=shakefields2txt(FileName,PathName,handles)

Nshake = length(handles.shakefield);
str    = handles.pop_scenario.String;

ButtonName = questdlg('Sort Data by...', ...
    'Export', ...
    'Period', 'Site', 'Period');

for i=1:Nshake
    FILE   = [strrep(FileName,'.out',num2str(i)),'.out'];
    fname  = fullfile(PathName,FILE);
    header = str{i};
    writeshakefield(fname,header,handles.shakefield(i),handles.opt,handles.h,ButtonName);
end

function writeshakefield(fname,header,shakefield,opt,h,ButtonName)
Y      = exp(bsxfun(@plus,shakefield.mulogIM,shakefield.L*shakefield.Z));
IMs    = [opt.IM1;opt.IM2];
NIM    = length(IMs);
Nsites = size(h.p,1);

fid = fopen(fname,'w');
fprintf(fid,'Ground Shaking Field: %s\n',header);
spat = 'Per, Site, ';
for i=1:opt.NumSim
    spat = [spat,sprintf('R%-5g, ',i)]; %#ok<AGROW>
end
spat = [spat(1:end-2),'\n'];
fprintf(fid,spat);

switch ButtonName
    case 'Period'
        for i=1:NIM
            for j=1:Nsites
                k= j+Nsites*(i-1);
                spat = ['%3.2f, %3g, ',repmat('%5.4f, ',1,opt.NumSim),'\n'];
                spat = strrep(spat,', \n','\n');
                fprintf(fid,spat,IMs(i),j,Y(k,:));
            end
        end
    case 'Site'
        for j=1:Nsites
            for i=1:NIM
                k= j+Nsites*(i-1);
                spat = ['%3.2f, %3g, ',repmat('%5.4f, ',1,opt.NumSim),'\n'];
                spat = strrep(spat,', \n','\n');
                fprintf(fid,spat,IMs(i),j,Y(k,:));
            end
        end
end

fclose(fid);
if ispc,winopen(fname);end
