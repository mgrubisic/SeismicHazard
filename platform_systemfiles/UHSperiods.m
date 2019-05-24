function tlist=UHSperiods(handles)

Tmax = [];
for i=1:length(handles.model)
    for j=1:length(handles.model(i).source)
        str = func2str(handles.model(i).source(j).gmpe.handle);
        [~,T]=mGMPE_info(str);
        Tmax = [Tmax;max(T)]; %#ok<AGROW>
    end
end
Tmax  = min(Tmax);
Tcond = str2double(handles.Cond_Period.String);
tbase = logsp(0.01,10,30)';
tlist = unique([tbase;Tcond;Tmax]);
tlist(tlist>Tmax)=[];
tlist = tlist(:);