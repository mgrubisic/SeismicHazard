function dsha_bin(FileName,PathName,handles,etype)

switch etype
    case 'IS'
        FILE = fullfile(PathName,FileName);
        [rate,Y] =dsha_is(handles);
        fid = fopen(FILE,'w');
        fwrite(fid,[rate,Y],'double');
        fclose(fid);
        fprintf('done\n')
        
    case 'KM'
        FILE = fullfile(PathName,FileName);
        if isempty(handles.krate) && isempty(handles.kY)
            [rate,Y] =dsha_kmeans(handles,handles.optkm);
        else
            rate = handles.krate;
            Y    = handles.kY;
        end
        fid = fopen(FILE,'w');
        fwrite(fid,[rate,Y],'double');
        fclose(fid);
        fprintf('done\n')
end