function [dt,acc1,num]=readSIBERRISK(filename)

if exist(filename,'file')==0
    disp(['Filename ',filename,' not found on matlab''s path']);
    dt=[];
    acc1=[];
    return
end
fid  = fopen(filename);
line = fgetl(fid);
line = line(32:end-1);
if nargout==3
    num  = datevec(line);
end
line = fgetl(fid);
line = regexp(line,'\s+','split');
dt   = 1/str2double(line{5});
acc1 = textscan(fid,'%f','Headerlines',4);
acc1 =acc1{1}';
fclose(fid);