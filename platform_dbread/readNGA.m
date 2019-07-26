function[dt,acc]=readNGA(filename)

if exist(filename,'file')==0
    disp(['Filename ',filename,' not found on matlab''s path']);
    dt=[];
    acc=[];
    return
end

fid       = fopen(filename);
line3     = textscan(fid,'%f%f%s%s%s',1,'Headerlines',3);
acc       = textscan(fid,'%f','Headerlines',1);acc=acc{1}';
dt   = line3{2};
fclose(fid);
