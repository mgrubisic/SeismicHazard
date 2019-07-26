function[dt,acc]=read_record(filename)

[dt,acc] = readSIBERRISK(filename); 
acc = acc/9.81;