clc

% fun = 'SeismicHazard.m';
% % fun = 'GMselect.m';
% % fun = 'PSDA2.m';
% % fun = 'pshatoolbox_methods.m';
% list  = matlab.codetools.requiredFilesAndProducts(fun)';
% for i=1:size(list,1)
% copyfile(list{i},'C:\Users\Gabriel\Desktop\GMSelect')
% end

%% Pcode
clc
% cd('C:\Users\gcandia\Desktop\temp')
D=dir('*.m');
[~,B]=intersect({D.name},'Listdependencies.m');
D(B)=[];
for i=1:size(D,1)
    i
%     delete(D(i).name);
    %pcode(D(i).name);
end