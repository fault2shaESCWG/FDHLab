% questo codice prepara lo shapefile delle distanze r
clear all
clc
%%
pathtable = 'TABLE_db_20231026';
pathSHP ='SHAPEFILES'
if exist('pathSHP','dir') ==0
mkdir (pathSHP)
end
kin = 'Reverse';
%kin = 'Normal';
%%
T = readtable(fullfile(pathtable,strcat(kin,'_s_distance.txt')));

% 'IdE','Mw','IdO','latDR','lonDR','latPF','lonPF','RankDR','RankPF','ThrowDR','ThrowPFpoint','ThrowPFmean','distance','HWFW','kinR=1kinN=2'
%%
for i =1:size(T,1)

    Data(i).X = [T.lonDR(i),T.lonPF(i)];
    Data(i).Y = [T.latDR(i),T.latPF(i)];
    Data(i).IdE = [T.IdE(i)];
    Data(i).IdO = [T.IdO(i)];
    Data(i).ThrowDR = [T.ThrowDR(i)];
    Data(i).ThrowPFpoint = [T.ThrowPFpoint(i)];
    Data(i).ThrowPFmean = [T.ThrowPFmean(i)];
    Data(i).distance = [T.distance(i)];
    Data(i).HWFW = [T.HWFW(i)];

end   
%%

[Data(1:size(T,1)).Geometry] = deal('Line');
%%
shapewrite(Data, fullfile(pathSHP,strcat(kin,'_mappa_s.shp')));

%