% this codes prepare a shapefile of r-distances
clear all
clc
%%
pathtable = 'TABLE_db_20231026';
pathSHP ='SHAPEFILES'
if exist('pathSHP','dir') ==0
mkdir (pathSHP)
end

for kin = 1:2
    if kin ==1
Kin = 'Normal';
elseif kin ==2
Kin = 'Reverse';
    end
%%
T = readtable(fullfile(pathtable,strcat(Kin,'_r_distance_table.txt')));
%%
for i =1:size(T,1)

    Data(i).X = [T.lonR(i),T.lonPr(i)];
    Data(i).Y = [T.latR(i),T.latPr(i)];
    Data(i).IdE = [T.IdE(i)];
    Data(i).IdR = [T.IdR(i)];
    Data(i).Distance = [T.lengthR(i)];
    Data(i).RangoDR = [T.rankR(i)];
    Data(i).RangoPF = [T.rankPr(i)];

end   
%%

[Data(1:size(T,1)).Geometry] = deal('Line');
%%
shapewrite(Data, fullfile(pathSHP,strcat(Kin,'_mappa_r.shp')));

end