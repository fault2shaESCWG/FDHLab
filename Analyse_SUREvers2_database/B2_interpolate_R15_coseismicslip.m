% this code interpolates coseismic throw along the R1.5

clear all
clc
%%
pathoutFig = 'FigureThrowInterpolated';
if isempty(dir(pathoutFig))
mkdir(pathoutFig)
end
pathoutTable = 'TABLE_db_20231026';
if isempty(dir(pathoutTable))
mkdir(pathoutTable)
end

R = 1.5;
maxdiff = 10; % in meters
dati_point_all = readtable(fullfile('SURE-main','SURE2.0_Slip_Obs_matlab.xlsx'),'format','auto');
IdE = unique(dati_point_all.IdE);
    dati_rupture_all = [];
    dati_rupture_all = shaperead(fullfile('SURE-main','SURE2.0_ruptures','SURE2.0_ruptures.shp'));


%%
  if isnumeric([dati_rupture_all.Comp_rank]) == 0
        for dr = 1:size(dati_rupture_all,1)
        dati_rupture_all(dr).Comp_rank = str2num(dati_rupture_all(dr).Comp_rank);
        end
        end
        %%
        
   dati_rupture_allIdE = [dati_rupture_all.IdE]';
%% 

%%
out = [];
for id = 1:length(IdE)
    %%
% looks for measurements points associated with R1.5

rows_point = find(dati_point_all.IdE == IdE(id) & dati_point_all.Comp_rank==R);
%%
dati_point = [dati_point_all.Longitude(rows_point), dati_point_all.Latitude(rows_point), dati_point_all.T(rows_point)];
%%
% assign the value contained in the column "SH" of the database to points
% with no-value of Throw
    nv_nsub = find(isnan(dati_point_all.T) & (dati_point_all.SH>0));
    for nsub = 1:length(nv_nsub)
    dati_point_all.T(nv_nsub(nsub)) = dati_point_all.SH(nv_nsub(nsub));
    end
% remove points with no value in T    
    nv = find(~isnan(dati_point_all.T));
    dati_points = dati_point_all(nv,:);
dati_point(isnan(dati_point(:,3)),:)=[];
%
if ~isempty(dati_point) %if measurement points exist then calculate interpolation
%%

     dati_rupture = [];
     ind_id = find(dati_rupture_allIdE == IdE(id));
     dati_rupture = dati_rupture_all(ind_id,:);   

%% 
rows = [];dati = [];
rows_point = [];
tips_all = [];
dist =[];tips_dist =[];
dati_point_mod = [];

xq = [];yq = [];F = [];vq =[];
x = [];y = [];v =[];
%%
rows = find([dati_rupture.Comp_rank]' == R);
dati = dati_rupture(rows,:);


%%
% assign Throw = 0 m to the tips of each R.15 segment (note this is different from B1 for the PF)

for i = 1: size(dati,1)
tips_all =   [tips_all; dati(i).X(1),dati(i).Y(1);dati(i).X(end-1),dati(i).Y(end-1)] ;
end
% add points with T =0 to the tips of the R1.5
dati_point_mod = [dati_point;tips_all,repmat(0,size(tips_all,1),1)];
%%
% throw is interpolated at verteces of the shapefile
xq = [];
yq = [];  
vect_IdS = [];
for i = 1: size(dati,1)
    
    IdStemp = dati(i).IdS;
    tempVert=[];
    tempVert=[(dati(i).X)',(dati(i).Y)'];
    tempVert(isnan(tempVert(:,1)),:) = [];
    [xutm,yutm,datum] = ll2utm(tempVert(:,2),tempVert(:,1)); 
    if length(datum) > 1
    % in case of traces cross 2 datum
     yresutm = [];xresutm = [];datum_temp=[];
    for ld = 1: length(datum)
        
    [yresutm_1, xresutm_1] = interpm(yutm(datum==datum(ld)),xutm(datum==datum(ld)),maxdiff);
    datum_temp = [datum_temp;repmat(datum(ld),length(yresutm_1),1)];
    
    yresutm = [yresutm;yresutm_1]; xresutm = [xresutm;xresutm_1];
    end
     [yres,xres]=utm2ll(xresutm, yresutm,datum_temp);
    else
        [yresutm, xresutm] = interpm(yutm,xutm,maxdiff);
        [yres,xres]=utm2ll(xresutm, yresutm,datum);
    end
vect_IdS = [vect_IdS;repmat(IdStemp,length(xres),1)];  
xq = [xq;xres];
yq = [yq;yres]; 
end

x = dati_point_mod(:,1);
y = dati_point_mod(:,2);
v = dati_point_mod(:,3);
%%
% Use scatteredInterpolant to perform interpolation on a 2-D or 3-D data set of scattered data. 
% scatteredInterpolant returns the interpolant F for the given data set. You can evaluate F at a set of query points, such as (xq,yq) in 2-D, to produce interpolated values vq = F(xq,yq).

F = scatteredInterpolant(x,y,v);
F.Method = 'linear';
%%
vq = F(xq,yq);
vq(vq<0)=0;
%%
figure(id)
hold on
plot([dati.X]',[dati.Y]');

scatter(xq,yq,30,vq,'filled')
scatter(x,y,50,v,'Marker','s')
text(x,y,num2str(v),'FontSize',6)
colorbar
saveas(id,fullfile(pathoutFig,strcat(num2str(IdE(id)),'_rank15.png')),'png')
%%
close(id)

%%

out = [out;xq,yq,round(abs(vq),2),repmat(IdE(id),length(xq),1),vect_IdS,repmat(R,length(xq),1)];
   
%%
Tout = array2table(out);
Tout.Properties.VariableNames = {'lon','lat','Throw','IdE','IdS','Rank'};


writetable(Tout,fullfile(pathoutTable,strcat(num2str(IdE(id)),'_R15_vd_interp.txt')));%%%%

end
end
