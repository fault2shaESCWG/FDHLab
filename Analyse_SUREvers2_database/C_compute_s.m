% this code calculates the s-distance between DR and PF, R1.5
% it uses the excel SURE file database in the folder DATABASE

clear all
close all
%% 
for kin = 1:2
    if kin ==1
        event  = load('list_Reverse.txt'); name = 'Reverse';
    elseif kin == 2; 
        event  = load('list_Normal.txt');name = 'Normal';
    end
%%
buffer_from_R15 = 1000;
e=referenceEllipsoid('earth');
% output paths
pathoutTable = 'TABLE_db_20231026';
if isempty(dir(pathoutTable))
mkdir(pathoutTable)
end

%% read the database of points
dati_point_all = readtable(fullfile('SURE-main','SURE2.0_Slip_Obs_matlab.xlsx'),'format','auto');

% Change names to numbers
dati_point_all.HW_FW_PF = categorical(dati_point_all.HW_FW_PF);
dati_point_all.HW_FW_PF = renamecats(dati_point_all.HW_FW_PF,{'HW','FW'},{'1','-1'});
dati_point_all.HW_FW_PF = str2double(string(dati_point_all.HW_FW_PF));
dati_point_all.HW_FW_near = categorical(dati_point_all.HW_FW_near);
dati_point_all.HW_FW_near = renamecats(dati_point_all.HW_FW_near,{'HW','FW'},{'1','-1'});
dati_point_all.HW_FW_near = str2double(string(dati_point_all.HW_FW_near));
%%   
% assign SH to points with no-value in T
    nv_nsub = find(isnan(dati_point_all.T) & (dati_point_all.SH>0));
    for nsub = 1:length(nv_nsub)
    dati_point_all.T(nv_nsub(nsub)) = dati_point_all.SH(nv_nsub(nsub));
    end
% remove points with no value in T or T ==0    
    nv = find(~isnan(dati_point_all.T) & abs((dati_point_all.T) > 0));
    dati_points = dati_point_all(nv,:);

%%
% select PF and DR       
    indR1       = find(dati_points.Comp_rank == 1); %ranking 1==principal fault otherwhise off-fault
    indDR       = find(dati_points.Comp_rank > 1); % DR rupture ranks
    R1points_all    =   dati_points(indR1,:);
    DRpoints_all   =   dati_points(indDR,:);
%%
out = [];
id_list = event(:,:)
%%
    dati_rupture_all = [];
    dati_rupture_all = shaperead(fullfile('SURE-main','SURE2.0_ruptures','SURE2.0_ruptures.shp'));

    
        if isnumeric([dati_rupture_all.Comp_rank]) == 0
        for dr = 1:size(dati_rupture_all,1)
        dati_rupture_all(dr).Comp_Rank = str2num(dati_rupture_all(dr).Comp_rank);
        end
        end
dati_rupture_allIdE = [dati_rupture_all.IdE]';
%%

for i = 1:size(id_list,1)
%%
id       =  id_list(i,1)
Mw_event =  id_list(i,2)
%%
% select measurements points from the SURE excel corresponding to the
% earthquake ID and classify for Rank
% R1.5,21,22 e 3 have distances from R1
% R2 can have disatnces from R1 or R1.5
puntiDR = DRpoints_all(DRpoints_all.IdE == id,:);
puntiDRcomp = puntiDR(puntiDR.Comp_rank ~= 2,:);
puntiDRsimple = puntiDR(puntiDR.Comp_rank == 2,:);

%% shapefile of the event
ind_id = find(dati_rupture_allIdE == id);
dati_rupture = dati_rupture_all(ind_id,:); 
  R = [dati_rupture.Comp_rank]';
 id_feature = [dati_rupture.IdS]';
  ispresent = isfield(dati_rupture,'IdS_PF')
    if ispresent ==1
    id_link  = [dati_rupture.IdS_PF]';
    elseif ispresent ==0
        id_link (1:(size(dati_rupture)),1)=0;
    end
%%
% open resampled points R1 e R15
inputnameR1=fullfile(pathoutTable,strcat(num2str(id),'_R1_vd_interp.txt'));
inputnameR15=fullfile(pathoutTable,strcat(num2str(id),'_R15_vd_interp.txt'));
if exist(inputnameR1,'file') 
R1_vd_pointsTable = readtable(inputnameR1);
R1_vd_pointsTable.Throw(isnan(R1_vd_pointsTable.Throw))=0;
if exist(inputnameR15,'file') 
R15_vd_pointsTable = readtable(inputnameR15);
R15_vd_pointsTable.Throw(isnan(R15_vd_pointsTable.Throw))=0;

else
    R15_vd_pointsTable = [];
end
 %%   puntiDRcomp
    if isempty(puntiDRcomp)
        disp(strcat('no points Rank1.5,21,22,3 for eq id:',num2str(id)))
    else
        
        for j =1:size(puntiDRcomp,1) %% start of the loop for secondary points in the excels
 %%
        f_dist =[]; az = []; points_for_virtualD = []; point_id = [];
        dist = []; min_dist = [];  values = []; ordine = [];
     
        IdS = puntiDRcomp.IdS(j);
    % check if the segment to which the point belongs to, is linked to a particular feature      
    h_link =   id_link(id_feature == IdS);
if h_link > 0
    R1_vd_points = R1_vd_pointsTable(R1_vd_pointsTable.IdS == h_link,:);
    
else 
    R1_vd_points = R1_vd_pointsTable;
   
end
% looks for the point that minimize the distance
temp = distance(puntiDRcomp.Latitude(j),puntiDRcomp.Longitude(j),R1_vd_points.lat,R1_vd_points.lon,e);
dist = min(temp);
h_dist = find(temp == dist);
coord_point = [R1_vd_points.lat(h_dist),R1_vd_points.lon(h_dist)];
VD_PF =  R1_vd_points.Throw(h_dist);
% looks for the points of R1 that are inside a given searching- radius
raggio = dist/2;
temp = [];
temp = distance(coord_point(1,1),coord_point(1,2),R1_vd_pointsTable.lat,R1_vd_pointsTable.lon,e);
h_dist_incircle = find(temp <= raggio);
meanVD_PF =  mean(R1_vd_pointsTable.Throw(h_dist_incircle));
        if isnan(puntiDRcomp.HW_FW_PF(j))
        HWFW = puntiDRcomp.HW_FW_near(j); 
        else
        HWFW = puntiDRcomp.HW_FW_PF(j); 
        end 
          

  % save the outputs
out = [out;id, Mw_event, puntiDRcomp.IdO(j),...
        puntiDRcomp.Latitude(j),puntiDRcomp.Longitude(j),...
        coord_point(1,:),...
        puntiDRcomp.Comp_rank(j),1,puntiDRcomp.T(j),mean(VD_PF),meanVD_PF,...
        dist,HWFW,kin];
        end
    end
%%  loop for R2
    if isempty(puntiDRsimple)
        disp(strcat('no points Rank 2 for eq id:',num2str(id)))
    else
        
        for j =1:size(puntiDRsimple,1) %% start of the loop for secondary points in the excels
 %%   
        IdS = puntiDRsimple.IdS(j);

        %%
       h_link =   id_link(id_feature == IdS);
if h_link > 0
    R1_vd_points = R1_vd_pointsTable(R1_vd_pointsTable.IdS == h_link,:);
   
    if isempty(R1_vd_points) & ~isempty(R15_vd_pointsTable)
    R15_vd_points = R15_vd_pointsTable(R15_vd_pointsTable.IdS == h_link,:);
    
    end
    
    if isnan(puntiDRsimple.HW_FW_near(j))
          HWFW = puntiDRsimple.HW_FW_PF(j); 
    else
          HWFW = puntiDRsimple.HW_FW_near(j); 
    end
    
else 
    R1_vd_points = R1_vd_pointsTable;
    R15_vd_points = R15_vd_pointsTable;
end
% looks for the point that minimize the distance
temp1 = distance(puntiDRsimple.Latitude(j),puntiDRsimple.Longitude(j),R1_vd_points.lat,R1_vd_points.lon,e);
dist1 = min(temp1);
if ~isempty(R15_vd_points)
temp2 = distance(puntiDRsimple.Latitude(j),puntiDRsimple.Longitude(j),R15_vd_points.lat,R15_vd_points.lon,e);
dist2 = min(temp2);
else
  dist2 = NaN;
end
 
disp([dist1,dist2])

if (dist1 <= dist2) | (dist2 > buffer_from_R15) | isnan(dist2)
    disp('Rank 2 associated to Rank 1')
h_dist = find(temp1 == dist1);
coord_point = [R1_vd_points.lat(h_dist),R1_vd_points.lon(h_dist)];
VD_PF =  R1_vd_points.Throw(h_dist);   
% looks for the points of R1 that are inside a given searching- radius
raggio = dist1/2;
temp = [];
temp = distance(coord_point(1,1),coord_point(1,2),R1_vd_pointsTable.lat,R1_vd_pointsTable.lon,e);
h_dist_incircle = find(temp <= raggio);
meanVD_PF =  mean(R1_vd_pointsTable.Throw(h_dist_incircle));
HWFW = puntiDRsimple.HW_FW_PF(j);

    
  % save the outputs
out = [out;id, Mw_event, puntiDRsimple.IdO(j),...
        puntiDRsimple.Latitude(j),puntiDRsimple.Longitude(j),...
        coord_point(1,:),...
        puntiDRsimple.Comp_rank(j),1,puntiDRsimple.T(j),mean(VD_PF),meanVD_PF,...
        dist1,HWFW,kin];
elseif (dist1 > dist2) & (dist2 <= buffer_from_R15)
    disp('Rank 2 associated to Rank 1.5')
 
h_dist = find(temp2 == dist2);
coord_point = [R15_vd_points.lat(h_dist),R15_vd_points.lon(h_dist)];
VD_PF =  R15_vd_points.Throw(h_dist);   
% looks for the points of R1 that are inside a given searching- radius
raggio = dist2/2;
temp = [];
temp = distance(coord_point(1,1),coord_point(1,2),R15_vd_pointsTable.lat,R15_vd_pointsTable.lon,e);
h_dist_incircle = find(temp <= raggio);
meanVD_PF =  mean(R15_vd_pointsTable.Throw(h_dist_incircle));
   if isnan(puntiDRsimple.HW_FW_near(j))
          HWFW = puntiDRsimple.HW_FW_PF(j); 
   else
          HWFW = puntiDRsimple.HW_FW_near(j);
   end
  % save the outputs
out = [out;id, Mw_event, puntiDRsimple.IdO(j),...
        puntiDRsimple.Latitude(j),puntiDRsimple.Longitude(j),...
        coord_point(1,:),...
        puntiDRsimple.Comp_rank(j),1.5,puntiDRsimple.T(j),mean(VD_PF),meanVD_PF,...
        dist2,HWFW,kin];
end
        end
    end
end
end

Tout = array2table(out);
Tout.Properties.VariableNames = {'IdE','Mw','IdO','latDR','lonDR','latPF','lonPF','RankDR','RankPF','ThrowDR','ThrowPFpoint','ThrowPFmean','distance','HWFW','kinR=1kinN=2'};
writetable(Tout,fullfile(pathoutTable,strcat(name,'_s_distance.txt')));%%%%

end


