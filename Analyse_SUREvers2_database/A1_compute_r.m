% this code calculates the r-distance between DR and PF
% it uses the shapefile in the DATABASE with data from SURE 2.0 
% user has to set: (i) id_list; (ii)maxdiff that is the resampling spacing;
% (iii) buffer_dist_R2_R15


% OUTPUTS are tables in TABLE defined in pathout1
%%
 clc
 clear all
 close all
%% user input
event_rev  = load('list_Reverse.txt');
event_nor  = load('list_Normal.txt');
id_all = [event_rev(:,1);event_nor(:,1)]

maxdiff = 10; % in meters
buffer_dist_R2_R15 = 1000; % in meters
%% output paths
pathout1 = 'TABLE_db_20231026';

if isempty(dir(pathout1))
mkdir(pathout1)
end

 %%
    dati_rupture_all = [];
    dati_rupture_all = shaperead(fullfile('SURE-main','SURE2.0_ruptures','SURE2.0_ruptures.shp'));
%%
    
        if isnumeric([dati_rupture_all.Comp_rank]) == 0
        for dr = 1:size(dati_rupture_all,1)
        dati_rupture_all(dr).Comp_Rank = str2num(dati_rupture_all(dr).Comp_rank);
        end
        end
        
        
   dati_rupture_allIdE = [dati_rupture_all.IdE]';
%%   
 for i = 1:size(id_all,1)
     
     T = [];
    id = num2str(id_all(i,1));
     dati_rupture = [];
     ind_id = find(dati_rupture_allIdE == str2num(id));
     dati_rupture = dati_rupture_all(ind_id,:);        
        
%% prepare variables
fault_res_x = []; fault_res_y = [];
fault_res_id = [];fault_res_Rank = [];
fault2_res_x = [];fault2_res_y = [];
fault2_res_id = [];fault2_res_Rank = [];
fault2_res_hwfwPF = []; fault2_res_hwfwNEAR = []; fault2_res_linked = [];
fault15_res_x = [];fault15_res_y = [];
fault15_res_id = [];fault15_res_Rank = [];
fault15_res_hwfwPF = [];fault15_res_linked = [];
fault3_res_x = [];fault3_res_y = [];
fault3_res_id = [];fault3_res_Rank = [];
fault3_res_hwfwPF = [];fault3_res_linked = [];

FINAL_Length = [];
e = referenceEllipsoid('earth');
tempout = []; tempLength=[];
%% start calculation

prog = 0; % this counts in sequence all the  point
R = [dati_rupture.Comp_rank]';
    
    %% values in fields of the shape are converted in vector
    
    L = [dati_rupture.Length]';
    id_feature = [dati_rupture.IdS]';
    HW_FW_PF = [dati_rupture.HW_FW_PF]';
    %% check that fields in the shapefiles are in the right format
    ispresent = isfield(dati_rupture,'IdS_PF')
    if ispresent ==1
    id_link  = [dati_rupture.IdS_PF]';
    elseif ispresent ==0
        id_link (1:(size(dati_rupture)),1)=0;
    end
    %%
    HWFW_PF(1:length(R),1)=0;
    HWFW_PF_input=[{dati_rupture.HW_FW_PF}]';
    HWind = find((strcmp(HWFW_PF_input,'HW')==1));
    FWind = find((strcmp(HWFW_PF_input,'FW')==1));
    HWFW_PF(HWind,1) = 1;
    HWFW_PF(FWind,1) =-1;
    
    HWFW_near(1:length(R),1)=0;
    HWFW_near_input=[{dati_rupture.HW_FW_near}]';
    HWind = find((strcmp(HWFW_near_input,'HW')==1));
    FWind = find((strcmp(HWFW_near_input,'FW')==1));
    HWFW_near(HWind,1) = 1;
    HWFW_near(FWind,1) =-1;
    
 
   %% saves the length, converts to UTM and resamples with maxdiff, returns to WGS
   disp(['eq ID:',id])
   disp('code is saving the length of DRs, converting to UTM and resampling with maxdiff')
   p1 = 1; p2 = 1; p3 = 1; p4 = 1;
  for j = 1:size(dati_rupture,1) % 
   size(dati_rupture,1)-j ;
      
    tempLength=[tempLength;str2num(id),j,R(j), L(j),HWFW_PF(j),id_feature(j)];
    tempVert=[];
    tempVert=[(dati_rupture(j).X)',(dati_rupture(j).Y)'];
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
    
    if R(j)==1
        fault_res(p1).x = xres;
        fault_res(p1).y = yres;
        fault_res(p1).id = repmat(id_feature(j),length(xres),1);
        fault_res(p1).Rank = repmat(R(j),length(xres),1);
        fault_res_x = [fault_res_x;fault_res(p1).x];
        fault_res_y = [fault_res_y;fault_res(p1).y];
        fault_res_id = [fault_res_id;fault_res(p1).id];
        fault_res_Rank = [fault_res_Rank; fault_res(p1).Rank];
       p1=p1+1; 
     
    elseif R(j)==1.5
        fault15_res(p2).x =xres;
        fault15_res(p2).y =yres;
        fault15_res(p2).id =repmat(id_feature(j),length(xres),1);
        fault15_res(p2).Rank = repmat(R(j),length(xres),1);
        fault15_res(p2).hwfwPF = repmat(HWFW_PF(j),length(xres),1);
        fault15_res(p2).linked = repmat(id_link(j),length(xres),1);
        fault15_res_x = [fault15_res_x; fault15_res(p2).x];
        fault15_res_y = [fault15_res_y; fault15_res(p2).y];
        fault15_res_id = [fault15_res_id; fault15_res(p2).id];
        fault15_res_Rank = [fault15_res_Rank; fault15_res(p2).Rank];
        fault15_res_hwfwPF = [fault15_res_hwfwPF; fault15_res(p2).hwfwPF];
        fault15_res_linked = [fault15_res_linked; fault15_res(p2).linked];
       p2 = p2+1;
        
    elseif R(j)==2
        fault2_res(p3).x =xres;
        fault2_res(p3).y =yres;
        fault2_res(p3).id =id_feature(j);
        fault2_res(p3).Rank = repmat(R(j),length(xres),1);
        fault2_res(p3).hwfwPF = repmat(HWFW_PF(j),length(xres),1);
        fault2_res(p3).hwfwNEAR = repmat(HWFW_near(j),length(xres),1);
        fault2_res(p3).linked =repmat(id_link(j),length(xres),1);
        fault2_res_x = [fault2_res_x; fault2_res(p3).x];
        fault2_res_y = [fault2_res_y; fault2_res(p3).y];
        fault2_res_id = [fault2_res_id; fault2_res(p3).id];
        fault2_res_Rank = [fault2_res_Rank; fault2_res(p3).Rank];
        fault2_res_hwfwPF = [fault2_res_hwfwPF; fault2_res(p3).hwfwPF];
        fault2_res_hwfwNEAR = [fault2_res_hwfwNEAR; fault2_res(p3).hwfwNEAR];
        fault2_res_linked = [fault2_res_linked; fault2_res(p3).linked];
       p3 = p3+1;
        
    elseif R(j)>=3
        fault3_res(p4).x =xres;
        fault3_res(p4).y =yres;
        fault3_res(p4).id =id_feature(j);
        fault3_res(p4).Rank = repmat(R(j),length(xres),1);
        fault3_res(p4).hwfwPF = repmat(HWFW_PF(j),length(xres),1);
        fault3_res(p4).linked =repmat(id_link(j),length(xres),1);
        fault3_res_x = [fault3_res_x; fault3_res(p4).x];
        fault3_res_y = [fault3_res_y; fault3_res(p4).y];
        fault3_res_id = [fault3_res_id; fault3_res(p4).id];
        fault3_res_Rank = [fault3_res_Rank; fault3_res(p4).Rank];
        fault3_res_hwfwPF = [fault3_res_hwfwPF; fault3_res(p4).hwfwPF];
        fault3_res_linked = [fault3_res_linked; fault3_res(p4).linked];
       p4 = p4+1;
       
    end 

  end
FINAL_Length=[FINAL_Length; tempLength];
T2 = array2table(FINAL_Length,'VariableNames',{'IdE', 'IdR', 'Rank', 'length', 'HWFW', 'id_feature'});
writetable(T2,fullfile(pathout1,strcat(id,'_length_segments_table_20231018.txt')));


  
%%
disp('code is calculating distances of each resampled point')
disp('first loop for rank 2 DR')
 for j = 1:length(fault2_res_x) 
      prog = prog+1;

%% CHECK: if the DR segment is associated to a particular PF Rank 1 or a Rank 1.5 fault  	
  if fault2_res_linked(j) > 0  % distance must be computed from a particular feature 

      tempY = []; tempX = [];   
      pos_id_feature = find(fault_res_id ==  fault2_res_linked(j));
      rango = 1;
      tempY = fault_res_y(pos_id_feature) ;
      tempX = fault_res_x(pos_id_feature) ;
     
      if  isempty(pos_id_feature)
      pos_id_feature = find(fault15_res_id ==  fault2_res_linked(j));
      rango = 1.5;
      tempY = fault15_res_y(pos_id_feature) ;
      tempX = fault15_res_x(pos_id_feature) ;
      end

        dist = [];
        dist = distance(fault2_res_y(j),fault2_res_x(j),tempY,tempX,e);
        pos_mindist = find(dist == min(dist),1,'first') ;
        attributes_mindist = [];
        attributes_mindist = [fault2_res_y(j),fault2_res_x(j), tempY(pos_mindist),tempX(pos_mindist),fault2_res_Rank(j),min(dist)];
        if isnan(fault2_res_hwfwNEAR(j))
        tempout = [tempout;str2num(id), prog, attributes_mindist(1:end-1), attributes_mindist(end)*fault2_res_hwfwPF(j),rango];
        elseif ~isnan(fault2_res_hwfwNEAR(j))
         tempout = [tempout;str2num(id), prog, attributes_mindist(1:end-1), attributes_mindist(end)*fault2_res_hwfwNEAR(j),rango];    
        end
  else  % if the DR segment is not associated to a particular PF Rank 1 or a Rank 1.5 fault  	
      
  
    dist1 = [];
    dist1 = distance(fault2_res_y(j),fault2_res_x(j), fault_res_y,fault_res_x,e);
    
    temp_pos_mindist1 = find(dist1 == min(dist1),1,'first') ;
    temp_attributes_mindist1 = [fault2_res_y(j),fault2_res_x(j),...
                                 fault_res_y(temp_pos_mindist1),fault_res_x(temp_pos_mindist1),fault2_res_Rank(j),min(dist1)];
   
    % if R1.5 exists then compute distances also from them.
    % we save the minimum distance between  R2-R1 and R2-R1.5
    % if the "distance R2- R1.5  > buffer" then the DR R2 is associated with the R1
    if ~isempty(fault15_res_y)
  
    dist15 = [];
    dist15 = distance(fault2_res_y(j),fault2_res_x(j),fault15_res_y,fault15_res_x,e);
    
    temp_pos_mindist15 = find(dist15 == min(dist15),1,'first') ;
    temp_attributes_mindist15 = [fault2_res_y(j),fault2_res_x(j),...
                                 fault15_res_y(temp_pos_mindist15),fault15_res_x(temp_pos_mindist15),fault2_res_Rank(j),min(dist15)];
    
        if min(dist15) > buffer_dist_R2_R15  % se la distanza da R1.5 e' > buffer, viene scartata
        temp_attributes_mindist15 = [];
        end
    
    else
    temp_attributes_mindist15 = [];
    end
    
    if isempty(temp_attributes_mindist15) | (temp_attributes_mindist1(end) <= temp_attributes_mindist15(end))
 
    tempout = [tempout;str2num(id), prog, temp_attributes_mindist1(1:end-1), temp_attributes_mindist1(end)*fault2_res_hwfwPF(j),1];
    
    elseif (temp_attributes_mindist1(end) > temp_attributes_mindist15(end))
  
        tempout = [tempout;str2num(id), prog, temp_attributes_mindist15(1:end-1), temp_attributes_mindist15(end)*fault2_res_hwfwNEAR(j),1.5];
    end
    
    end
 end
 %% 
 disp('code is calculating distances of each resampled point')
 disp('second loop for rank 1.5')
 for j = 1:length(fault15_res_x) 
      prog = prog+1;
        length(fault15_res_x) -j ;
  	
  if fault15_res_linked(j) >0 % CHECK: if the R1.5 is associated to a particular segment R1 

      tempY = []; tempX = [];   
      pos_id_feature = find(fault_res_id ==  fault15_res_linked(j));
      rango = 1;
    
     tempY = fault_res_y(pos_id_feature) ;
     tempX = fault_res_x(pos_id_feature) ;
     

  dist = [];
  dist = distance(fault15_res_y(j),fault15_res_x(j),...
                  tempY,tempX,e);
  pos_mindist = find(dist == min(dist),1,'first') ;
  attributes_mindist = [];
  attributes_mindist = [fault15_res_y(j),fault15_res_x(j), tempY(pos_mindist),tempX(pos_mindist),fault15_res_Rank(j),min(dist)];
    
       
   tempout = [tempout;str2num(id), prog, attributes_mindist(1:end-1), attributes_mindist(end)*fault15_res_hwfwPF(j),rango];

       
    
  else  
   

    dist1 = [];
    dist1 = distance(fault15_res_y(j),fault15_res_x(j),...
        fault_res_y,fault_res_x,e);
    
    temp_pos_mindist1 = find(dist1 == min(dist1),1,'first') ;
    temp_attributes_mindist1 = [fault15_res_y(j),fault15_res_x(j),...
                                 fault_res_y(temp_pos_mindist1),fault_res_x(temp_pos_mindist1),fault15_res_Rank(j),min(dist1)];
   

    tempout = [tempout;str2num(id), prog, temp_attributes_mindist1(1:end-1), temp_attributes_mindist1(end)*fault15_res_hwfwPF(j),1];
    

    
  end
 end
 %% 
disp('code is calculating distances of each resampled point')
 disp('third loop for rank 3') 
 for j = 1:length(fault3_res_x) 
      prog = prog+1;
        length(fault3_res_x) -j ;
  	
  if fault3_res_linked(j) >0  % CHECK: if the R3 is associated to a particular PF Rank 1

      tempY = []; tempX = [];   
      pos_id_feature = find(fault_res_id ==  fault3_res_linked(j));
      rango = 1;
    
     tempY = fault_res_y(pos_id_feature) ;
     tempX = fault_res_x(pos_id_feature) ;
     
  dist = [];
  dist = distance(fault3_res_y(j),fault3_res_x(j),...
                  tempY,tempX,e);
  pos_mindist = find(dist == min(dist),1,'first') ;
  attributes_mindist = [];
  attributes_mindist = [fault3_res_y(j),fault3_res_x(j), tempY(pos_mindist),tempX(pos_mindist),fault3_res_Rank(j),min(dist)];
    
       
  tempout = [tempout;str2num(id), prog, attributes_mindist(1:end-1), attributes_mindist(end)*fault3_res_hwfwPF(j),rango];

       
    
  else  
   

    dist1 = [];
    dist1 = distance(fault3_res_y(j),fault3_res_x(j),...
        fault_res_y,fault_res_x,e);
    
    temp_pos_mindist1 = find(dist1 == min(dist1),1,'first') ;
    temp_attributes_mindist1 = [fault3_res_y(j),fault3_res_x(j),...
                                 fault_res_y(temp_pos_mindist1),fault_res_x(temp_pos_mindist1),fault3_res_Rank(j),min(dist1)];
   
    tempout = [tempout;str2num(id), prog, temp_attributes_mindist1(1:end-1), temp_attributes_mindist1(end)*fault3_res_hwfwPF(j),1];     

  end
    
    end
T=array2table(tempout,...
    'VariableNames',{'IdE', 'IdR', 'latR', 'lonR',...
    'latPr', 'lonPr', 'rankR', 'lengthR', 'rankPr'});
writetable(T,fullfile(pathout1,strcat(id,'_r_distance_table_20231018.txt')));



end
%%

