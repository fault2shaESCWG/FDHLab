% this code calculates the ratio between DR and PF lengths
% OUTPUTS are in FIGURE ratio
clc
clear all
close all

%%
pathin = '../Analyse_SUREvers2_database/TABLE_db_20231026';
pathin2 = '../Analyse_SUREvers2_database/';
% output paths
path2 = fullfile('FIGURE','ratio');
if isempty(dir(path2))
mkdir(path2)
end

%%   user defined inputs
Across_strike_site_dim = [10,20,50,100,200,500]
step_slices = 10;% place slices every wstep meters
bufferfault = 5; %meters
maxdistance = 5000; 
nearfield = 200; % this is the near/far field threshold in meters
SoF_label = {'Reverse','Normal'};
SoF = 1; % Reverse
%SoF = 2; % Normal

if SoF ==1
event_id  = load(fullfile([pathin2,'list_Reverse.txt']));
input_data_rdist  = readtable(fullfile(pathin,'Reverse_r_distance_table.txt')); % input data

else SoF == 2
event_id  = load(fullfile([pathin2,'list_Normal.txt']));
input_data_rdist  = readtable(fullfile(pathin,'Normal_r_distance_table.txt')); % input data

end
id_all = [event_id(1:end,:)];

R_dist_table = [input_data_rdist];
%%

x_hw        = bufferfault:step_slices:(maxdistance+bufferfault);
x_fw        = -(maxdistance+bufferfault):step_slices:-bufferfault;      


    RATIO_HWsim = [];
    RATIO_FWsim = [];

    hwsim   = []; 
    fwsim   = [];  
    length_hwsim = [];length_fwsim = [];

R2_lengths_HW = [];
R2_lengths_FW = [];    
R1_length = [];    

for i = 1:size(id_all,1)
id = id_all(i,1);

input_data  = readtable(fullfile(pathin,[num2str(id),'_length_segments_table_20231018.txt'])); % input data
f_r2 = find(input_data.Rank==2);
R2_lengths = input_data(f_r2,:);
R2_lengths_HW = [R2_lengths_HW;R2_lengths(R2_lengths.HWFW==1,:)];
R2_lengths_FW = [R2_lengths_FW;R2_lengths(R2_lengths.HWFW==-1,:)];    
label(i,:) = [id];

f_r1 = find(input_data.Rank==1);
R1_length = [R1_length;input_data(f_r1,:)];
   LMAIN_fromtable(i,1) = id_all(i,3)*1000; LMAIN(i,1) = LMAIN_fromtable(i,1);
    LHW2(i,1)  = sum(R2_lengths_HW.length(R2_lengths_HW.IdE==id));
    LFW2(i,1)  = sum(R2_lengths_FW.length(R2_lengths_FW.IdE==id));

    
    HWsim = []; FWsim=[];
    HWsim   = R_dist_table(R_dist_table.IdE== id &  R_dist_table.rankR == 2 & R_dist_table.lengthR >0,:);
    FWsim   = R_dist_table(R_dist_table.IdE== id &  R_dist_table.rankR == 2 & R_dist_table.lengthR <0,:);

    
    hwsim(:,i)   = histcounts(HWsim.lengthR,x_hw);
    fwsim(:,i)   = histcounts(-1*FWsim.lengthR,x_hw); 

  
    length_hwsim(:,i) = hwsim(:,i)./sum(hwsim(:,i))*LHW2(i,1);
    length_fwsim(:,i) = fwsim(:,i)./sum(fwsim(:,i))*LFW2(i,1);
    RATIO_HWsim = [RATIO_HWsim,length_hwsim(:,i)/ LMAIN(i,1)];
    RATIO_FWsim = [RATIO_FWsim,length_fwsim(:,i)/ LMAIN(i,1)];
    

    
end

%% here we compute the ratio between DR-lengths and PF lenghts
% we compute the mean of the percentiles (50th, 84th, 97.5th) of this ratio
% for each earthquake

msize = 10; % dimension of points in the figure

i_wstep = 0; % to identify the column in the table, it has the same dimension of the length of Across_strike_site_dim vector
% assd = take the ith Across_strike_site_dim defined in the inputs
for assd = Across_strike_site_dim % spacing in meters of the site (number of slices)
i_wstep = i_wstep+1;
num_slices = round(assd/step_slices);

site_hw_R2 = [];
site_fw_R2 = [];
DdataHW = [];

 
    for i = 1:num_slices:(length(x_hw)-1)
    ending = min([i+num_slices-1,(length(x_hw)-1)]);
   DdataHW = [DdataHW;mean([x_hw(i), x_hw(ending)])];
  

    site_hw_R2 = [site_hw_R2;sum(RATIO_HWsim(i:ending,:),1)];
    site_fw_R2 = [site_fw_R2;sum(RATIO_FWsim(i:ending,:),1)];

 

   end

DdataFW = -1*DdataHW;



if assd > nearfield
    nearfield = nearfield + num_slices;
end
temphw_near = site_hw_R2(DdataHW <= nearfield,:);
tempfw_near = site_fw_R2(DdataFW >= -nearfield,:);
temphw_far = site_hw_R2(DdataHW > nearfield,:);
tempfw_far = site_fw_R2(DdataFW < -nearfield,:);

pcts_1HWsim(:,i_wstep) = nanmean(prctile(temphw_near,[50, 84, 97.5],1),2);
pcts_1FWsim(:,i_wstep) = nanmean(prctile(tempfw_near,[50, 84, 97.5],1),2);

pcts_2HWsim(:,i_wstep) = nanmean(prctile(temphw_far,[50, 84, 97.5],1),2);
pcts_2FWsim(:,i_wstep) = nanmean(prctile(tempfw_far,[50, 84, 97.5],1),2);




figure(i_wstep)
hold on
title(['SoF:',SoF_label{SoF},', Across strike S.D:',num2str(assd)])
for i = 1:size(id_all,1)
plot(DdataHW,site_hw_R2(:,i),'.','MarkerSize',msize,'color',[0.4 0.4 0.4])
plot(DdataFW,site_fw_R2(:,i),'.','MarkerSize',msize,'color',[0.4 0.4 0.4])
end
line([5 nearfield],[pcts_1HWsim(3,i_wstep) pcts_1HWsim(3,i_wstep)],'color','k');
line([-nearfield -5],[pcts_1FWsim(3,i_wstep) pcts_1FWsim(3,i_wstep)],'color','k');
line([nearfield 5005],[pcts_2HWsim(3,i_wstep) pcts_2HWsim(3,i_wstep)],'color','k');
line([-5005 -nearfield],[pcts_2FWsim(3,i_wstep) pcts_2FWsim(3,i_wstep)],'color','k');
grid on
xlim([-2000 2000])
ylim([0 1])
%legend(hleg)
set(gca,'YScale','log')
xlabel( 'Ratio' )
ylabel( 'Distance from the PF (m)' )
saveas(i_wstep,fullfile(path2,['SoF_',SoF_label{SoF},'_Ratio',num2str(assd),'meters_',date,'.pdf']),'pdf')

end
%%
writematrix([Across_strike_site_dim;pcts_1HWsim],['TABLE_outputs/',SoF_label{SoF},'RATIOpcts_nearfaultHWsim.txt'])
writematrix([Across_strike_site_dim;pcts_1FWsim],['TABLE_outputs/',SoF_label{SoF},'RATIOpcts_nearfaultFWsim.txt'])
writematrix([Across_strike_site_dim;pcts_2HWsim],['TABLE_outputs/',SoF_label{SoF},'RATIOpcts_farfaultHWsim.txt'])
writematrix([Across_strike_site_dim;pcts_2FWsim],['TABLE_outputs/',SoF_label{SoF},'RATIOpcts_farfaultFWsim.txt'])
%%
T_paper = [pcts_1HWsim(3,:);pcts_2HWsim(3,:);pcts_1FWsim(3,:);pcts_2FWsim(3,:)];
%% R2_histogram_lengths 
R2_histogram_length = [];
all_hwfw =[R2_lengths_HW;R2_lengths_FW];
R2_histogram_length = [prctile(R2_lengths_HW.length, [16 84]); 
                        prctile(R2_lengths_FW.length, [16 84]);
                        prctile(all_hwfw.length, [16 84])];
R2_histogram_length = array2table(R2_histogram_length);
R2_histogram_length.Properties.VariableNames = {'16thPCT','84thPCT'};
R2_histogram_length.Properties.RowNames = {'Hangingwall','Footwall','all'} ;        
%%
writetable(R2_histogram_length,['TABLE_outputs/',SoF_label{SoF},'_R2_histogram_length.txt'],'WriteRowNames',true);