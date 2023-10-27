% this code calculates the ratio between DR and PF lengths
% OUTPUTS are in FIGURE ratio
clc
clear all
close all

%%
pathin = '../Analyse_SUREvers2_database/TABLE_db_20231026';
% output paths
path2 = fullfile('FIGURE','ratio');
if isempty(dir(path2))
mkdir(path2)
end

%%   user defined inputs
site_dim = [10,20,50,100,200,500]
wstep = 10;
bufferfault = 5; %metri
maxdistance = 5000;  
event_rev  = load('list_Reverse.txt');
event_nor  = load('list_Normal.txt');
id_all = [event_rev(1:end,:);event_nor(:,:)];
%%    
input_data_rdistnorm  = readtable(fullfile('Table_inputs','Normal_r_distance_table.txt')); % input data
input_data_rdistrev  = readtable(fullfile('Table_inputs','Reverse_r_distance_table.txt')); % input data
R_dist_table = [input_data_rdistnorm;input_data_rdistrev];
%%

x_hw        = bufferfault:wstep:(maxdistance+bufferfault);
x_fw        = -(maxdistance+bufferfault):wstep:-bufferfault;      

                
    tabella = [];

    PROPORZIONE_HWsim = [];
    PROPORZIONE_FWsim = [];

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

    
   %LMAIN(i,1) = sum(R1_length.length(R1_length.IdE==id));
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
    PROPORZIONE_HWsim = [PROPORZIONE_HWsim,length_hwsim(:,i)/ LMAIN(i,1)];
    PROPORZIONE_FWsim = [PROPORZIONE_FWsim,length_fwsim(:,i)/ LMAIN(i,1)];
    

    
end

%%
msize = 10;
col = colormap('jet');
pos_col = round(linspace(1,size(col,1),size(id_all,1)));
col_sof = col(pos_col,:);
i_wstep = 0;
for wstep = site_dim% spacing in meters of the site (number of slices)
i_wstep = i_wstep+1;
num_slices = round(wstep/10);

site_hw_R2 = [];
site_fw_R2 = [];
DdataHW = [];

 
    for i = 1:num_slices:(length(x_hw)-1)
    ending = min([i+num_slices-1,(length(x_hw)-1)]);
   DdataHW = [DdataHW;mean([x_hw(i), x_hw(ending)])];
  

    site_hw_R2 = [site_hw_R2;sum(PROPORZIONE_HWsim(i:ending,:),1)];
    site_fw_R2 = [site_fw_R2;sum(PROPORZIONE_FWsim(i:ending,:),1)];

 

   end

DdataFW = -1*DdataHW;


nearfield = 200;
if wstep > nearfield
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
for i = 1:size(id_all,1)
plot(DdataHW,site_hw_R2(:,i),'.','MarkerSize',msize,'color',[0.4 0.4 0.4])
plot(DdataFW,site_fw_R2(:,i),'.','MarkerSize',msize,'color',[0.4 0.4 0.4])
% in case you want to color points using ID uncomment and use the following
% two lines instead of the previous two
%hleg(i) = plot(DdataHW,site_hw_R2(:,i),'.','MarkerSize',msize,'color',col_sof(i,:),'Display',num2str(id_all(i,1)));
%plot(DdataFW,site_fw_R2(:,i),'.','MarkerSize',msize,'color',col_sof(i,:))
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
saveas(i_wstep,fullfile(path2,['Ratio',num2str(site_dim(i_wstep)),'meters_',date,'.pdf']),'pdf')

end
%%
writematrix([site_dim;pcts_1HWsim],'TABLE_outputs/RATIOpcts_nearfaultHWsim.txt')
writematrix([site_dim;pcts_1FWsim],'TABLE_outputs/RATIOpcts_nearfaultFWsim.txt')
writematrix([site_dim;pcts_2HWsim],'TABLE_outputs/RATIOpcts_farfaultHWsim.txt')
writematrix([site_dim;pcts_2FWsim],'TABLE_outputs/RATIOpcts_farfaultFWsim.txt')

