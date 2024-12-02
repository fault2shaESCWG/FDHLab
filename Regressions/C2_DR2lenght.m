% this code calculates the ratio between DR and PF lengths
% OUTPUTS are in FIGURE ratio
clc
clear all
close all
%%
outputname = 'R2_histogram_length_NoSoF';
%%
pathin = '../Analyse_SUREvers2_database/TABLE_db_20231026';
pathin2 = '../Analyse_SUREvers2_database/';
% output paths
path2 = fullfile('FIGURE','R2DRlenghts/');
if isempty(dir(path2))
mkdir(path2)
end

SoF_label = {'Reverse','Normal'};
%SoF = 1; % Reverse
SoF = 2; % Normal

if SoF ==1
event_id  = load(fullfile([pathin2,'list_Reverse.txt']));
input_data_rdist = readtable(fullfile(pathin,'Reverse_r_distance_table.txt')); % input data

else SoF == 2
event_id  = load(fullfile([pathin2,'list_Normal.txt']));
input_data_rdist  = readtable(fullfile(pathin,'Normal_r_distance_table.txt')); % input data

end
id_all = [event_id(1:end,:)];
%%  count percentage on the total table of R2-R1 and R2-R1.5
hr2 = find(input_data_rdist.rankR ==2);
hr2r1 = find(input_data_rdist.rankR ==2 & input_data_rdist.rankPr==1);
hr2r15 = find(input_data_rdist.rankR ==2 & input_data_rdist.rankPr>1);
perc_r2r1 = (length(hr2r1)/length(hr2))*100
perc_r2r15 = (length(hr2r15)/length(hr2))*100
%%
R2_lengths_HW = [];
R2_lengths_FW = [];    

for i = 1:size(id_all,1)
id = id_all(i,1);

input_data  = readtable(fullfile(pathin,[num2str(id),'_length_segments_table_20231018.txt'])); % input data
f_r2 = find(input_data.Rank==2);
R2_lengths = input_data(f_r2,:);
R2_lengths_HW = [R2_lengths_HW;R2_lengths(R2_lengths.HWFW==1,:)];
R2_lengths_FW = [R2_lengths_FW;R2_lengths(R2_lengths.HWFW==-1,:)];    
label(i,:) = [id];
   
end

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
[param_hw] = lognfit(R2_lengths_HW.length);
[param_fw] = lognfit(R2_lengths_FW.length);
[param] = lognfit(all_hwfw.length);
LOGN_MU_SIGMA = [param_hw(1),param_hw(2); param_fw(1),param_fw(2); param(1), param(2)];
LOGN_MU_SIGMA = array2table(LOGN_MU_SIGMA);
LOGN_MU_SIGMA.Properties.VariableNames = {'mu','sigma'};
LOGN_MU_SIGMA.Properties.RowNames = {'Hangingwall','Footwall','all'} ;        


%%
figure(1)
subplot(1,2,1)
hold on
title([SoF_label{SoF},'- hangingwall'])
histogram(R2_lengths_HW.length,'BinWidth',10,'Normalization','pdf')
line([R2_histogram_length{1,1} R2_histogram_length{1,1}],[0 0.02],'color','r')
line([R2_histogram_length{1,2} R2_histogram_length{1,2}],[0 0.02],'color','r')
xlabel('length (m)')
ylabel ('pdf')

subplot(1,2,2)
hold on
title([SoF_label{SoF},'- footwall'])
histogram(R2_lengths_FW.length,'BinWidth',10,'Normalization','pdf')
line([R2_histogram_length{2,1} R2_histogram_length{2,1}],[0 0.02],'color','r')
line([R2_histogram_length{2,2} R2_histogram_length{2,2}],[0 0.02],'color','r')
xlabel('length (m)')
ylabel ('pdf')

saveas(1,[path2,SoF_label{SoF},'_R2R1_histogram_length.png'],'png');
%%
writetable(R2_histogram_length,[path2,SoF_label{SoF},'_R2R1_histogram_length.txt'],'WriteRowNames',true);
writetable(LOGN_MU_SIGMA,[path2,SoF_label{SoF},'LOGN_MU_SIGMA.txt'],'WriteRowNames',true);