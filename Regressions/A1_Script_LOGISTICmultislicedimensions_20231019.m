% this code calculates the parametres of the logistic regressions

%% OUTPUTS are table in TABLE 
clear all
clc
close all
output_i_wstep_BC1 = [];
output_i_wstep_BC2 = [];
output_i_wstep_BC3 = [];

%%
% USER DEFINED INPUTS
pathin = '../Analyse_SUREvers2_database/TABLE_db_20231026';
pathin2 = '../Analyse_SUREvers2_database/';

for kin = 1:2
    if kin ==1
Kin = 'Normal';
elseif kin ==2
Kin = 'Reverse';
    end

sizes = [10,20,50,100,200,500];
bufferfault = 5; % we discard the first meters  from the main fault
maxdistance = 20000;

% output paths
pathout1 ='TABLE_outputs';
if isempty(dir(pathout1))
mkdir(pathout1)
end

input_data  = readtable(fullfile(pathin,[Kin,'_r_distance_table.txt'])); % input data
event_mw_l  = load(fullfile([pathin2,'list_',Kin,'.txt']));
i_wstep = 0;

%%
disp(['calculating logistic for:',Kin])
%%
    
  
x_bin        = bufferfault:10:(maxdistance+bufferfault); % edges of the classes
x_centre     = x_bin + 10/2; % centre of the classes of distance for the plot
outputname = strcat(Kin,'_',date); % used in the output name

%% 
    count_eq =size(event_mw_l,1)
    % inizialize variables
    hwall   = [];fwall   = []; 
    hw_R2_R1   = [];fw_R2_R1   = [];
    hw_R2_R15   = [];fw_R2_R15   = [];
    hw_Rcomp_R1  = [];fw_Rcomp_R1  = [];
    MWdata   = [];IDdata   = [];Ddata   = [];
    Xdata = []; flagHW = [];flagFW=[];
%%  selection of the data from the input table
    % we have many categories: 
    %HW vs FW  ALL, rank2 vs rank1,rank1.5 vs rank1, TRIGGERED (rank =1.5,21,22,3) vs Rank1
    HWall = input_data(input_data.rankR > 1 & input_data.lengthR >= bufferfault,:);
    FWall = input_data(input_data.rankR > 1 & input_data.lengthR <= -1*bufferfault,:);
    FWall.lengthR=abs(FWall.lengthR);  % positive distances, before the conversion the footwall values were negative
    
    HW_R2_R1data   = HWall(HWall.rankR == 2 & HWall.rankPr == 1,:);
    FW_R2_R1data   = FWall(FWall.rankR == 2 & FWall.rankPr == 1,:);
    HW_R2_R15data   = HWall(HWall.rankR == 2 & HWall.rankPr == 1.5,:);
    FW_R2_R15data   = FWall(FWall.rankR == 2 & FWall.rankPr == 1.5,:);
    
    HWcompdata  = HWall(HWall.rankR == 1.5 |...
                    HWall.rankR == 3 |...
                    HWall.rankR == 21 |...
                    HWall.rankR == 22,:);
                
                
    FWcompdata  = FWall(FWall.rankR == 1.5 |...
                    FWall.rankR == 3 |...
                    FWall.rankR == 21 |...
                    FWall.rankR == 22,:);
 %% counts the ruptures per bin - creates a vector column for each
  % category and for ID, MW and the distance (Ddata)
    for i =1: count_eq
 
    IDdata     = [IDdata;repmat(event_mw_l(i,1),length(x_centre)-1,1)];
    MWdata     = [MWdata;repmat(event_mw_l(i,2),length(x_centre)-1,1)];

    hwall      = [hwall;histcounts(HWall.lengthR(HWall.IdE==event_mw_l(i,1)),x_bin)'];
    fwall      = [fwall;histcounts(FWall.lengthR(FWall.IdE==event_mw_l(i,1)),x_bin)'];
    
    hw_R2_R1      = [hw_R2_R1;histcounts(HW_R2_R1data.lengthR(HW_R2_R1data.IdE==event_mw_l(i,1)),x_bin)'];
    fw_R2_R1      = [fw_R2_R1;histcounts(FW_R2_R1data.lengthR(FW_R2_R1data.IdE==event_mw_l(i,1)),x_bin)']; 

    hw_R2_R15      = [hw_R2_R15;histcounts(HW_R2_R15data.lengthR(HW_R2_R15data.IdE==event_mw_l(i,1)),x_bin)'];
    fw_R2_R15      = [fw_R2_R15;histcounts(FW_R2_R15data.lengthR(FW_R2_R15data.IdE==event_mw_l(i,1)),x_bin)']; 

    hw_Rcomp_R1      = [hw_Rcomp_R1;histcounts(HWcompdata.lengthR(HWcompdata.IdE==event_mw_l(i,1)),x_bin)'];
    fw_Rcomp_R1      = [fw_Rcomp_R1;histcounts(FWcompdata.lengthR(FWcompdata.IdE==event_mw_l(i,1)),x_bin)']; 

    end 
    clear i
%%
% if there's at least 1 rupture per bin = 1, otherwise = 0
    hwall(hwall>0)=1;
    fwall(fwall>0)=1;
    hw_R2_R1(hw_R2_R1>0)=1;
    fw_R2_R1(fw_R2_R1>0)=1;
    hw_R2_R15(hw_R2_R15>0)=1;
    fw_R2_R15(fw_R2_R15>0)=1;
    hw_Rcomp_R1(hw_Rcomp_R1>0)=1;
    fw_Rcomp_R1(fw_Rcomp_R1>0)=1;
%%    
for wstep = sizes% spacing in meters of the site (number of slices)
i_wstep = i_wstep+1;
num_slices = round(wstep/10);

site_hw_R2_R1 = [];
site_fw_R2_R1 = [];
site_hw_R2_R15 = [];
site_fw_R2_R15 = [];
site_hw_Rcomp_R1 = [];
site_fw_Rcomp_R1 = [];
Ddata = [];

flagHW = [];
flagFW = [];


   for j =1: count_eq
 
    t = find(IDdata == event_mw_l(j,1));

% for cycle for summing occurrences in the number of consecutive slices

for i = 1:(length(x_bin)-1)
    ending = (t(i)-i) + min([i+num_slices-1,(length(x_bin)-1)]);
    ending_ddata =  min([i+num_slices-1,(length(x_bin)-1)]);
    site_hw_R2_R1 = [site_hw_R2_R1;sum(hw_R2_R1(t(i):ending))];
    site_fw_R2_R1 = [site_fw_R2_R1;sum(fw_R2_R1(t(i):ending))];
    site_hw_R2_R15 = [site_hw_R2_R15;sum(hw_R2_R15(t(i):ending))];
    site_fw_R2_R15 = [site_fw_R2_R15;sum(fw_R2_R15(t(i):ending))];
    site_hw_Rcomp_R1 = [site_hw_Rcomp_R1;sum(hw_Rcomp_R1(t(i):ending))];
    site_fw_Rcomp_R1 = [site_fw_Rcomp_R1;sum(fw_Rcomp_R1(t(i):ending))];

    Ddata = [Ddata;mean([x_bin(i), x_bin(ending_ddata)])];
end
   end
   
    site_hw_R2_R1(site_hw_R2_R1>0)=1;
    site_fw_R2_R1(site_fw_R2_R1>0)=1;
    site_hw_R2_R15(site_hw_R2_R15>0)=1;
    site_fw_R2_R15(site_fw_R2_R15>0)=1;
    site_hw_Rcomp_R1(site_hw_Rcomp_R1>0)=1;
    site_fw_Rcomp_R1(site_fw_Rcomp_R1>0)=1;


flagHW(1:(size(site_hw_R2_R1,1)),1) = 0;
flagFW(1:(size(site_fw_R2_R1,1)),1) = 1;
%%
Xdata = [MWdata,(Ddata),flagHW;
    MWdata,(Ddata),flagFW];
Xtable = array2table(Xdata);

%%
YdataC1 = [hw_R2_R1;fw_R2_R1];
YdataC2_alldist = [hw_R2_R15;fw_R2_R15];
YdataC3 = [hw_Rcomp_R1;fw_Rcomp_R1];

X = ([Xtable.Xdata1 Xtable.Xdata2 double(Xtable.Xdata3) ]);
X_C2 = Xdata(Xdata(:,2)<=2000,:);
YdataC2 = YdataC2_alldist(Xdata(:,2)<=2000,:);
%% logistic fit
disp('logistic fit')
[BC1,devC1,statsC1] = mnrfit(X,categorical(YdataC1))
[BC2,devC2,statsC2] = mnrfit(X_C2,categorical(YdataC2))
[BC3,devC3,statsC3] = mnrfit(X,categorical(YdataC3))

tempC1 = [wstep,NaN;BC1,statsC1.p];
tempC2 = [wstep,NaN;BC2,statsC2.p];
tempC3 = [wstep,NaN;BC3,statsC3.p];

output_i_wstep_BC1 = [output_i_wstep_BC1,tempC1];
output_i_wstep_BC2 = [output_i_wstep_BC2,tempC2];
output_i_wstep_BC3 = [output_i_wstep_BC3,tempC3];

end


save(strcat(pathout1,'/parameters_logistic_multisizeC1_',Kin,'.txt'),'output_i_wstep_BC1','-ascii')
save(strcat(pathout1,'/parameters_logistic_multisizeC2_',Kin,'.txt'),'output_i_wstep_BC2','-ascii')
save(strcat(pathout1,'/parameters_logistic_multisizeC3_',Kin,'.txt'),'output_i_wstep_BC3','-ascii')

end
