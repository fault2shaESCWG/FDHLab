% this code calculates the P of DR rank 2 occurrence at s ite using a montecarlo approach 
clc
clear all
close all
%%
rng('default');
% output paths 
path1 = 'TABLE_outputs';
path2 = fullfile('FIGURE','simulations');
if isempty(dir(path1))
mkdir(path1)
end
if isempty(dir(path2))
mkdir(path2)
end
%%
% USER DEFINED INPUTS
% sampling of legths 1 for uniform 2 for lognormal
SL = 1;
% length of the fault (in meter)
%Fault_length = 5800;
Fault_length = 34000;
%dr_over_tip_distance = 0;%(in meter) >0 to place DR over the tips
along_strike_dimension = 100; % meters
site_distance = 7600;% meters from the PF
%SoF = ['Reverse'];
SoF = ['Normal'];
HWFW = 'HW'; % HW = Hanging wall; FW = Footwall location of the site
Simulations = 10^4 ;% number of simulations
Simulations_figure = 100;% number of simulations placed in the figure
%%
% code
%SpaceDRLength = dr_over_tip_distance+Fault_length+dr_over_tip_distance ; % space where DR can occurr
SpaceDRLength = Fault_length; % space where DR can occurr

R2_lengths_PCTS = readtable(['TABLE_outputs/',SoF,'_R2_histogram_length.txt'],'VariableNamingRule', 'preserve');
R2_LOGN = readtable(['TABLE_outputs/',SoF,'LOGN_MU_SIGMA.txt'],'VariableNamingRule', 'preserve');

if strcmp(HWFW,'HW')==1
        DRlengths_min_max = [round(R2_lengths_PCTS{1,2}),round(R2_lengths_PCTS{1,3})]; % minimum and maximum length of DR to be simulated
    logn_MU = R2_LOGN{1,2};logn_SIGMA = R2_LOGN{1,3}; 
elseif strcmp(HWFW,'FW')==1
        DRlengths_min_max = [round(R2_lengths_PCTS{2,2}),round(R2_lengths_PCTS{2,3})]; % minimum and maximum length of DR to be simulated
    logn_MU = R2_LOGN{2,2};logn_SIGMA = R2_LOGN{2,3};     
end


Site_pos = [round(Fault_length/2)-along_strike_dimension/2 round(Fault_length/2)+along_strike_dimension/2] ;% starting and ending position of the site respct to the left-tip
zoom_x = [Site_pos(1)-along_strike_dimension*2 Site_pos(2)+along_strike_dimension*2]; % used to reduce dimensions in the figure
%zoom_x = [0 Fault_length];
% ratio between DRlength and PFlength
if site_distance > 200
fraction = load(fullfile('TABLE_outputs',[SoF,'RATIOpcts_farfault',char(HWFW),'sim.txt']))
elseif site_distance <= 200
fraction = load(fullfile('TABLE_outputs',[SoF,'RATIOpcts_nearfault',char(HWFW),'sim.txt']))
end

fcol = find(fraction(1,:)==along_strike_dimension);
F = fraction(4,fcol);%ratio between DRlength and PFlength
total_DR_lenght = Fault_length*F;% (in meters)


%% inizialize matrix and variables
Points_Sim_uniform(1:Simulations,1:SpaceDRLength)=NaN;
Points_Sim_exp(1:Simulations,1:SpaceDRLength)=NaN;
Total_cases = Simulations ;
Rupture_yes_case_unif =  0 ;
Rupture_yes_case_cluster =  0 ;

%% random sampling of the lenghts of DR
if SL ==1
mu = mean([DRlengths_min_max]);
t1 = min(DRlengths_min_max); 
t2 = max(DRlengths_min_max); 
truncPD = truncate(makedist('Normal','mu',mu,'sigma',(t2-t1)/2),t1,t2);
ypdf = truncPD.pdf(t1:t2);
elseif SL ==2
t1 = min(DRlengths_min_max); 
t2 = max(DRlengths_min_max); 
truncPD = truncate(makedist('Lognormal','mu',logn_MU,'sigma',logn_SIGMA),t1,t2);
ypdf = truncPD.pdf(t1:t2);
end
% cumulative lengths
for k = 1:Simulations
DRlengths(:,k) = [round(randsample(t1:t2, 1000, true, ypdf))]';
DRcumlengths(:,k) = cumsum(DRlengths(:,k),1);
DRsemilengths(:,k) = DRlengths(:,k)./2; % used later to keep tips separated
DRcumsemilengths(:,k) = cumsum(DRsemilengths(:,k),1);
f_numDR(:,k) = find(DRcumlengths(:,k) >= total_DR_lenght,1,'first'); %to fit total_DR_lenght
end




%% place rupture along the strike of the fault
for i = 1 : Simulations
    
 
    NumSeg = f_numDR(1,i);
    semiLSeg = DRlengths(1:NumSeg,i)/2;
    %% uniform distribution
    Centro_Seg_Simulati_unif = randperm(SpaceDRLength,NumSeg);
    Centro_Seg_Simulati_unif = sort(Centro_Seg_Simulati_unif);
    for g = 2 : length(Centro_Seg_Simulati_unif)
        if (Centro_Seg_Simulati_unif(g) - Centro_Seg_Simulati_unif(g-1)) < DRcumsemilengths(g-1,i)
            Centro_Seg_Simulati_unif(g) = Centro_Seg_Simulati_unif(g) + DRcumsemilengths(g-1,i);
        end
    end
    Punti_Simulati_unif = [];
    
    for j = 1:NumSeg
        Punti_Simulati_unif = [Punti_Simulati_unif,(Centro_Seg_Simulati_unif(j)-semiLSeg(j)):1:(Centro_Seg_Simulati_unif(j)+semiLSeg(j))];
    
    end
    Points_Sim_uniform(i,1:(length(Punti_Simulati_unif))) = [Punti_Simulati_unif];
    %% exponential distribution
    mean_distance = SpaceDRLength/NumSeg; % mean distance between centres
    starting_pos_of_DR = randi(SpaceDRLength,1); 
    interdistance_Seg_Simulati = exprnd(mean_distance,NumSeg,1);

    Punti_Simulati_cluster = [];
    ini_s = 0;
    ini_s = ini_s + starting_pos_of_DR;
    for j = 1:NumSeg
        
        end_s = ini_s + DRlengths(j,i);
        if ini_s > SpaceDRLength | end_s > SpaceDRLength
            ini_s = ini_s- SpaceDRLength;
            end_s = end_s - SpaceDRLength;
        end
        Punti_Simulati_cluster = [Punti_Simulati_cluster,ini_s:1:end_s];
        ini_s = end_s+ round(interdistance_Seg_Simulati(j));
    end
    Points_Sim_exp(i,1:(length(Punti_Simulati_cluster))) = [Punti_Simulati_cluster];
    
%% check if at aleast a point of DR is within the site
    if sum( Punti_Simulati_unif >= Site_pos(1) & Punti_Simulati_unif <= Site_pos(2) ) > 0
         Rupture_yes_case_unif = Rupture_yes_case_unif + 1 ;
    end
    
    if sum( Punti_Simulati_cluster >= Site_pos(1) & Punti_Simulati_cluster <= Site_pos(2) ) > 0
         Rupture_yes_case_cluster = Rupture_yes_case_cluster + 1 ;
    end
end
% calculate probabilities 
Prob_unif = Rupture_yes_case_unif / Total_cases;
Prob_cluster = Rupture_yes_case_cluster / Total_cases;
%%
Pname = fullfile('TABLE_outputs',[SoF,'_P_montecarlo_SiteDim',num2str(along_strike_dimension),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.txt']);
Pout = table(Prob_unif,Prob_cluster,(Prob_unif+Prob_cluster)/2,'VariableNames',{'Punif','Pexp','Pmean'});
writetable(Pout,Pname);
%%
figure(1)

for i = 1:Simulations_figure
    subplot(1,2,1)
    hold on
    plot(Points_Sim_uniform(i,:),repmat(i,1,size(Points_Sim_uniform,2)),'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','none')

    subplot(1,2,2)
    hold on
    plot(Points_Sim_exp(i,:),repmat(i,1,size(Points_Sim_exp,2)),'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','none')
end

subplot(1,2,1)
line([0 Fault_length],[0 0],'color','r','LineWidth',4)
line([Site_pos(1) Site_pos(1)] ,[0 i],'color','b','LineStyle',':','LineWidth',2)
line([Site_pos(2) Site_pos(2)] ,[0 i],'color','b','LineStyle',':','LineWidth',2)
hold off
xlim([zoom_x(1) zoom_x(2)])
set(gca,'XTick',[],'YTick',0:10:Simulations_figure)
xlabel('Distance along the strike of the PF (m)')
ylabel('simulations')

subplot(1,2,2)
line([0 Fault_length],[0 0],'color','r','LineWidth',4)
line([Site_pos(1) Site_pos(1)] ,[0 i],'color','b','LineStyle',':','LineWidth',2)
line([Site_pos(2) Site_pos(2)] ,[0 i],'color','b','LineStyle',':','LineWidth',2)
hold off
xlim([zoom_x(1) zoom_x(2)])
set(gca,'XTick',[],'YTick',0:10:Simulations_figure)
xlabel('Distance along the strike of the PF (m)')
ylabel('simulations')
%%
saveas(1,fullfile(path2,[SoF,'_simulations_momntecarlo_SiteDim',num2str(along_strike_dimension),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.pdf']),'pdf')
