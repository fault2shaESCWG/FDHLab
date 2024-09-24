% calculate logistic 
clear all
pathout1 = 'FIGURE/Logistic';
mkdir(pathout1)
%%   user defined inputs
nameforfigure = 'CombB_Normal_m7_dist500_dim100';
% calculate logistic for a given combination in "param_logistic"
m = 7;
x = 5:10:8000;
site_dim = 100;
site_distance = 500;% meters from the PF
HWFW = 'HW'; % HW = Hanging wall; FW = Footwall location of the site
%SoF = 'Reverse';
SoF = 'Normal';
param_logistic = load(fullfile('../Regressions/TABLE_outputs',['parameters_logistic_multisizeC2_',char(SoF),'.txt']));

%%
Pname = fullfile('TABLE_outputs',[SoF,'_P_montecarlo_SiteDim',num2str(site_dim),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.txt']);
P_montecarlo_table = readtable(Pname)
P_montecarlo = P_montecarlo_table.Punif;
%%
i = find(param_logistic(1,:) ==site_dim);
if strcmp(HWFW,'HW')==1
y = param_logistic(2,i)+ param_logistic(3,i)*m + param_logistic(4,i).*x + param_logistic(5,i)*0;
elseif strcmp(HWFW,'FW')==1
y = param_logistic(2,i)+ param_logistic(3,i)*m + param_logistic(4,i).*x + param_logistic(5,i)*1;
end  
P =1-  exp(y)./(1+exp(y));
%%
condP = P.* P_montecarlo;
%%
p1 = find(x >= site_distance, 1,'first');
p2 = find(x >= (site_distance+site_dim), 1,'first');

figure(1)
hold on
plot(x,condP,'-k')
hleg1(1) =line([site_distance site_distance],[0 condP(p1)],'color','r','LineStyle','-','LineWidth',2,'display','nearest distance site PF (P used for PFDHA)');
hleg1(2) =line([site_distance+site_dim site_distance+site_dim],[0 condP(p2)],'color','b','LineStyle','--','LineWidth',2,'display','most distant part of the site from PF');
grid on
xlabel('distance (m)')
ylabel('Probability')
legend(hleg1)
saveas(1,fullfile(pathout1,[nameforfigure,'.png']),'png');
%%
output = table(x',condP','VariableNames',{'distance','Logistic'});
writetable(output,fullfile('TABLE_outputs',['LOGISTIC_CombB_Mw',num2str(m),'_SiteDim',num2str(site_dim),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.txt']))
%%


