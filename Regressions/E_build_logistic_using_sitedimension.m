% calculate logistic 
clear all

%%   user defined inputs
m = 6.5;
x = 5:10:5000;
site_dim = 500;
site_distance = 100;% meters from the PF
HWFW = 'HW'; % HW = Hanging wall; FW = Footwall location of the site
kin = 'Reverse';
%kin = 'Normal';
%%
Pname = fullfile('TABLE_outputs',['P_montecarlo_SiteDim',num2str(site_dim),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.txt']);
P_montecarlo_table = readtable(Pname)
P_montecarlo = P_montecarlo_table.Punif;
%%
param_logistic = load(fullfile('TABLE_outputs',['parameters_logistic_multisizeC1_',char(kin),'.txt']));
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
output = table(x',condP','VariableNames',{'distance','Logistic'});
writetable(output,fullfile('TABLE_outputs',['LOGISTIC_Mw',num2str(m),'_SiteDim',num2str(site_dim),'_SiteDist',num2str(site_distance),'_',char(HWFW),'.txt']))
%%


