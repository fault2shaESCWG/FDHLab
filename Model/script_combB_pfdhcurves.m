% this code plot the decision-tree case study

clear all
clc
%% 
% output paths
pathout1 = 'DecisionTree/';
if isempty(dir(pathout1))
mkdir(pathout1)
end
%%
VDL=0.01:0.01:3; % vertical displacement levels meter
eps=3;

%% Scenario
Mw_scenario     = [7];
VDmain_scenario = [1.42];
SiteDim = 100;
dist = [500];

hwfw=2; %1 = footwall, 2 = hangingwall
SoF = 2; % 2 = normal
T_logistic = readtable('../Regressions/TABLE_outputs/LOGISTIC_CombB_Mw7_SiteDim100_SiteDist500_HW.txt');

T_regression = readtable('../Regressions/TABLE_outputs/coefficients_throw.txt');
StdDev = load('../Regressions/TABLE_outputs/sigma.txt');
% attenuation of Vertical Displacement
% Name	Estimate
Interceptcoeff = T_regression.value(1);
Mwcoeff= T_regression.value(2);
HWFW_2coeff = T_regression.value(3);
SoF_2coeff = T_regression.value(4);
combination_2coeff = T_regression.value(5);
combination_3coeff = T_regression.value(6);
distanceLNcoeff = T_regression.value(7);
ThrowPFmeanLNcoeff = T_regression.value(8);

%% Prob logistic
Plogistic = T_logistic{find(T_logistic{:,1} <= dist,1,'last'),2};

%% Prob Vertical Displacement
dummy = ones(1,8);
dummy(8) = 0; % comb B
coeff = [Interceptcoeff(1),Mwcoeff(1),distanceLNcoeff(1),ThrowPFmeanLNcoeff(1),HWFW_2coeff(1),SoF_2coeff(1),combination_2coeff(1),combination_3coeff(1)];
    if hwfw == 2 
        dummy(5) = 0;
    end
    if SoF == 1 % Reverse
        dummy(6) = 0;
    end

%%
HAZARD_CURVES=[]; COND_PROB_CURVES=[];eccedenza=[];COND_PROB=[];
x = 0.01:0.01:(VDL(end)*2);



lnY = coeff(1) + coeff(2)*Mw_scenario + coeff(3)*log(dist) + coeff(4)*log(VDmain_scenario) + sum(coeff(5:8).*dummy(5:8));
exp(lnY)
pdY=[];truncY=[];
pdY = makedist('normal','mu',(lnY),'sigma', StdDev);
truncY =truncate(pdY,-eps*(StdDev)+lnY,eps*(StdDev)+lnY);
PY=cdf(truncY,log(x));

    for j=1: length(VDL)
    eccedenza(1,j)=1-PY(find(log(x)<=log(VDL(j)),1,'last'));  
    end


HAZARD_CURVES(1,:) = [dist,Mw_scenario,VDmain_scenario, eccedenza];
COND_PROB_CURVES(1,:) =[dist,Mw_scenario,VDmain_scenario,Plogistic(1,1),Plogistic(1,1).*eccedenza];


%%       


COND_PROB_CURVES(COND_PROB_CURVES<0)=0;
COND_PROB = 1 - prod(1-COND_PROB_CURVES(:,5:end),1);

figure(1)
hold on

plot(VDL,COND_PROB,'-','LineWidth',2,'color',rand(1,3),'display',...
    strcat('comb:B',',Mw:',num2str(Mw_scenario),...
    ',VD:',num2str(VDmain_scenario)))

grid minor
xlabel('vd (m)')
ylabel('cond probability of exceedance')
set(gca,'XScale','log','YScale','log','FontSize',12)
legend('show','location','NorthEast')
out(:,1) = [COND_PROB'];



writematrix([VDL',out],'Valori_probability_combB.txt')

 saveas(1,'COMB_B.pdf','pdf')
