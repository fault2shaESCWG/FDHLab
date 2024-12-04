%% plot 3 combinations from script_Comb A, B, C 

clear all
clc
%%
T_CombA = readtable('Valori_probability_combA.txt');
T_CombB = readtable('Valori_probability_combB.txt');
T_CombC = readtable('Valori_probability_combC.txt');
%%
case_1 = 1-(1-T_CombA{:,2}).*(1-T_CombB{:,2}).*(1-T_CombC{:,2});
case_2 = 1-(1-T_CombA{:,2}).*(1-T_CombB{:,2});
case_3 = T_CombA{:,2};
CombA =  T_CombA{:,2};
CombB =  T_CombB{:,2};
CombC =  T_CombC{:,2};
x =  T_CombA{:,1};

figure(1)
hold on
plot(x,case_1,'-','LineWidth',2,'color','r','display','case 1')
plot(x,case_2,'-','LineWidth',2,'color','k','display','case 2')
plot(x,case_3,'-','LineWidth',2,'color','m','display','case 3')

plot(x,CombA,'-','LineWidth',1,'color','b','display','CombA')
plot(x,CombB,'-','LineWidth',1,'color','c','display','CombB')
plot(x,CombC,'-','LineWidth',1,'color','g','display','CombC')


grid minor
xlabel('throw (m)')
ylabel('cond probability of exceedance')
set(gca,'XScale','log','YScale','log','FontSize',12)
legend('show','location','NorthEast')
xlim([0.01 10])
ylim([1e-5 1])
axis square


%saveas(1,'DecisionTree_Paper.pdf','pdf')
