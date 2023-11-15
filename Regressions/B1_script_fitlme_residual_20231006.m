%% this code computes  multiple regression to estimate vd on DR

clc
clear all
close all

%%

pathout1 ='../Analyse_SUREvers2_database/TABLE_db_20231026';
path2 = fullfile('FIGURE','residual');

if isempty(dir(path2))
mkdir(path2)
end

pathinputs = '../Analyse_SUREvers2_database/TABLE_db_20231026'
T1=readtable(fullfile(pathinputs,'Reverse_s_distance.txt'));
T2=readtable(fullfile(pathinputs,'Normal_s_distance.txt'));
%%
T = [T1;T2];

rmv = find(T.ThrowPFmean==0);
T(rmv,:) = [];

%% combinations

% rank DR-PF 
C1 = find(T.RankDR == 2 & T.RankPF == 1);
C2 = find(T.RankDR == 2 & T.RankPF == 1.5);
C3 = find(T.RankDR == 1.5); % always vs rankPF ==1
C4 = find(T.RankDR > 2);% always v rankPF ==1
%%
% we use 3 combinations
T.combination(C1) = 1; 
T.combination(C2) = 2; 
T.combination(C3) = 3;
T.combination(C4) = 3; 

%%
T.combination = categorical(T.combination);
T.HWFW(T.HWFW  == -1) = 2; % footwall

T.HWFW = categorical(T.HWFW);
T = renamevars(T,["kinR_1kinN_2"], ["SoF"]);
T.SoF = categorical(T.SoF);
%% convert in log
T.distanceLN = log(T.distance);
T.ThrowPFmeanLN = log(T.ThrowPFmean);
T.ThrowDRLN = log(T.ThrowDR);

%%
lme1 = fitlme(T,'ThrowDRLN ~ Mw');
lme2 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN');
lme3 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN');
lme4 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN + HWFW');
lme5 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN + HWFW + combination');
lme6 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN + HWFW + SoF + combination');
lme7 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN + HWFW + SoF + combination + (1|IdE)');
lme8 = fitlme(T,'ThrowDRLN ~ Mw + distanceLN + ThrowPFmeanLN + combination + (1|IdE)');
lme9 = fitlme(T,'ThrowDRLN ~ Mw +  ThrowPFmeanLN  + (distanceLN:SoF:HWFW)+ combination + (1|IdE)');

%%
AIC_model  =[];
AIC_model  = {'lme1',lme1.ModelCriterion.AIC;'lme2',lme2.ModelCriterion.AIC;
    'lme3',lme3.ModelCriterion.AIC;'lme4',lme4.ModelCriterion.AIC;
    'lme5',lme5.ModelCriterion.AIC;'lme6',lme6.ModelCriterion.AIC;
    'lme7',lme7.ModelCriterion.AIC;'lme8',lme8.ModelCriterion.AIC;
    'lme9',lme9.ModelCriterion.AIC};
%%
coeff = table(lme7.Coefficients.Name,lme7.Coefficients.Estimate,lme7.Coefficients.pValue,'VariableNames',{'name','value','p-value'});
sigma_value = sqrt(lme7.MSE);
writetable(coeff,fullfile(pathout1,'coefficients_throw.txt'));
save(fullfile(pathout1,'sigma.txt'),'sigma_value','-ascii');
%% residuals

T.lme7residuals = lme7.residuals;
T.lme7predict=lme7.predict;

%% convert in log
T.distanceLN = log(T.distance);
T.ThrowPFmeanLN = log(T.ThrowPFmean);
T.ThrowDRLN = log(T.ThrowDR);

%%  fit curve
polydata = polyfit(T.distance, T.lme7residuals,1);
fitcurve = polydata(1).*T.distance + polydata(2);
polydata2 = polyfit(T.distance, (T.ThrowDR-exp(T.lme7predict)),1);
fitcurve2 = polydata2(1).*T.distance + polydata2(2);

hleg = [];

figure(3)

subplot(3,2,1)
hold on
histogram(T.lme7residuals,'FaceColor',[.5 .5 .5])
xlim([-4 4])
xlabel('residuals')
ylabel('count')
hold off

subplot(3,2,2)
hold on
plot(T.ThrowDRLN,T.lme7predict,'.','MarkerEdgeColor',[.5 .5 .5])
grid on
ylim([-7 2])
xlim([-7 2])
axis square
ylabel('ln predicted')
xlabel('ln observed')
hold off

subplot(3,2,3:4)
hold on
hleg(1) = plot(T.distance, T.lme7residuals,'ok','MarkerFaceColor',[.5 .5 .5],'Display','all data');
hleg(2) = plot(T.distance,fitcurve,'-k','Display','fit');
grid on
ylim([-4 4])
legend([hleg],'location','best')
ylabel('ln residual (data – fit)')
xlabel ('distance (m)')
hold off

subplot(3,2,5:6)
hold on
plot(T.distance, (T.ThrowDR-exp(T.lme7predict)) ,'ok','MarkerFaceColor',[.5 .5 .5],'Display','all data');
plot(T.distance,fitcurve2,'-k','Display','fit');
grid on
ylim([-5 5])

ylabel('residual (data – fit)')
xlabel ('distance (m)')
hold off
saveas(3,fullfile(path2,'Residuals_vs_distance.pdf'),'pdf')

%%
col = colormap('jet');

for s = 1:2

       
        Ts = T(T.SoF == num2str(s),:);
        
        
        ide = unique(Ts.IdE);
        ide = flipud(ide);
        pos_col = round(linspace(1,size(col,1),length(ide)));
        col_sof = col(pos_col,:);
        
        figure(s)
        hold on
        if s ==1
        title('Reverse')
        elseif s ==2
        title('Normal') 
        end
        polydata = polyfit(Ts.distance, Ts.lme7residuals,1);
        fitcurve = polydata(1).*Ts.distance + polydata(2);
        subplot(3,2,1:2)
        hold on
        for  i =1:length(ide)
            Ttemp = [];
            Ttemp = Ts(Ts.IdE == ide(i),:);
        
        plot(Ttemp.distance, Ttemp.lme7residuals,'ok','MarkerFaceColor',col_sof(i,:),'Display',num2str(ide(i)))
        end
        plot(Ts.distance,fitcurve,'-k')
        grid on
        ylim([-4 4])
        ylabel('ln residual (data – fit)')
        xlabel ('distance (m)')
        legend('show','location','best','fontsize',4)
        hold off
        
        
        for combi = 1:3
        Tcombi = [];
        Tcombi = Ts(Ts.combination == num2str(combi),:);
        polydata = polyfit(Tcombi.distance, Tcombi.lme7residuals,1);
        fitcurve = polydata(1).*Tcombi.distance + polydata(2);
        
            for  i =1:length(ide)
            Tide = [];
            Tide = Tcombi(Tcombi.IdE == ide(i),:);
        
            subplot(3,2,combi+2)
            hold on
            title(strcat('combination:',num2str(combi)))
            plot(Tide.distance, Tide.lme7residuals,'ok','MarkerFaceColor',col_sof(i,:),'Display',num2str(ide(i)))
            end
            
        plot(Tcombi.distance,fitcurve,'-k')
        grid on
        ylim([-4 4])
       % legend('show','location','best')
        ylabel('ln residual (data – fit)')
        xlabel ('distance (m)')
   end

end
saveas(1,fullfile(path2,'REVERSE_combination_vs_distanza.pdf'),'pdf')
saveas(2,fullfile(path2,'NORMAL_combination_vs_distanza.pdf'),'pdf')