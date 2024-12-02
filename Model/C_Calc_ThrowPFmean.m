%% this code uses empirical regressions to calculate ThrowPFmean
clear all
clc
%%
% input_ length of PF and magnitude, position and distance of the site
Length_PF = 40000; %m 
magnitude_PF = 7;
position = 0.5; % this is l/L
distance = 2; %km
dip = 60; %degree
kin = 'Normal'; 
%kin = 'Reverse';

pathout1 = 'TABLE_outputs';
pathout2 = 'Figure/displacement';
%% empirical scale regressions
regression = {'W&C94all','T&A17dip','Leo10dip'};

Md_wc94all = 10.^( -5.46 + 0.82 .* magnitude_PF); 
sigma_Md_wc94all = 0.42;

Ad_wc94all = 10.^( -4.8 + 0.69 .* magnitude_PF); 
sigma_Md_wc94all = 0.36;

%  Thingbaijam et al 2017
if strcmp(kin,'Reverse')==1
Ad_ta17 = 10.^(0.451*magnitude_PF - 3.156); % this is the AVERAGE
sigma_ad_ta17 = 0.149;

elseif strcmp(kin,'Normal')==1
Ad_ta17 = 10.^(0.693*magnitude_PF - 4.967); % this is the AVERAGE
sigma_ad_ta17 = 0.195;
end

Md_ta17 = 2*Ad_ta17;

%  Leonard, 2010 use square meteres for area
logMo = 9.1 +1.5.*magnitude_PF;
area_dip = (logMo - 6.1)/1.5;
area_dip_low = (logMo - 6.60)/1.5;
area_dip_upp = (logMo - 5.69)/1.5;
Ad_le10dip= 10.^(0.5 * (area_dip) -4.42 ); %this is the AVERAGE
Ad_le10dip_low=0.5 * (area_dip_upp) -4.82;
Ad_le10dip_upp=0.5 * (area_dip_low) -3.92;
diff = Ad_le10dip_upp- Ad_le10dip_low;
sigma_ad_le10dip = 0.5*mean(diff);
Md_le10dip = Ad_le10dip*2;
%% values of Max and AVG Displacement in meters
MaxD =[Md_wc94all;Md_ta17;Md_le10dip];
AD =[Ad_wc94all;Ad_ta17;Ad_le10dip];

MaxD = MaxD*100;
AD = AD*100;
%% values of the displacement from triang and tapered
x = [0:0.00001:0.5];
for i = 1:length(MaxD)
% triangular profile case
displacement_triang(i,:) = (x .* MaxD(i))/0.5;
% tapered-slip profile case
displacement_tapered(i,:) = 1.311.*AD(i).*sqrt((sin(pi.*x)));
end
% create a complete profile
X = [x,1-(fliplr(x(1:end-1)))];
D_triang = [displacement_triang,fliplr(displacement_triang(:,1:end-1))];
D_tapered = [displacement_tapered,fliplr(displacement_tapered(:,1:end-1))];

%% displacement at point
d_site_pos = find(X==position);
D_site = [D_triang(:,d_site_pos);D_tapered(:,d_site_pos)];
D_mean = mean([D_triang;D_tapered]);

%% radius of search for ThrowPFmean
l = Length_PF/2;
p = distance/2;

r = (p/l)*0.5;
left_point = round(position-r,5);
left_point(left_point<0)=0;

right_point = round(position+r,5);
right_point(right_point>1)=1;
%%
f1 = find(X==left_point);
if isempty(f1)
  f1=find(X<=left_point,1,'last');
end
f2 = find(X==right_point);
if isempty(f2)
 f2=find(X>=right_point,1,'first');
end

%%  calculate Throw and ThrowPFmean
T_triang = D_triang*sind(dip);
T_tapered = D_tapered*sind(dip);
T_mean = D_mean*sind(dip);
T_values = T_mean(f1:f2);
ThrowPFmean = mean(T_values)
save(fullfile(pathout1,[char(kin),'_',num2str(magnitude_PF),'_',num2str(position),'_',num2str(distance),'.txt']),'ThrowPFmean','-ascii')

%% figure
figure(1)
hold on
for i = 1:length(MaxD)

plot(X,T_triang(i,:),'LineWidth',1,'display',char(regression(i)));
plot(X,T_tapered(i,:),'LineWidth',1,'display',char(regression(i)));

end
plot(X,T_mean,'color','k','LineWidth',2,'display','mean profile');

line([position position],[0, max(D_site)], 'display','site position')
line([left_point left_point],[0, max(D_site)], 'display','left limit')
line([right_point right_point],[0, max(D_site)], 'display','right limit')
line([left_point right_point],[ThrowPFmean, ThrowPFmean],'color','k', 'display',['ThrowPFmean:',num2str(round(ThrowPFmean,2))])

xlabel('x/L')
ylabel ('Throw [cm]')

legend('show')
grid on
set (gca,'fontsize',12)
saveas(1,fullfile(pathout2,[char(kin),'_',num2str(magnitude_PF),'_',num2str(position),'_',num2str(distance),'.png']),'png')
