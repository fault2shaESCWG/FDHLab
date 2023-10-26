%% questo codice fa il merge delle tabelle ottenute dallo script_A
clear all
clc
%%
pathtable = 'TABLE_db_20231026';
%%
event_rev  = load('list_Reverse.txt');
event_nor  = load('list_Normal.txt');
id_reverse = [event_rev(:,1)];
id_normal = [event_nor(:,1)];
%% r tables
header = {'IdE','IdR','latR','lonR','latPr','lonPr','rankR','lengthR','rankPr'};
%%
if exist('id_reverse','var') 
T_reverse = []
for i = 1:length(id_reverse)
    id_list = num2str(id_reverse(i,1))
   
    T=[];
    T = readtable(fullfile(pathtable,strcat(char(id_list),'_r_distance_table_20231018.txt')));
    T_reverse = [T_reverse;T];
end
writetable(T_reverse,fullfile(pathtable,'Reverse_r_distance_table.txt'))

end

%%
if exist('id_normal','var')
T_normal = []
for i = 1:length(id_normal)
    id_list = num2str(id_normal(i))
   
    T=[];
    T = readtable(fullfile(pathtable,strcat(char(id_list),'_r_distance_table_20231018.txt')));
    T_normal = [T_normal;T];
end

writetable(T_normal,fullfile(pathtable,'Normal_r_distance_table.txt'))
end
%%
