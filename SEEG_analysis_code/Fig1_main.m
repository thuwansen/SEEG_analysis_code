% Code for Fig1
clc,clear,close all;
addpath('D:\software\toolbox\fieldtrip-20191213');
addpath(genpath('D:\software\toolbox\nature_walk_wave'));
addpath('D:\software\toolbox\BrainNetViewer_20191031')
code_path = 'D:\SEEG_preoject\BIDS_data\code\research_functions';

addpath(genpath([code_path,'\research_functions']));
run([code_path,'\Global_variable_define_LOC2ROC.m']);
fs =500;
% Define the colors for CLA, OB and FI.
color_group = [243/255,164/255,0/255;
               209/255,21/255,71/255;
               39/255,65/255,146/255];
      
savepath = '.\Fig1';
if ~exist(savepath)
    mkdir(savepath);
end

%% Fig 1b
subject_info = readtable(subject_info_path);
age = subject_info.Var2;
age_mean = nanmean(age);
age_std = nanstd(age);
disp(['The age of all of the subjests is:', num2str(age_mean),'±',num2str(age_std)]);

%% Fig 1c
Fig1c_path = [savepath,'\Fig1c'];
if ~exist(Fig1c_path)
    mkdir(Fig1c_path);
end
all_subjects_coordinate = [];
for idx =1:length(subs_name)
    if ~isempty(SEEG_coordinate_path{idx,1})
        tmp_coordinate = load(SEEG_coordinate_path{idx,1});   
        tmp_mni = tmp_coordinate.SEEG_COORDINATE.MNI_coordinate;
        [chan_num,~]=size(tmp_mni);
        final_coordinate = zeros(chan_num,6);
        final_coordinate(:,1:3)=tmp_mni;
        final_coordinate(:,4)=idx;
        final_coordinate(:,5)=idx;
        final_coordinate(:,6)=idx;
        all_subjects_coordinate=[all_subjects_coordinate; final_coordinate];
    end
end
save([Fig1c_path,'/all_subjects_coordinate.node'], 'all_subjects_coordinate','-ascii');
% Plot using Brainnet Viewer 

%% Fig 1d
Fig1d_path = [savepath,'\Fig1d'];
if ~exist(Fig1d_path)
    mkdir(Fig1d_path);
end
subjects = dir([dataset_path,'\sub*']);
leads_in_AAL=[];
for i=1:length(subjects)
    if exist([dataset_path,'\',subjects(i).name,'\ieeg\SEEG_COORDINATE.mat'])
        tmp = load([dataset_path,'\',subjects(i).name,'\ieeg\SEEG_COORDINATE.mat']);
        leads_in_AAL = [leads_in_AAL;tmp.SEEG_COORDINATE.AAL_index];
    end
end
leads_in_AAL(leads_in_AAL==0)=[];
AAL_statistics = tabulate(leads_in_AAL(:));
% Calculate the brain area (AAL) coverage of electrode contacts
Cerebrum_AAL_coverage=0;
for i=1:90
    if AAL_statistics(i,2)~=0
        Cerebrum_AAL_coverage=Cerebrum_AAL_coverage+1;
    end
end
Cerebrum_AAL_coverage_ratio = Cerebrum_AAL_coverage/90; 
% Top 15 electrode contact numbers for each brain area.
top_15_numbers = zeros(15,4);
AAL_statistics_sorted = sortrows(AAL_statistics,2,'descend');
flag = 0;
for i=1:length(AAL_statistics_sorted)
    if flag==15
        break;
    end
    AAL_idx = AAL_statistics_sorted(i);
    AAL_idx_left_brain =2*round(AAL_idx/2)-1;
    AAL_idx_right_brain =2*round(AAL_idx/2);
    
    if ~ismember(AAL_idx_left_brain,top_15_numbers)&&~ismember(AAL_idx_right_brain,top_15_numbers)
        flag=flag+1;
        top_15_numbers(flag,1)=AAL_idx_left_brain;
        top_15_numbers(flag,3)=AAL_idx_right_brain;
        tmp_idx = find(AAL_statistics_sorted(:,1)==AAL_idx_left_brain);
        top_15_numbers(flag,2)=AAL_statistics_sorted(tmp_idx,2);
        tmp_idx = find(AAL_statistics_sorted(:,1)==AAL_idx_right_brain);
        top_15_numbers(flag,4)=AAL_statistics_sorted(tmp_idx,2);       
    end
end

top_15_subject_number=zeros(15,4);
top_15_subject_number(:,[1,3])=top_15_numbers(:,[1,3]);
for i=1:length(subjects)
    if exist([dataset_path,'\',subjects(i).name,'\ieeg\SEEG_COORDINATE.mat'])
        tmp = load([dataset_path,'\',subjects(i).name,'\ieeg\SEEG_COORDINATE.mat']);
        tmp_leads = tmp.SEEG_COORDINATE.AAL_index;
        for j=1:15
            if ismember(top_15_subject_number(j,1),tmp_leads)
                top_15_subject_number(j,2)=top_15_subject_number(j,2)+1;
            end
            if ismember(top_15_subject_number(j,3),tmp_leads)
                top_15_subject_number(j,4)=top_15_subject_number(j,4)+1;
            end
        end
    end
end
AAL_IDX =  top_15_subject_number(:,1);
left_brain = top_15_numbers(:,2);
left_brain_subject = top_15_subject_number(:,2);
right_brain = top_15_numbers(:,4);
right_brain_subject = top_15_subject_number(:,4);
total_brain = left_brain+right_brain;
total_brain_subject=left_brain_subject+right_brain_subject;
save([Fig1d_path,'\Fig1d_AAL_statistic.mat'],'AAL_IDX','left_brain','left_brain_subject',...
    'right_brain','right_brain_subject','total_brain','total_brain_subject');

figure(), set(gcf,'position',[100,100,600,800]);
bar(total_brain);
hold on;bar(left_brain);
xticklabels(AAL_IDX);
ylabel('Numbers of electrode contacts');
xlabel('AAL index');
view(-90,-90);
set(gca,'FontName','Arial','FontSize',15,'LineWidth',1);
saveas(gcf,[Fig1d_path,'/Fig1d_lead_contact.tif']);

figure(), set(gcf,'position',[100,100,600,800]);
bar(total_brain_subject);
hold on;bar(left_brain_subject);
xticklabels(AAL_IDX);
ylabel('Numbers of subjects');
xlabel('AAL index');
view(90,90);
set(gca,'FontName','Arial','FontSize',15,'LineWidth',1);
saveas(gcf,[Fig1d_path,'/Fig1d_subject.tif']);

%% Fig 1e
Fig1e_path =[savepath '\Fig1e'];
if ~exist(Fig1e_path)
    mkdir(Fig1e_path);
end
subject_info = readtable(subject_info_path);
t_start_stopmedicine = subject_info.Var5- subject_info.Var4;
t_stopmedicine_ROC = subject_info.Var6- subject_info.Var5;
t_ROC_CON = subject_info.Var7- subject_info.Var6;
time_list = [t_start_stopmedicine;t_stopmedicine_ROC;t_ROC_CON];
hue = cell(length(time_list),1);
hue(:,1) = cellstr(['GA']);
x= cell(length(time_list),1);
x(1:length(t_start_stopmedicine),1)=cellstr(['Phase I']);
x(length(t_start_stopmedicine)+1:length(t_start_stopmedicine)+length(t_stopmedicine_ROC),1)=cellstr(['Phase II']);
x(length(t_start_stopmedicine)+length(t_stopmedicine_ROC)+1:length(time_list),1)=cellstr(['Phase III']);
figure(), 
boxplot([t_start_stopmedicine,t_stopmedicine_ROC,t_ROC_CON]);
xticklabels({'Phase I';'Phase II';'Phase III'});
ylabel('Time (s)');
ylim([0,2000]);
set(gca,'FontName','Arial','FontSize',15,'LineWidth',1);
saveas(gcf,[Fig1e_path,'/Fig1e_time_statistic_matlab.tif']);
save([Fig1e_path,'/Fig1e_4_python_plot.mat'],'time_list','hue','x');

nanmean(t_start_stopmedicine(1:39,1)),nanstd(t_start_stopmedicine(1:39,1))
nanmean(t_stopmedicine_ROC(1:39,1)),nanstd(t_stopmedicine_ROC(1:39,1))
nanmean(t_ROC_CON(1:39,1)),nanstd(t_ROC_CON(1:39,1))
