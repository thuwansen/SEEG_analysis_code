% Code for Fig4 and Supplementary Fig3
clc,clear,close all;
addpath('D:\software\toolbox\fieldtrip-20191213');
addpath(genpath('D:\software\toolbox\nature_walk_wave'));
addpath('D:\software\toolbox\BrainNetViewer_20191031');
addpath(genpath('D:\software\PTE-variants-master'));
code_path = 'D:\SEEG_preoject\BIDS_data\code\research_functions';
addpath(genpath([code_path,'\research_functions']));
addpath(genpath('D:\software\toolbox\PAC\PACmeg-master'));
run([code_path,'\Global_variable_define_LOC2ROC.m']);
fs =500;
% Define the colors for CLA, OB and FI.
color_group = [243/255,164/255,0/255;
               209/255,21/255,71/255;
               39/255,65/255,146/255];
      
savepath = '.\Fig4';
if ~exist(savepath)
    mkdir(savepath);
end  

%% Fig. 4a top
Fig4a_path =[savepath '\Fig4a'];
if ~exist(Fig4a_path)
    mkdir(Fig4a_path);
end
idx = 2;
t_win = 30;
t_s_ana = 250;
load(data_MO_name_ana{idx});
data_ana = SEEG;   
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq = [0.3,1]; 
cfg.channel     = 'all';
cfg.bpfilttype = 'firws';
data_SO = ft_preprocessing(cfg, data_ana); 
% select and plot cla ob  fi channel
channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
channel_se_idx = zeros(3,1);
for i = 1:length(channel_se)
    channel_se_idx(i,1) = find(ismember(SEEG.label,channel_se{i}));
end
data_ana_se = SEEG.trial{1}(channel_se_idx,fs*t_s_ana:fs*(t_s_ana+t_win));
data_ana_se_SO = data_SO.trial{1}(channel_se_idx,fs*t_s_ana:fs*(t_s_ana+t_win));
figure(),set(gcf,'position',[100,100,800,400]);
set(gca,'position',[0 0 1 1])
plot(data_ana_se(1,:)-1000,'color',[243/255,164/255,0/255]);xlim([0,length(data_ana_se)]);
hold on; plot(data_ana_se_SO(1,:)-1000,'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
hold on;plot(data_ana_se(2,:)-500,'color',[209/255,21/255,71/255]);xlim([0,length(data_ana_se)]);
hold on;plot(data_ana_se_SO(2,:)-500,'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
hold on;plot(data_ana_se(3,:),'color',[39/255,65/255,146/255]);xlim([0,length(data_ana_se)]);
hold on;plot(data_ana_se_SO(3,:),'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
axis off;
saveas(gcf,[Fig4a_path,'/Fig4a_top.tif']);

%% Fig.4a bottom
SO_mask = zeros(size(data_ana_se_SO));
t=zeros(size(data_ana_se_SO));
for i=1:3
    [s_SO_mask_tmp,SO_de] = SO_detect(data_ana_se_SO(i,:),[0.4,5],0.5,fs);
    SO_mask(i,:)=s_SO_mask_tmp;
    t_temp = 1:length(s_SO_mask_tmp);
    t_temp(s_SO_mask_tmp==0)=nan;
    t(i,:) = t_temp;
end
% calculate the overlap of SO
s_SO_mask_overlap = SO_mask(1,:).*SO_mask(2,:).*SO_mask(3,:);
s_SO_mask_union = double(SO_mask(1,:)|SO_mask(2,:)|SO_mask(3,:));
s_SO_mask_overlap_union = zeros(size(s_SO_mask_union));
st = 0;
ed = 0;
for i=1:length(s_SO_mask_union)-1
    if(s_SO_mask_union(1,i)==0)&&(s_SO_mask_union(1,i+1)==1)
        st = i;
    end
    if(s_SO_mask_union(1,i)==1)&&(s_SO_mask_union(1,i+1)==0)
        ed = i;
    end
    if st<ed
        if sum(s_SO_mask_overlap(1,st:ed))>0
            s_SO_mask_overlap_union(1,st:ed)=1;
        end
    end
end
s_SO_mask_selected = SO_mask;
t_overlap=zeros(size(SO_mask));
for i=1:3
    tmp_mask = SO_mask(i,:).*s_SO_mask_overlap_union;
    t_temp = 1:length(s_SO_mask_overlap);
    t_temp(tmp_mask==0)=nan;
    t_overlap(i,:)=t_temp;
end
color2=[255/255,0/255,128/255];
color1=[22/255,214/255,250/255];
figure(),set(gcf,'position',[100,100,800,400]);
set(gca,'position',[0 0 1 1])
plot(data_ana_se_SO(3,:),'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
hold on; plot(zeros(size(data_ana_se_SO(3,:))),'--','color',[128/255,128/255,128/255,],'LineWidth',1.5);
hold on;plot(t(3,:),data_ana_se_SO(3,:),'color',color1,'LineWidth',1.5);
hold on; plot(t_overlap(3,:),data_ana_se_SO(3,:),'color',color2,'LineWidth',1.5);
shift1=500;
plot(data_ana_se_SO(2,:)-shift1,'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
hold on; plot(zeros(size(data_ana_se_SO(2,:)))-shift1,'--','color',[128/255,128/255,128/255],'LineWidth',1.5);
hold on;plot(t(2,:),data_ana_se_SO(2,:)-shift1,'color',color1,'LineWidth',1.5);
hold on; plot(t_overlap(2,:),data_ana_se_SO(2,:)-shift1,'color',color2,'LineWidth',1.5);
shift2=1000;
plot(data_ana_se_SO(1,:)-shift2,'color','k','LineWidth',1.5);xlim([0,length(data_ana_se)]);
hold on; plot(zeros(size(data_ana_se_SO(1,:)))-shift2,'--','color',[128/255,128/255,128/255],'LineWidth',1.5);
hold on;plot(t(1,:),data_ana_se_SO(1,:)-shift2,'color',color1,'LineWidth',1.5);
hold on; plot(t_overlap(1,:),data_ana_se_SO(1,:)-shift2,'color',color2,'LineWidth',1.5);
axis off;
xlim([7000,13000]);
saveas(gcf,[Fig4a_path,'/Fig4a_bottom.tif']);

%% Prepare data for group analysis 
data_SO_4_analysis = cell(length(data_MO_name_ana),3);
% 读取并存储数据
for idx = 1:length(data_MO_name_ana)
    t_win = 100;
    t_s_ana = ana_time{idx,1};
    load(data_MO_name_ana{idx});
    data_ana = SEEG;   
    cfg = [];
    cfg.bpfilter  = 'yes';
    cfg.bpfreq = [0.3,1]; 
    cfg.channel     = 'all';
    cfg.bpfilttype = 'firws';
    data_SO = ft_preprocessing(cfg, data_ana); 
    channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
    for i = 1:length(channel_se)
        if ~isempty(channel_se{i})
            tmp_channel = find(ismember(SEEG.label,channel_se{i}));
            t_s = data_SO.trial{1}(tmp_channel,t_s_ana*fs:(t_s_ana+t_win)*fs);
            data_SO_4_analysis(idx,i) =mat2cell(t_s,1);
        end
    end 
end
save([savepath,'\data_SO_4_analysis.mat'],'data_SO_4_analysis');

%% Fig. 4b
Fig4b_path = [savepath,'\Fig4b'];
if ~exist(Fig4b_path)
    mkdir(Fig4b_path);
end
corr_coefficient = nan(length(data_MO_name_ana),2);
corr_lags = nan(length(data_MO_name_ana),2);
for i=1:length(data_MO_name_ana)
    s_cla =  data_SO_4_analysis{i,1};
    for j=2:3
        if ~isempty(data_SO_4_analysis{i,j})
            s_tmp = data_SO_4_analysis{i,j};
            [c,lags] = xcorr(s_cla(1:20000),s_tmp(1:20000),400,'normalized');
            [max_coefficient,loc]=max(c);
            tmp_lag = lags(loc);
            corr_coefficient(i,j-1)=max_coefficient;
            corr_lags(i,j-1) = tmp_lag;
        end
    end
end
figure(),
h=boxplot(corr_coefficient,'Labels',{'CLA&OB','CLA&FI'});
set(h,'LineWidth',1.5);
ylabel({'Cross-correlation','coefficient of SO'})
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2);
ylim([0,1]);
box off;
saveas(gcf,[Fig4b_path,'/Fig4b.tif']);
save([Fig4b_path,'/Fig4b_corr_coefficient.mat'],'corr_coefficient');

% prepare data for python plot
corr_coefficient_in_line = corr_coefficient(:);
x = cell(length(corr_coefficient_in_line),1);
hue = cell(length(corr_coefficient_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(corr_coefficient),1)=cellstr(['CLA&OB']);
x(length(corr_coefficient)+1:length(corr_coefficient)*2,1)=cellstr(['CLA&FI']);
save([Fig4b_path,'/Fig4b_corr_coefficient_4_python.mat'],'x','corr_coefficient_in_line','hue');

%% Fig. 4c 
% SO detect
Fig4c_path = [savepath,'\Fig4c'];
if ~exist(Fig4c_path)
    mkdir(Fig4c_path);
end
load ([savepath,'\data_SO_4_analysis.mat']);
SO_peaks = cell(length(data_MO_name_ana),1);
for idx=1:length(data_MO_name_ana)
    SO_all = data_SO_4_analysis{idx,1};
    if (~isempty(data_SO_4_analysis{idx,2}))&&(~isempty(data_SO_4_analysis{idx,3}))
        SO_all=[SO_all;data_SO_4_analysis{idx,2};data_SO_4_analysis{idx,3}];
    else 
        continue;
    end
    SO_mask = zeros(size(SO_all));
    for i=1:3
        [s_SO_mask_tmp,SO_de] = SO_detect(SO_all(i,:),[0.4,5],0.5,fs);
        SO_mask(i,:)=s_SO_mask_tmp;
    end
    s_SO_mask_overlap = SO_mask(1,:).*SO_mask(2,:).*SO_mask(3,:);
    s_SO_mask_union = double(SO_mask(1,:)|SO_mask(2,:)|SO_mask(3,:));
    s_SO_mask_overlap_union = zeros(size(s_SO_mask_union));
    st = 0;
    ed = 0;
    for i=1:length(s_SO_mask_union)-1
        if(s_SO_mask_union(1,i)==0)&&(s_SO_mask_union(1,i+1)==1)
            st = i;
        end
        if(s_SO_mask_union(1,i)==1)&&(s_SO_mask_union(1,i+1)==0)
            ed = i;
        end
        if st<ed
            if sum(s_SO_mask_overlap(1,st:ed))>0
                s_SO_mask_overlap_union(1,st:ed)=1;
            end
        end
    end
    
    SO_peak_time = [];
    st = 0;
    ed = 0;
    for i=1:length(s_SO_mask_overlap_union)-1
        if(s_SO_mask_overlap_union(1,i)==0)&&(s_SO_mask_overlap_union(1,i+1)==1)
            st = i;
        end
        if(s_SO_mask_overlap_union(1,i)==1)&&(s_SO_mask_overlap_union(1,i+1)==0)
            ed = i;
            [~,loc]=max(SO_all(1,st:ed));
            SO_peak_time = [SO_peak_time;st+loc];
        end
    end
%     figure(),
%     subplot(2,1,1),plot(s_SO_mask_overlap_union);
%     subplot(2,1,2),plot(SO_all(1,:),'b'); hold on;  plot(SO_all(1,:).*s_SO_mask_overlap_union(1,:),'r');   
    tmp_s = {};
    win_length = 1;
    if length(SO_peak_time)>=1
        for i=1:length(SO_peak_time)
            if(SO_peak_time(i)-win_length*fs>=1&&SO_peak_time(i)+win_length*fs<=length(SO_all))
                    tmp_single_SO=SO_all(:,SO_peak_time(i)-win_length*fs:SO_peak_time(i)+win_length*fs);
                    tmp_s = [tmp_s,tmp_single_SO];
            end
        end
    end
    SO_peaks{idx,1} = tmp_s;
end
save([Fig4c_path,'/SO_peaks.mat'],'SO_peaks');
SO_time_lag = nan(length(data_MO_name_ana),2);
for idx =1:length(data_MO_name_ana)
    if ~isempty(SO_peaks{idx,1})
        tmp_ss= SO_peaks{idx,1};
    else
        continue;
    end
    SO_num = length(tmp_ss);
    [x_dim,y_dim]=size(tmp_ss{1});
    tmp_all_s = zeros(x_dim,y_dim,SO_num);
    for i=1:SO_num
        tmp_all_s(:,:,i)=tmp_ss{i};
    end
    tmp_all_s_mean=mean(tmp_all_s,3);
    t_SO_time =1:1:length(tmp_all_s_mean(1,:));
    t_SO_time=t_SO_time*2;
    color_CI = [243/255,164/255,0/255; 209/255 21/255 71/255; 39/255 65/255 146/255];
    color_average = color_CI;
    plot_average(tmp_all_s,t_SO_time,color_average,color_CI)
    saveas(gcf,[Fig4c_path,'/averaged_',subs_name{idx},'.tif']);
    close all;
    [~,loc_cla]=max(tmp_all_s_mean(1,:));
    [~,loc_ob]=max(tmp_all_s_mean(2,:));
    [~,loc_fi]=max(tmp_all_s_mean(3,:));
    SO_time_lag(idx,1) = (loc_ob-loc_cla)*2;
    SO_time_lag(idx,2) = (loc_fi-loc_cla)*2;
end

%% Fig. 4d 
% Group analysis of SO time lag 
Fig4d_path = [savepath,'\Fig4d'];
if ~exist(Fig4d_path)
    mkdir(Fig4d_path);
end
figure(),
h=boxplot(SO_time_lag,'Labels',{'CLA-OB','CLA-FI'});
set(h,'LineWidth',2);
ylabel({'Time lag (ms)'})
set(gca,'FontName','Arial','FontSize',25,'LineWidth',2);
box off;
ylim([-150,150]);
saveas(gcf,[Fig4d_path,'/Fig4d.tif']);
save([Fig4d_path,'/SO_time_lag.mat'],'SO_time_lag');

% prepare data for python plot
SO_time_lag_in_line = SO_time_lag(:);
x = cell(length(SO_time_lag_in_line),1);
hue = cell(length(SO_time_lag_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(SO_time_lag),1)=cellstr(['CLA&OB']);
x(length(SO_time_lag)+1:length(SO_time_lag)*2,1)=cellstr(['CLA&FI']);
save([Fig4d_path,'/Fig4d_SO_time_lag_4_python.mat'],'x','SO_time_lag_in_line','hue');

%% Fig. 4e
Fig4e_path =[savepath, '\Fig4e'];
if ~exist(Fig4e_path)
    mkdir(Fig4e_path);
end
phase_diff_group = cell(length(data_MO_name_ana),2);
phase_diff_average = nan(length(data_MO_name_ana),2);
for idx=1:length(data_MO_name_ana)
    s_cla = data_SO_4_analysis{idx,1};
    for j=1:2
        if ~isempty(data_SO_4_analysis{idx,j+1})
            s_tmp = data_SO_4_analysis{idx,j+1};       
            phase_cla = angle(generalized_phase_vector(s_cla',fs, 0 ))';
            phase_tmp= angle(generalized_phase_vector(s_tmp', fs, 0 ))';
            phase_diff = phase_cla-phase_tmp;
            phase_diff_group{idx,j}=phase_diff;
            if j==1
                z_plot_radar(phase_diff,Fig4e_path,[subs_name{idx},'_CLA_OB.tif'],[209/255,21/255,71/255]); 
            elseif j==2
                z_plot_radar(phase_diff,Fig4e_path,[subs_name{idx},'_CLA_FI.tif'],[39/255,65/255,146/255]);  
            end
            [mu, ~, ~] = circ_mean(phase_diff);
            phase_diff_average(idx,j) =rad2deg(mu);
            close all;
        end
    end
end

phase_diff_group_CLA_OB=[];
phase_diff_group_CLA_FI=[];
for idx=1:length(data_MO_name_ana)
    if ~isempty(phase_diff_group{idx,1})
        phase_diff_group_CLA_OB=[phase_diff_group_CLA_OB,phase_diff_group{idx,1}];
    end
    if ~isempty(phase_diff_group{idx,2})
        phase_diff_group_CLA_FI=[phase_diff_group_CLA_FI,phase_diff_group{idx,2}];
    end
end

z_plot_radar(phase_diff_group_CLA_OB,Fig4e_path,['group_avg_CLA_OB.tif'],[209/255,21/255,71/255]); 
z_plot_radar(phase_diff_group_CLA_FI,Fig4e_path,['group_avg_CLA_FI.tif'],[39/255,65/255,146/255]); 
save([Fig4e_path,'/Fig4e_phase_difference.mat'],'phase_diff_group','phase_diff_group_CLA_OB','phase_diff_group_CLA_FI');

%% Fig. 4f
Fig4f_path =[savepath, '\Fig4f'];
if ~exist(Fig4f_path)
    mkdir(Fig4f_path);
end
figure(),
h=boxplot(phase_diff_average,'Labels',{'CLA-OB','CLA-FI'},'FullFactors','on');
set(h,'LineWidth',1.5);
ylabel({'Phase difference of SO(°)'})
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2);
box off;
ylim([-20,40]);
saveas(gcf,[Fig4f_path,'/Fig4f.tif']);
save([Fig4f_path,'/Fig4f.mat'],'phase_diff_average');

% prepare data for python plot
phase_diff_average_in_line = phase_diff_average(:);
x = cell(length(phase_diff_average_in_line),1);
hue = cell(length(phase_diff_average_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(phase_diff_average),1)=cellstr(['CLA&OB']);
x(length(phase_diff_average)+1:length(phase_diff_average)*2,1)=cellstr(['CLA&FI']);
save([Fig4f_path,'/Fig4f_phase_diff_average_4_python.mat'],'x','phase_diff_average_in_line','hue');

 %% Fig 4g 
Fig4g_path =[savepath,'\Fig4g'];
if ~exist(Fig4g_path)
    mkdir(Fig4g_path);
end
load ([savepath,'\data_SO_4_analysis.mat']);
fs=500;
trial_win = 10;
TE_raw = cell(length(data_MO_name_ana),2);
for idx=1:length(data_MO_name_ana)
    for j=1:2
        if ~isempty(data_SO_4_analysis{idx,j+1})
            s_cla = data_SO_4_analysis{idx,1};
            s_frontal = data_SO_4_analysis{idx,j+1};
            trial_num = floor(length(s_cla)/trial_win/fs);
            TE_single_subject = nan(trial_num,2);
            tau = 6;h = tau;k=4; nsur=10;m = 16;
            parfor k=1:trial_num
                cla_phase_tmp = s_cla(1,(k-1)*trial_win*fs+1:k*trial_win*fs);
                frontal_phase_tmp = s_frontal(1,(k-1)*trial_win*fs+1:k*trial_win*fs);
                [TEM,surM] = TE(cla_phase_tmp',frontal_phase_tmp',h,m,tau,k,nsur);
                TE_single_subject(k,:) = TEM;
            end
            TE_raw{idx,j}=TE_single_subject;
        end
        
    end
    disp(idx);
end
save([Fig4g_path,'/TE_raw.mat'],'TE_raw');

TE_avrage_4_single_subject = nan(length(data_MO_name_ana),4); % column 1~4:CLA-OB,OB-CLA,CLA-FI,FI~CLA

for idx=1:length(data_MO_name_ana)
    for j=1:2
        if ~isempty(TE_raw{idx,j})
           if j==1
                figure();            
                TE_tmp=TE_raw{idx,j};
                h=boxplot(TE_tmp,'Labels',{'CLA-OB','OB-CLA'},'FullFactors','on');
                set(h,'LineWidth',1.5);
                ylabel({'Transfer Entropy'});
                set(gca,'FontName','Arial','FontSize',20,'LineWidth',1.5);
                [~,p]=ttest2(TE_tmp(:,1),TE_tmp(:,2));
                title(num2str(p));
                saveas(gcf,[Fig4g_path,'\',subs_name{idx},'_CLA_OB.tif']);
                close all;
                TE_avrage_4_single_subject(idx,1:2) = mean(TE_raw{idx,j});
           end
           if j==2
                figure();            
                TE_tmp=TE_raw{idx,j};
                h=boxplot(TE_tmp,'Labels',{'CLA-FI','FI-CLA'},'FullFactors','on');
                set(h,'LineWidth',1.5);
                ylabel({'Transfer Entropy'});
                set(gca,'FontName','Arial','FontSize',20,'LineWidth',1.5);
                [~,p]=ttest2(TE_tmp(:,1),TE_tmp(:,2));
                title(num2str(p));
                saveas(gcf,[Fig4g_path,'\',subs_name{idx},'_CLA_FI.tif']);
                close all;
                TE_avrage_4_single_subject(idx,3:4) = mean(TE_raw{idx,j});
           end
        end
    end
end


% prepare data for python plot
% TE_4_python = [TE_avrage_4_single_subject(:,1:2);TE_avrage_4_single_subject(:,3:4)];
x_labl={'CLA-OB','OB-CLA','CLA-FI','FI-CLA'};
x = cell(1,40);
x(1,1:10)= x_labl(1);x(1,11:20)= x_labl(2);x(1,21:30)= x_labl(3);x(1,31:40)= x_labl(4);
hue = cell(1,40);
hue(1,1:20)=cellstr(['CLA&OB']);
hue(1,21:40)=cellstr(['CLA&FI']);
x=x';
hue = hue';
TE_CLA_OB = TE_avrage_4_single_subject(:,1:2);
TE_CLA_FI = TE_avrage_4_single_subject(:,3:4);
x_bar = zeros(10,2);
x_bar(:,1) = 1;
x_bar(:,2) = 2;
x_bar2 = zeros(10,2);
x_bar2(:,1) = 4;
x_bar2(:,2) = 5;
save([Fig4g_path,'/k_TE_4_python.mat'],'TE_CLA_OB','TE_CLA_FI','x_bar','x_bar2');

figure();
subplot(1,2,1),
boxplot(TE_avrage_4_single_subject(:,1:2));
subplot(1,2,2),
boxplot(TE_avrage_4_single_subject(:,3:4));

%% Fig 4i&h
Fig4h_path=[savepath,'./Fig4h'];
if ~exist(Fig4h_path)
    mkdir(Fig4h_path);
end
data_spindle_4_analysis = cell(length(subs_name),3);
data_SO_4_analysis = cell(length(subs_name),3);

for idx = 1:length(subs_name)
    t_win = 100;
    t_s_ana = ana_time{idx,1};
    load(data_MO_name_ana{idx});
    data_ana = SEEG;   
    cfg = [];
    cfg.bpfilter  = 'yes';
    cfg.bpfreq = [10,16]; 
    cfg.channel     = 'all';
    cfg.bpfilttype = 'firws';
    data_spindle = ft_preprocessing(cfg, data_ana); 
    
    cfg = [];
    cfg.bpfilter  = 'yes';
    cfg.bpfreq = [0.3,1]; 
    cfg.channel     = 'all';
    cfg.bpfilttype = 'firws';
    data_SO = ft_preprocessing(cfg, data_ana); 
    
    channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
    for i = 1:length(channel_se)
        if ~isempty(channel_se{i})
            tmp_channel = find(ismember(SEEG.label,channel_se{i}));
            t_s = data_spindle.trial{1}(tmp_channel,t_s_ana*fs:(t_s_ana+t_win)*fs);
            data_spindle_4_analysis(idx,i) =mat2cell(t_s,1);
            t_s_SO = data_SO.trial{1}(tmp_channel,t_s_ana*fs:(t_s_ana+t_win)*fs);
            data_SO_4_analysis(idx,i) =mat2cell(t_s_SO,1);
        end
    end 
end
save([Fig4h_path,'/data_4_analysis.mat'],'data_spindle_4_analysis','data_SO_4_analysis');

load([Fig4h_path,'/data_4_analysis.mat']);
corr_SO_spindle_envelop = nan(10,3);
MI_tort = nan(length(subs_name),3);
for i=1:length(data_SO_4_analysis)
    for j=1:3
        if ~isempty(data_SO_4_analysis{i,j})
            tmp_SO_s = data_SO_4_analysis{i,j}(1,:);
            tmp_SO_phase =   angle(generalized_phase_vector(tmp_SO_s', fs, 0 ))';
            tmp_spindle_s = data_spindle_4_analysis{i,j}(1,:);
            tmp_spindle_envelop = cal_spindle_envelop(tmp_spindle_s);
            corr_SO_spindle_envelop(i,j)=corr(tmp_SO_s',tmp_spindle_envelop');
            nbin=18;
            [MI]=calc_MI_tort(tmp_SO_phase,tmp_spindle_envelop,nbin);
            MI_tort(i,j)=MI;
%             figure(),
%             plot(tmp_SO_phase,tmp_spindle_envelop,'.');
%             xticks([-pi,0,pi]);
%             xlim([-pi,pi]);
%             title(['Subject: ',num2str(i),', ROI: ',num2str(j),', MI:',num2str(MI)]);
%             set(gca,'FontSize',14,'FontName','Arial');
%             saveas(gcf,[Fig4h_path,'/sub',num2str(i),'_ROI',num2str(j),'.tif']);
%             close all;
        end
    end
end
figure();
subplot(1,2,1),boxplot(corr_SO_spindle_envelop);ylim([-1,1]);
subplot(1,2,2),boxplot(MI_tort);ylim([0,0.01]);
saveas(gcf,[Fig4h_path,'/Fig_4i.tif']);

load([Fig4h_path,'/data_4_analysis.mat']);
ts = 500;
te = ts+10000;
i=1;j=1;
tmp_SO_s = data_SO_4_analysis{i,j}(1,ts:te);
tmp_SO_phase =   angle(generalized_phase_vector(tmp_SO_s', fs, 0 ))';
tmp_spindle_s = data_spindle_4_analysis{i,j}(1,ts:te);
tmp_spindle_envelop = cal_spindle_envelop(tmp_spindle_s);

figure(),set(gcf,'position',[100,100,1200,300]);
h4 = cline([1:length(tmp_SO_phase)], tmp_SO_s, [],tmp_SO_phase);
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
hold on;
h1=plot([1:length(tmp_SO_phase)],tmp_spindle_s+tmp_SO_s,'color',[100/256,100/256,100/256]);
set(gca,'child',[h4 h1])
xlim([0,length(tmp_SO_phase)]);
axis off;
saveas(gcf,[Fig4h_path,'/example_signal.tif']);
save([Fig4h_path,'/Fig4h_corr_SO_spindle_envelop.mat'],'corr_SO_spindle_envelop');

% prepare data for python plot
corr_SO_spindle_envelop_in_line = corr_SO_spindle_envelop(:);
MI_in_line = MI_tort(:);
x = cell(length(corr_SO_spindle_envelop_in_line),1);
hue = cell(length(corr_SO_spindle_envelop_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(corr_SO_spindle_envelop),1)=cellstr(['CLA']);
x(length(corr_SO_spindle_envelop)+1:length(corr_SO_spindle_envelop)*2,1)=cellstr(['OB']);
x(length(corr_SO_spindle_envelop)*2+1:length(corr_SO_spindle_envelop)*3,1)=cellstr(['FI']);
save([Fig4h_path,'/Fig4h_4_python.mat'],'x','corr_SO_spindle_envelop_in_line','MI_in_line','hue');
