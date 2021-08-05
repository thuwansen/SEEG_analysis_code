% Code for Fig3 and Supplementary Fig 2
clc,clear,close all;
addpath('D:\software\toolbox\fieldtrip-20191213');
addpath(genpath('D:\software\toolbox\nature_walk_wave'));
addpath('D:\software\toolbox\BrainNetViewer_20191031')
addpath(genpath('D:\software\PTE-variants-master'));
code_path = 'D:\SEEG_preoject\BIDS_data\code\research_functions';
addpath(genpath([code_path,'\research_functions']));
run([code_path,'\Global_variable_define_LOC2ROC.m']);
fs =500;
% Define the colors for CLA, OB and FI.
color_group = [243/255,164/255,0/255;
               209/255,21/255,71/255;
               39/255,65/255,146/255];
      
savepath = '.\Fig3';
if ~exist(savepath)
    mkdir(savepath);
end
           
           
%% Fig.3a
Fig3a_path=[savepath, '.\Fig3a'];
if ~exist(Fig3a_path)
    mkdir(Fig3a_path);
end
idx = 2;
load(data_MO_name_ana{idx});
data_ana = SEEG;   
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq = [10,16]; 
cfg.channel     = 'all';
cfg.bpfilttype = 'firws';
data_spindle = ft_preprocessing(cfg, data_ana); 
channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
channel_se_idx = zeros(3,1);
for i = 1:length(channel_se)
    channel_se_idx(i,1) = find(ismember(SEEG.label,channel_se{i}));
end
t_s_ana_fig_f = 300;
data_ana_se = data_spindle.trial{1}(channel_se_idx,fs*t_s_ana_fig_f:fs*(t_s_ana_fig_f+10));
data_ana_se_spindle_std = std(data_ana_se,1,2);
data_ana_se = data_ana_se./repmat(data_ana_se_spindle_std,[1,length(data_ana_se)]);

for i=1:3
     figure; set(gcf,'position',[100,100,900,150]);
     set(gca,'position',[0 0 1 1]);
     plot(data_ana_se(i,:),'LineWidth',2,'color',color_group(i,:));xlim([0,length(data_ana_se)]); ylim([-5,5]);
     axis off;
     saveas(gcf,[Fig3a_path,'/spindle_raw_',num2str(i),'.tif']);
end

% time-frequency
dt =1/fs;timestart=0;timeend=10;
t=(0:(timeend-timestart)/dt-1)*dt+timestart;L=length(t);
for i=1:3
    z=data_ana_se(i,:);
    z2=wextend(1,'sym',z,round(length(z)/2));
    wlen=512;hop=1;
    z2=wkeep1(z2,L+1*wlen);
    h=hamming(wlen);
    f=0.3:0.2:25;
    [tfr2,f,t2]=spectrogram(z2,h,wlen-hop,f,fs);
    tfr2=tfr2*2/wlen*2;
    figure(),set(gcf,'position',[100,100,900,100]);
    set(gca,'position',[0 0 1 1]);
    imagesc(t2+timestart-wlen/fs/2,f,abs(tfr2),[0,2]);
    colormap jet;axis off;axis xy;
    saveas(gcf,[Fig3a_path,'/spindle_time_frequency',num2str(i),'.tif']);
end

%% Fig3e
% plot exam signal with phase
Fig3e_path=[savepath, '.\Fig3e'];
if ~exist(Fig3e_path)
    mkdir(Fig3e_path);
end
st = 2020;
ed = st+150;
figure; set(gcf,'position',[100,100,800,350]);set(gca,'position',[0 0 1 1]);
for i=1:3
     plot(data_ana_se(i,:)-10*i,'LineWidth',2,'color',color_group(i,:));xlim([st,ed]);
     tmp_phase =   angle(generalized_phase_vector(data_ana_se(i,:)', fs, 0 ))';
     h4 = cline([1:length(tmp_phase)], data_ana_se(i,:)-10*i, [],tmp_phase);
     set( h4, 'linestyle', '-', 'linewidth', 3); 
     map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
     hold on;
end
grid off; axis off;
saveas(gcf,[Fig3e_path,'/spindle_phase.tif']);

%% Prepare data for group analysis 
data_spindle_4_analysis = cell(length(subs_name),3);
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
    channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
    for i = 1:length(channel_se)
        if ~isempty(channel_se{i})
            tmp_channel = find(ismember(SEEG.label,channel_se{i}));
            t_s = data_spindle.trial{1}(tmp_channel,t_s_ana*fs:(t_s_ana+t_win)*fs);
            data_spindle_4_analysis(idx,i) =mat2cell(t_s,1);
        end
    end 
end
save([savepath,'\data_spindle_4_analysis.mat'],'data_spindle_4_analysis');

%% Fig.3b
% correlation of spindles envelop
Fig3b_path=[savepath,'\Fig3b'];
if ~exist(Fig3b_path)
    mkdir(Fig3b_path);
end
load([savepath,'\data_spindle_4_analysis.mat']);
spindle_envelop_corr = nan(length(data_spindle_4_analysis),2);
spindle_envelop_corr_loc = nan(length(data_spindle_4_analysis),2);

t_win = 100;
for idx=1:length(data_spindle_4_analysis)
    tmp_CLA = data_spindle_4_analysis{idx,1}(1,1:fs*t_win);
    tmp_CLA_envelop = abs(hilbert(tmp_CLA'))';
    for j=1:2
        if ~isempty(data_spindle_4_analysis{idx,j+1})
            tmp_frontal= data_spindle_4_analysis{idx,j+1}(1,1:fs*t_win);
            tmp_frontal_envelop = abs(hilbert(tmp_frontal'))';
            [c,lag]=xcorr(tmp_CLA_envelop',tmp_frontal_envelop',100,'normalized');
            spindle_envelop_corr(idx,j)=max(c);
            [~,loc]=max(c);
            spindle_envelop_corr_loc(idx,j)=lag(loc)*2;
        end
    end
end

figure(), set(gcf,'position',[100,100,500,600]);
h=boxplot(spindle_envelop_corr,'Labels',{'CLA&OB','CLA&FI'});
set(h,'LineWidth',1.5);
% ylabel({'Cross-correlation coefficient of spindles envelop'})
ylabel({'Max cross-correlation coefficient'})
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2);
ylim([0,1]);
box off;
saveas(gcf,[Fig3b_path,'/Fig3b.tif']);
save([Fig3b_path,'/Fig3b_crosscorrelation_coefficient.mat'],'spindle_envelop_corr');

% prepare data for python plot
spindle_envelop_corr_in_line = spindle_envelop_corr(:);
x = cell(length(spindle_envelop_corr_in_line),1);
hue = cell(length(spindle_envelop_corr_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(spindle_envelop_corr),1)=cellstr(['CLA&OB']);
x(length(spindle_envelop_corr)+1:length(spindle_envelop_corr)*2,1)=cellstr(['CLA&FI']);
save([Fig3b_path,'/Fig3b_4_python.mat'],'x','spindle_envelop_corr_in_line','hue');

%% Fig.3c
% cross-correlation of spindles
Fig3c_path=[savepath,'/Fig3c'];
if ~exist(Fig3c_path)
    mkdir(Fig3c_path);
end
load([savepath,'\data_spindle_4_analysis.mat']);
idx=2;
s_4_fig3c = [data_spindle_4_analysis{idx,1};data_spindle_4_analysis{idx,2};data_spindle_4_analysis{idx,3}];
spindle_4_fig3c = s_4_fig3c(:,11900:12400);
[c_cla_ob,lags_cla_ob]=xcorr(spindle_4_fig3c(1,:),spindle_4_fig3c(2,:),250,'normalized');
[c_cla_fi,lags_cla_fi]=xcorr(spindle_4_fig3c(1,:),spindle_4_fig3c(3,:),250,'normalized');
figure(), set(gcf,'position',[100,100,800,500]);
plot(lags_cla_ob*2,c_cla_ob,'color',color_group(2,:),'LineWidth',1.5);
hold on;
plot(lags_cla_fi*2,c_cla_fi,'color',color_group(3,:),'LineWidth',1.5);
ylim([-1,1]);
xlim([-500,500]);
xticks([-500,-250,0,250,500]);
% xlabel('Time (ms)');
box off;
set(gca,'FontName','Arial','FontSize',25,'LineWidth',2);
saveas(gcf,[Fig3c_path,'/Fig3_up.tif']);
[~,loc1] = max(c_cla_ob); lags_cla_ob(loc1)
[~,loc2] = max(c_cla_fi); lags_cla_ob(loc2)

figure(), set(gcf,'position',[100,100,800,500]);
set(gca,'position',[0 0 1 1])
tmp_win=20;
plot(lags_cla_ob(1,250-tmp_win:250+tmp_win)*2,c_cla_ob(1,250-tmp_win:250+tmp_win),'color',color_group(2,:),'LineWidth',2);
hold on;
plot(lags_cla_fi(1,250-tmp_win:250+tmp_win)*2,c_cla_fi(1,250-tmp_win:250+tmp_win),'color',color_group(3,:),'LineWidth',2);
ylim([-1,1]);
xlim([-tmp_win*2,tmp_win*2]);
% xticks([-500,-250,0,250,500]);
% xlabel('Time (ms)');
axis off;
set(gca,'FontName','Arial','FontSize',25,'LineWidth',2);
saveas(gcf,[Fig3c_path,'/Fig3_bottom.tif']);
save([Fig3c_path,'/Time_lags.mat'],'lags_cla_fi','lags_cla_ob');

%% Fig3d
% cross-correlation of spindles
Fig3d_path=[savepath,'\Fig3d'];
if ~exist(Fig3d_path)
    mkdir(Fig3d_path);
end
load([savepath,'\data_spindle_4_analysis.mat']);
spindle_cross_corr = nan(length(data_spindle_4_analysis),2);
spindle_max_corr_loc = nan(length(data_spindle_4_analysis),2);
% time lags calculated for every 10s signal
% outliers were removed and mean of time lags calculated 
spindle_cross_corr = cell(length(data_spindle_4_analysis),2);
spindle_max_corr_loc = cell(length(data_spindle_4_analysis),2);
t_win = 100;
win_num = 10;
for idx=1:length(data_spindle_4_analysis)
    tmp_CLA = data_spindle_4_analysis{idx,1}(1,1:fs*t_win);
    for j=1:2
        if ~isempty(data_spindle_4_analysis{idx,j+1})
            tmp_frontal= data_spindle_4_analysis{idx,j+1}(1,1:fs*t_win);
            
            cross_corr_group = zeros(10,1);
            max_corr_loc = zeros(10,1);
            for i=1:win_num
                [c,lag]=xcorr(tmp_CLA(1,1+(i-1)*fs*10:i*fs*10)',tmp_frontal(1,1+(i-1)*fs*10:i*fs*10)',50,'normalized');
                cross_corr_group(i,1)=max(c);
                [~,loc]=max(c);
                max_corr_loc(i,1)=lag(loc)*2;
            end
            spindle_cross_corr{idx,j}=cross_corr_group;
            spindle_max_corr_loc{idx,j} = max_corr_loc;
        end
    end
end

spindle_time_lags = nan(length(data_spindle_4_analysis),2);
for idx = 1:length(data_spindle_4_analysis)
    for j=1:2
          if ~isempty(spindle_max_corr_loc{idx,j})  
                tmp_timelags = spindle_max_corr_loc{idx,j};
                B = rmoutliers(tmp_timelags); % remove outliers
                spindle_time_lags(idx,j)=mean(B);
          end
    end
end


figure(),set(gcf,'position',[100,100,500,600]);
h=boxplot(spindle_time_lags,'Labels',{'CLA&OB','CLA&FI'});
set(h,'LineWidth',1.5);
ylabel({'Time lags of spindles'})
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2);
ylim([-5,15]);
box off;
saveas(gcf,[Fig3d_path,'/Fig3d.tif']);
save([Fig3d_path,'/Fig3d_time_lags.mat'],'spindle_time_lags');

% prepare data for python plot
spindle_time_lags_in_line = spindle_time_lags(:);
x = cell(length(spindle_time_lags_in_line),1);
hue = cell(length(spindle_time_lags_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(spindle_time_lags),1)=cellstr(['CLA&OB']);
x(length(spindle_time_lags)+1:length(spindle_time_lags)*2,1)=cellstr(['CLA&FI']);
save([Fig3d_path,'/Fig3d_4_python.mat'],'x','spindle_time_lags_in_line','hue');

%% Fig. 3e right part and Supplementary Fig 2b
% phase difference
load([savepath,'\data_spindle_4_analysis.mat']);
phase_diff_group = cell(length(data_spindle_4_analysis),2);
phase_diff_average = nan(length(data_spindle_4_analysis),2);
for idx=1:length(data_spindle_4_analysis)
    s_cla = data_spindle_4_analysis{idx,1};
    for j=1:2
        if ~isempty(data_spindle_4_analysis{idx,j+1})
            s_tmp = data_spindle_4_analysis{idx,j+1};       
            phase_cla = angle(generalized_phase_vector(s_cla',fs, 0 ))';
            phase_tmp= angle(generalized_phase_vector(s_tmp', fs, 0 ))';
            phase_diff = phase_cla-phase_tmp;
            phase_diff_group{idx,j}=phase_diff;
            if j==1
                z_plot_radar(phase_diff,Fig3e_path,[subs_name{idx},'_CLA_OB.tif'],[209/255,21/255,71/255]); 
            elseif j==2
                z_plot_radar(phase_diff,Fig3e_path,[subs_name{idx},'_CLA_FI.tif'],[39/255,65/255,146/255]);  
            end
            [mu, ~, ~] = circ_mean(phase_diff);
            phase_diff_average(idx,j) =rad2deg(mu);
            close all;
        end
    end
end
phase_diff_group_CLA_OB=[];
phase_diff_group_CLA_FI=[];
for idx=1:length(data_spindle_4_analysis)
    if ~isempty(phase_diff_group{idx,1})
        phase_diff_group_CLA_OB=[phase_diff_group_CLA_OB,phase_diff_group{idx,1}];
    end
    if ~isempty(phase_diff_group{idx,2})
        phase_diff_group_CLA_FI=[phase_diff_group_CLA_FI,phase_diff_group{idx,2}];
    end
end

z_plot_radar(phase_diff_group_CLA_OB,Fig3e_path,['group_avg_CLA_OB.tif'],[209/255,21/255,71/255]); 
z_plot_radar(phase_diff_group_CLA_FI,Fig3e_path,['group_avg_CLA_FI.tif'],[39/255,65/255,146/255]); 
save([Fig3e_path,'/Fig3e_phase_difference.mat'],'phase_diff_group','phase_diff_group_CLA_OB','phase_diff_group_CLA_FI');

%% Fig. 3f
% phase difference
Fig3f_path=[savepath,'\Fig3f'];
if ~exist(Fig3f_path)
    mkdir(Fig3f_path);
end
figure(),set(gcf,'position',[100,100,500,600]);
h=boxplot(phase_diff_average,'Labels',{'CLA-OB','CLA-FI'},'FullFactors','on');
set(h,'LineWidth',1.5);
ylabel({'Phase difference of spindles(°)'})
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2);
box off;
ylim([-60,20]);
saveas(gcf,[Fig3f_path,'/Fig3f.tif']);
save([Fig3f_path,'/Fig3f.mat'],'phase_diff_average');
% prepare data for python plot
phase_diff_average_in_line = phase_diff_average(:);
x = cell(length(phase_diff_average_in_line),1);
hue = cell(length(phase_diff_average_in_line),1);
hue(:,1)=cellstr(['CLA&OB']);
x(1:length(phase_diff_average),1)=cellstr(['CLA&OB']);
x(length(phase_diff_average)+1:length(phase_diff_average)*2,1)=cellstr(['CLA&FI']);
save([Fig3f_path,'/Fig3f_phase_diff_average_4_python.mat'],'x','phase_diff_average_in_line','hue');

%% Fig. 3g
% phase difference
Fig3g_path=[savepath,'\Fig3g'];
if ~exist(Fig3g_path)
    mkdir(Fig3g_path);
end

load ([savepath,'\data_spindle_4_analysis.mat']);
trial_win = 10;
TE_raw = cell(length(subs_name),2);
for idx=1:length(subs_name)
    for j=1:2
        if ~isempty(data_spindle_4_analysis{idx,j+1})
            s_cla = data_spindle_4_analysis{idx,1};
            s_frontal = data_spindle_4_analysis{idx,j+1};
            trial_num = floor(length(s_cla)/trial_win/fs);
            TE_single_subject = nan(trial_num,2);
            tau = 10;h = tau;k=4; nsur=10;m = 16;
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
save([Fig3g_path,'/TE_raw.mat'],'TE_raw');
TE_avrage_4_single_subject = nan(length(subs_name),4); % column 1~4:CLA-OB,OB-CLA,CLA-FI,FI~CLA
TE_raw_group=nan(length(subs_name),4);
for i=1:length(TE_raw)
    if ~isempty(TE_raw{i,1})
        tmp_mean=mean(TE_raw{i,1});
        TE_raw_group(i,1:2)=tmp_mean;
    end
    if ~isempty(TE_raw{i,2})
        tmp_mean=mean(TE_raw{i,2});
        TE_raw_group(i,3:4)=tmp_mean;
    end    
end

% prepare data for python plot
% TE_4_python = [TE_avrage_4_single_subject(:,1:2);TE_avrage_4_single_subject(:,3:4)];
x_labl={'CLA-OB','OB-CLA','CLA-FI','FI-CLA'};
x = cell(1,36);
x(1,1:9)= x_labl(1);x(1,10:18)= x_labl(2);x(1,19:27)= x_labl(3);x(1,28:36)= x_labl(4);
hue = cell(1,40);
hue(1,1:18)=cellstr(['CLA&OB']);
hue(1,19:36)=cellstr(['CLA&FI']);
x=x';
hue = hue';
TE_CLA_OB = TE_raw_group(:,1:2);
TE_CLA_FI = TE_raw_group(:,3:4);
x_bar = zeros(9,2);
x_bar(:,1) = 1;
x_bar(:,2) = 2;
x_bar2 = zeros(9,2);
x_bar2(:,1) = 4;
x_bar2(:,2) = 5;
save([Fig3g_path,'/k_TE_4_python.mat'],'TE_CLA_OB','TE_CLA_FI','x_bar','x_bar2');

%% Fig. 3h
Fig3h_path=[savepath,'\Fig3h'];
if ~exist(Fig3h_path)
    mkdir(Fig3h_path);
end
idx = 5;
load(SEEG_coordinate_path{idx});
channel_CLA = cla_selected{idx};
channel_OB = {'OB1';'OB2';'OB3';'OB4';'OB5';'OB6';'OB7';'OB8';
              'OB9';'OB10';'OB11';'OB12';'OB13';'OB14';'OB15'};
channel_CLA_idx = find(ismember(SEEG_COORDINATE.SEEG_leadlabel,channel_CLA));
channel_OB_idx = nan(size(channel_OB));
for i = 1:length(channel_OB)
    channel_OB_idx(i,1) = find(ismember(SEEG_COORDINATE.SEEG_leadlabel,channel_OB{i}));
end
CLA_coordinate = SEEG_COORDINATE.MNI_coordinate(channel_CLA_idx,:);
OB_coordinate = SEEG_COORDINATE.MNI_coordinate(channel_OB_idx,:);
CLA_OB_coordinate = zeros(length(channel_CLA)+length(channel_OB),6);
CLA_OB_coordinate(:,4:6)=1;
CLA_OB_coordinate(1:length(channel_OB),1:3)=OB_coordinate;
CLA_OB_coordinate(length(channel_OB)+1:end,1:3)=CLA_coordinate;
CLA_OB_coordinate(length(channel_OB)+1:end,4)=2;
save([Fig3h_path,'/CLA_OB.node'], 'CLA_OB_coordinate','-ascii');

%% Fig 3i
Fig3i_path=[savepath,'\Fig3i'];
if ~exist(Fig3i_path)
    mkdir(Fig3i_path);
end

idx=5;
t_win = 100;
t_s_ana = ana_time{idx,1};
load(data_MO_name_ana{idx});
load(SEEG_coordinate_path{idx});
data_ana = SEEG;   
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq = [10,16]; 
cfg.channel     = 'all';
cfg.bpfilttype = 'firws';
data_spindle = ft_preprocessing(cfg, data_ana); 
channel_CLA = cla_selected{idx};
channel_OB = {'OB1';'OB2';'OB3';'OB4';'OB5';'OB6';'OB7';'OB8';
              'OB9';'OB10';'OB11';'OB12';'OB13';'OB14';'OB15';'OB16';};
channel_CLA_idx = find(ismember(SEEG_COORDINATE.SEEG_leadlabel,channel_CLA));
channel_OB_idx = nan(size(channel_OB));
for i = 1:length(channel_OB)
    channel_OB_idx(i,1) = find(ismember(SEEG_COORDINATE.SEEG_leadlabel,channel_OB{i}));
end

% cal the distance between OB and CLA
CLA_OB_distance = nan(size(channel_OB_idx));
CLA_coordinate = SEEG_COORDINATE.MNI_coordinate(channel_CLA_idx,:);
for i=1:length(channel_OB_idx)
    OB_coordinate = SEEG_COORDINATE.MNI_coordinate(channel_OB_idx(i,1),:);
    CLA_OB_distance(i,1) = norm(CLA_coordinate-OB_coordinate);
end
% cal the phase difference
s_CLA = data_spindle.trial{1}(channel_CLA_idx,t_s_ana*fs:(t_s_ana+t_win)*fs);
s_OB= data_spindle.trial{1}(channel_OB_idx,t_s_ana*fs:(t_s_ana+t_win)*fs);
phase_diff_fig_h = cell(length(channel_OB),1);
phase_cla = angle(generalized_phase_vector(s_CLA',fs, 0 ))';
for j=1:length(channel_OB)
        phase_ob= angle(generalized_phase_vector(s_OB(j,:)', fs, 0 ))';
        phase_diff_mean = zeros(10,1);
        for k=1:10
            phase_diff = phase_cla(1,1+(k-1)*10*fs:k*10*fs)-phase_ob(1,1+(k-1)*10*fs:k*10*fs);
            [mu, ~, ~] = circ_mean(phase_diff);
            phase_diff_mean(k,1)=rad2deg(mu);
        end
        phase_diff_fig_h{j,1}=phase_diff_mean;
end

phase_diff_fig_i_mat =reshape(cell2mat(phase_diff_fig_h),10,16);
phase_diff_fig_i_mat_median = median(phase_diff_fig_i_mat);

figure(); set(gcf,'position',[100,100,1500,1200]);
subplot(2,1,1),
h=boxplot(-phase_diff_fig_i_mat(:,1:15),'Labels',channel_OB(1:15,1));
ylabel('Phase difference');
set(h,'LineWidth',1.5);
set(gca,'FontName','Arial','FontSize',20,'LineWidth',1.5);
box off;
subplot(2,1,2),
plot(CLA_OB_distance(1:15,1),'-ro','LineWidth',1.5);
% xticklabels(channel_OB(1:15,1));
saveas(gcf,[Fig3i_path,'/Fig3i.tif']);
save([Fig3i_path,'/Fig3i_phase_diff.mat'],'phase_diff_fig_i_mat','CLA_OB_distance');

% save data for python plot
x_lable=channel_OB;
x = cell(16*10,1);
hue = cell(16*10,1);
phase_diff_line = reshape(phase_diff_fig_i_mat,16*10,1);
hue(:,1)=cellstr(['OB1']);
for i=1:16
    x(1+(i-1)*10:i*10,1)=cellstr(x_lable{i});
end
save([Fig3i_path,'/Fig3i_4_python_plot.mat'],'x','phase_diff_line','hue','phase_diff_fig_i_mat_median','CLA_OB_distance');








