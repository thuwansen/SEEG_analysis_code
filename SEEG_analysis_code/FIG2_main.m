% Code for Fig2
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
      
savepath = '.\Fig2';
if ~exist(savepath)
    mkdir(savepath);
end
           
%% Fig 2a top 
% time frequency spectrum
Fig2a_path=[savepath,'\Fig2a'];
if ~exist(Fig2a_path)
    mkdir(Fig2a_path);
end
idx = 4;
load(data_MO_name_ana{idx});
data_ana = SEEG; 
channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
channel_se_idx = zeros(3,1);
for i = 1:length(channel_se)
    channel_se_idx(i,1) = find(ismember(SEEG.label,channel_se{i}));
end
slide_win = 5;
cfg = [];
cfg.length  = slide_win;
cfg.overlap = 2.5/slide_win;
data_segmented = ft_redefinetrial(cfg, SEEG);
cfg = [];
cfg.method     = 'mtmfft'
cfg.taper      = 'hanning'
cfg.foilim     = [0 30];
cfg.keeptrials = 'yes'
freq_segmented = ft_freqanalysis(cfg, data_segmented);
final_wave = freq_segmented.powspctrm;
temp2 =squeeze(final_wave(:,channel_se_idx(1),:));
close all;figure(),set(gcf,'position',[100 100 1200 200]);

imagesc(400:1526,freq_segmented.freq,log10(temp2(400:1526,:)'),[-4,4]);
axis xy
axis off;
colormap jet;
saveas(gcf,[Fig2a_path,'/Fig2a_top.png']);

hold on;plot([400,1526],[10,10],'k--','LineWidth',1);
hold on;plot([400,1526],[16,16],'k--','LineWidth',1);
hold on;plot([400,1526],[1,1],'k--','LineWidth',1);
saveas(gcf,[Fig2a_path,'/Fig2a_top_band.png']);

%% Fig 2a bottom
% raw signal
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq = [0.3,1]; 
cfg.channel     = 'all';
cfg.bpfilttype = 'firws';
data_SO = ft_preprocessing(cfg, data_ana); 
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq = [10,16]; 
cfg.channel     = 'all';
cfg.bpfilttype = 'firws';
data_spindle = ft_preprocessing(cfg, data_ana); 
t_s = 2140;
t_win =50;
b_signal_ori = SEEG.trial{1}(channel_se_idx,fs*t_s+1:fs*(t_s+t_win));
b_signal_SO = data_SO.trial{1}(channel_se_idx,fs*t_s+1:fs*(t_s+t_win));
b_signal_spindle = data_spindle.trial{1}(channel_se_idx,fs*t_s+1:fs*(t_s+t_win));
figure(),set(gcf,'position',[100,100,1200,400]);
plot(b_signal_ori(1,:),'color',color_group(1,:),'LineWidth',1.5);
hold on;plot(b_signal_SO(1,:)-400,'color',color_group(1,:),'LineWidth',1.5);
hold on;plot(b_signal_spindle(1,:)*2-800,'color',color_group(1,:),'LineWidth',1.5);
saveas(gcf,[Fig2a_path,'/Fig2a_bottom.tif']);
set(gca,'position',[0 0 1 1])
axis off;
saveas(gcf,[Fig2a_path,'/Fig2a_bottom_axis_off.tif']);

%% Fig 2b, d
Fig2b_path=[savepath,'\Fig2b_d'];
if ~exist(Fig2b_path)
    mkdir(Fig2b_path);
end
channel_se = {cla_selected{idx};ob_selected{idx};fi_selected{idx}};
channel_se_idx = zeros(3,1);
for i = 1:length(channel_se)
    channel_se_idx(i,1) = find(ismember(SEEG.label,channel_se{i}));
end
slide_win = 5;
cfg = [];
cfg.length  = slide_win;
cfg.overlap = 0.8;
data_segmented = ft_redefinetrial(cfg, SEEG);
cfg = [];
cfg.method     = 'mtmfft'
cfg.taper      = 'hanning'
cfg.foilim     = [0.5 200];
cfg.keeptrials = 'yes'
freq_segmented = ft_freqanalysis(cfg, data_segmented);
save([Fig2b_path,'/freq_segmented_1s.mat'],'freq_segmented','-v7.3');

pre_time =ana_time{idx};
post_time=con_time{idx};
epoch_num=100;
name = {'Fig2b_cla.tif', 'Fig2d_ob.tif','Fig2d_fi.tif'};

for i=1:3
    channel = channel_se_idx(i);
    ana_amp_squared = squeeze(freq_segmented.powspctrm(pre_time:pre_time+epoch_num,channel,:));%uV*uV
    ana_power = 20*log10(ana_amp_squared); % dBuV
    ana_fre =freq_segmented.freq; % Hz
    con_amp_squared = squeeze(freq_segmented.powspctrm(post_time:post_time+epoch_num,channel,:));%uV*uV
    con_power = 20*log10(con_amp_squared); % dBuV
    con_fre =freq_segmented.freq; % Hz
    figure(),
    set(gcf,'position',[100,100,250,300])
    [mu,sigma,muci,sigmaci]=normfit(ana_power);
    semilogx(ana_fre, mu,'-r');
    xconf = [ana_fre ana_fre(end:-1:1)];
    yconf = [muci(1,:) muci(2,end:-1:1)];
    hold on;
    p = fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColor为填充颜色，EdgeColor为边框颜�?
    hold on;
    a=semilogx(ana_fre, mu,'-r');

    [mu,sigma,muci,sigmaci]=normfit(con_power);
    semilogx(ana_fre, mu,'-k');
    xconf = [ana_fre ana_fre(end:-1:1)];
    yconf = [muci(1,:) muci(2,end:-1:1)];
    hold on;
    p = fill(xconf,yconf,'','FaceColor',[0.8 0.8 0.8],'EdgeColor','none');%FaceColor为填充颜色，EdgeColor为边框颜�?
    hold on;
    a=semilogx(ana_fre, mu,'-k');   
    ylabel('Power Density (dB/Hz)');
    xlabel('Frequency (Hz)');
    ylim([-100,100])
    xlim([0.7,100]);
    set(gca,'FontSize',12,'FontName','Arial','Linewidth',1.5);
    saveas(gcf,[Fig2b_path,'/',name{i}]);
    close all;
end

%% Fig 2c
Fig2c_path=[savepath,'\Fig2c'];
if ~exist(Fig2c_path)
    mkdir(Fig2c_path);
end

temp_OB =squeeze(final_wave(:,channel_se_idx(2),:));
close all;figure(),
imagesc(400:1526,freq_segmented.freq,log10(temp_OB(400:1526,:)'),[-4,4]);
axis xy
set(gcf,'position',[100 100 1200 200]);
axis off;
colormap jet;
hold on;plot([400,1526],[10,10],'k--','LineWidth',1);
hold on;plot([400,1526],[16,16],'k--','LineWidth',1);
hold on;plot([400,1526],[1,1],'k--','LineWidth',1);
saveas(gcf,[Fig2c_path,'\Fig2c_OB.png']);

temp_FI =squeeze(final_wave(:,channel_se_idx(3),:));
close all;figure(),
imagesc(400:1526,freq_segmented.freq,log10(temp_FI(400:1526,:)'),[-4,4]);
axis xy
set(gcf,'position',[100 100 1200 200]);
axis off;
colormap jet;
hold on;plot([400,1526],[10,10],'k--','LineWidth',1);
hold on;plot([400,1526],[16,16],'k--','LineWidth',1);
hold on;plot([400,1526],[1,1],'k--','LineWidth',1);
saveas(gcf,[Fig2c_path,'/Fig2c_FI.png']);

%% Fig 2e
Fig2e_path=[savepath,'\Fig2e'];
if ~exist(Fig2e_path)
    mkdir(Fig2e_path);
end
tmp_time_fre_log10_CLA = log10(temp2(400:1526,:)');
tmp_time_fre_log10_OB= log10(temp_OB(400:1526,:)');
tmp_time_fre_log10_FI= log10(temp_FI(400:1526,:)');
mean_fre_so_CLA = mean(tmp_time_fre_log10_CLA(1:6,:),1);
mean_fre_so_nor_CLA = mean_fre_so_CLA-mean(mean_fre_so_CLA(1000:1100));
mean_fre_spindle_CLA = mean(tmp_time_fre_log10_CLA(51:81,:),1);
mean_fre_spindle_nor_CLA = mean_fre_spindle_CLA-mean(mean_fre_spindle_CLA(1000:1100)); 
mean_fre_so_OB = mean(tmp_time_fre_log10_OB(1:6,:),1);
mean_fre_so_nor_OB = mean_fre_so_OB-mean(mean_fre_so_OB(1000:1100));
mean_fre_spindle_OB = mean(tmp_time_fre_log10_OB(51:81,:),1);
mean_fre_spindle_nor_OB = mean_fre_spindle_OB-mean(mean_fre_spindle_OB(1000:1100)); 
mean_fre_so_FI = mean(tmp_time_fre_log10_FI(1:6,:),1);
mean_fre_so_nor_FI = mean_fre_so_FI-mean(mean_fre_so_FI(1000:1100));
mean_fre_spindle_FI = mean(tmp_time_fre_log10_FI(51:81,:),1);
mean_fre_spindle_nor_FI = mean_fre_spindle_FI-mean(mean_fre_spindle_FI(1000:1100)); 

figure(),set(gcf,'position',[100 100 1200 200]);
t_time = 400:1:1526;
t_time = (t_time-400)*2.5;
hold on; c = smoothdata(mean_fre_so_nor_CLA,'gaussian',3);plot(t_time,c,'-','LineWidth',1.5,'color',color_group(1,:));
hold on; c = smoothdata(mean_fre_so_nor_OB,'gaussian',3)-1.7;plot(t_time,c,'-','LineWidth',1.5,'color',color_group(2,:));
hold on; c = smoothdata(mean_fre_so_nor_FI,'gaussian',3)-3.4;plot(t_time,c,'-','LineWidth',1.5,'color',color_group(3,:));
xlim([0,max(t_time)]);
ylim([-4.5,2]);
axis off;
saveas(gcf,[Fig2e_path,'\Power_dynamic_SO.png']);

figure(),set(gcf,'position',[100 100 1200 200]);
t_time = 400:1:1526;
t_time = (t_time-400)*2.5;
hold on; c = smoothdata(mean_fre_spindle_CLA,'gaussian',3);plot(t_time,c,'-','LineWidth',1.5,'color',color_group(1,:));
hold on; c = smoothdata(mean_fre_spindle_nor_OB,'gaussian',3)-0.5;plot(t_time,c,'-','LineWidth',1.5,'color',color_group(2,:));
hold on; c = smoothdata(mean_fre_spindle_nor_FI,'gaussian',3)-1;plot(t_time,c,'-','LineWidth',1.5,'color',color_group(3,:));
xlim([0,max(t_time)]);
axis off;
saveas(gcf,[Fig2e_path,'/Power_dynamic_spindle.png']);




