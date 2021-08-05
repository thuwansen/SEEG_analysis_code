function data_slice = k_func_slice_data(data,start_t,t_win,fs)
% This function slice
% This function slice the data to specific length of time.
%INPUTS
% - data         : data in fiedtrip method
% - start_t      : start time 
% - t_win        : time length
% - fs           : sampling rate
%OUTPUTS
% - data_slice   : sliced data

slice_time = 0:1/fs:t_win;
cfg = [];
cfg.latency = [start_t(1),start_t(1)+t_win];
data_tmp = ft_selectdata(cfg, data);
data_tmp.time{1,1}=slice_time;

if length(start_t)>=2
    for i=2:length(start_t)
        cfg = [];
        cfg.latency = [start_t(i),start_t(i)+t_win];
        temp_data = ft_selectdata(cfg, data);
        data_tmp.trial = [data_tmp.trial,temp_data.trial];
        data_tmp.time = [data_tmp.time,slice_time];
        data_tmp.sampleinfo =[data_tmp.sampleinfo; temp_data.sampleinfo];
    end
end
data_slice = data_tmp;