function s_ori = z_cut_data (data,t_s,t_win,fs,channel_name)
% This function cut the signal for given original data, channel name, start and end time, sampling rate.
%INPUTS
% - data         : original data in fieldtrip format
% - t_win        : time length
% - t_s          : start time 
% - fs           : sampling rate
% channel_name   : name of the target channel.
%OUTPUTS
% - s_ori        : target signal
%Example
%  s_ori=z_cut_data (data_MO,[3201,3211,3221],[10,10,10],500,'PI''2');
s_ori_idx = find(ismember(data.label,channel_name));
s_ori=[];
for i = 1:length(t_s)
    s_ori_temp = data.trial{1,1}(s_ori_idx,t_s(i)*fs:(t_s(i)+t_win(i))*fs);
    s_ori= [s_ori,s_ori_temp];
end

