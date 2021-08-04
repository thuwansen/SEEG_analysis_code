function [s_SW_mask,SW_de]= SO_detect(s_ori,win_length,amp_lim,fs)
% This function detect the SO for time serises.
%INPUTS
% - s_ori        : original spindle wave
% - win_length   : time window length for SO, generally [0.5,2](0.5s~2s)
% - amp_lim      : amplitude limitation for spindles detection 
% -fs            : sampling rate
%OUTPUTS
% - s_SW_mask    : mask for the detected SO
% - SW_de        : the start time point,end time point and peak time of the detected SOs
s_std = std(s_ori);
SW_de=[]; 
flag=0;
for i=1:length(s_ori)-1
    if s_ori(i)<=0&&s_ori(i+1)>=0
        flag = 1;
        temp(1)=i;
    end
    if flag ==1&&s_ori(i)>=0&&s_ori(i+1)<=0
        temp(2)=i;
        flag=2;
    end
    if flag==2
        try
            [pks,~]=findpeaks(s_ori(1,temp(1):temp(2)));
        catch
            pks = 1;
        end
        if temp(2)-temp(1)>win_length(1)*fs&&temp(2)-temp(1)<win_length(2)*fs&&max(s_ori(temp(1):temp(2)))>amp_lim*s_std && length(pks)<=1&&max(s_ori(temp(1):temp(2)))<=10*amp_lim*s_std% 将间距在0.25s~1s之外的去除，峰值小于std的去除
            [~,tmp_peak_loc] =max(s_ori(temp(1):temp(2)));
            temp(3) = temp(1)+tmp_peak_loc;
            SW_de = [SW_de;temp];            
        end
        temp=[];
        flag=0;
    end
end

s_SW_mask = zeros(size(s_ori));
[sw_num,~]=size(SW_de);
if sw_num>=1
    for i=1:sw_num
        s_SW_mask(1,SW_de(i,1):SW_de(i,2))=1; 
    end
end
