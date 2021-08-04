function [spindle_mask,spindle_de,s_envelop]= d_spindle_detect(s_ori,win_length,amp_lim,fs)
% This function detect the spindles for time serises.
%INPUTS
% - s_ori        : original spindle wave
% - win_length   : time window length for spindles, generally [0.5,2](0.5s~2s)
% -amp_lim       : amplitude limitation for spindles detection 
% -fs            : sampling rate
%OUTPUTS
% - spindle_mask : mask for the detected spindles
% - spindle_de   : the start time point,end time point and time window of the detected spindles
% - s_envelop    : the envelope of the spindle waves

s_cla_envelop = abs(hilbert(s_ori));
s_temp = s_cla_envelop;
s_temp(s_temp>=amp_lim)=1;
s_temp(s_temp~=1)=0;
spindle_de=[]; 
flag=0;
for i=1:length(s_temp)-1
    if s_temp(i)<=0&&s_temp(i+1)>=0
        flag = 1;
        temp(1)=i;
    end
    if flag ==1&&s_temp(i)>=0&&s_temp(i+1)<=0
        temp(2)=i;
        flag=2;
    end
    if flag==2
        if max(s_cla_envelop(1,temp(1):temp(2)))>=1.5*amp_lim && temp(2)-temp(1)>=win_length(1)*fs &&  temp(2)-temp(1)<=win_length(2)*fs
            spindle_de = [spindle_de;temp];
        end
        temp=[];
        flag=0;
    end
end
% generate spindle mask
spindle_mask =zeros(size(s_ori));
[a,b]=size(spindle_de);
if a~=0
    for i=1:a
        spindle_mask(1,spindle_de(i,1):spindle_de(i,2))=1;
    end
end
s_envelop = s_cla_envelop;


