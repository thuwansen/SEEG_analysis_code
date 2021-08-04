function s_en=cal_spindle_envelop(s_ori)
% This function calculate the envelope for time serises.
%INPUTS
% - s_ori        : time series, size: n x sample
%OUTPUTS
% - s_en         : envelope of time series, the same size with s_ori

[N,~]=size(s_ori);
s_en = zeros(size(s_ori));
for i=1:N
    s_en(i,:) = abs(hilbert(s_ori(i,:)));
end
