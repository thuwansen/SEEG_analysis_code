function s_norm = z_signal_norm(s_ori)
% This function normalize the time series.
%INPUTS
% - s_ori        : time series to be normalized, n x samples
%OUTPUTS
% - s_norm       : normalized time series, the same size with s_ori
[~,s_length]=size(s_ori);
s_mean = mean(s_ori,2);
s_std = std(s_ori,1,2);
s_norm = (s_ori-repmat(s_mean,[1,s_length]))./repmat(s_std,[1,s_length]);