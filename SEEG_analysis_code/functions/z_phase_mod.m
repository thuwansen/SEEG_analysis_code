function phase_normal = z_phase_mod(phase_diff)
% This function transform the phase difference between -180deg~180deg
%INPUTS
% - phase_diff   : original phase difference, in deg
%OUTPUTS
% - phase_normal : phase difference between -180deg~180deg

phase_normal = phase_diff;
phase_normal(phase_diff>180)=phase_normal(phase_diff>180)-360;
phase_normal(phase_diff<=-180)=phase_normal(phase_diff<=-180)+360;
