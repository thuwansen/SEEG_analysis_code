function z_plot_radar(phase_diff,savepath,savename,colorbar)
% This function plot radar map for the phase difference(in deg).
%INPUTS
% - phase_diff   : phase difference in deg format
% - savepaht     : path for the figure
% - savename     : name of the figure to be saved 
% - color bar    : color of the figure


sDist = z_phase_mod(rad2deg(phase_diff)); % generate sample, convert to deg
nBins = 36; % number of bins, makes bin size of 10 deg
% colorBar=[0/255.0,254/255.0,254/255.0];
obj1 = CircHist(sDist, nBins, 'colorBar',colorbar);
obj1.polarAxs.ThetaAxis.MinorTickValues = []; % remove dotted tick-lines
obj1.drawScale; % update scale bar
title('');
rl = rlim; % get current limits
rticks(0:rl(2)/2:rl(2)); 
obj1.scaleBarSide = 'right';
delete([obj1.scaleBar]); 
set(gca,'FontName','Arial','FontSize',25,'LineWidth',2);
saveas(gcf,[savepath,'/',savename]);
close all;

