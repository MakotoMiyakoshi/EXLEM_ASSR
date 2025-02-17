% 12/03/2024 Makoto. Used again.
% 10/31/2024 Makoto. Calculate delta-gamma power-power coupling.
clear
close all
clc

% Load delta-gamma coupling.
load('/srv/Makoto/ASSR/p0900_deltaGammaCoupling/rMatrix.mat')

% Load dummy data.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');

% Separate groups.
fxsData = rMatrix(:, fxsIdx);
hcsData = rMatrix(:, hcsIdx);

% Stats.
[H,P,CI,STATS] = ttest2(fxsData', hcsData');
nonsignificantIdx = find(P>=0.05);
tstats = STATS.tstat;
tstats(nonsignificantIdx) = 0;

% Plot.
figure
subplot(1,3,1)
topoplot(mean(fxsData,2), EEG.chanlocs, 'maplimits', [-0.05 0.05])
title('FXS')
colorbar

subplot(1,3,2)
topoplot(mean(hcsData,2), EEG.chanlocs, 'maplimits', [-0.05 0.05])
title('HCS')
colorbar

subplot(1,3,3)
topoplot(tstats, EEG.chanlocs, 'maplimits', [-3.2 3.2])
title('FXS-HCS where p<0.05')
colorbar

exportgraphics(gcf, '/srv/Makoto/ASSR/p0910_deltaGammaCouplingTopos/ERSP_withDelta.pdf', 'ContentType', 'vector')