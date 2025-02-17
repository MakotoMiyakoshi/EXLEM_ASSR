% 07/10/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools/frank-pk-DataViz-3.2.3.0/daboxplot')

% Load dummy EEG.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0100_upToDipfit/', 'loadmode', 'info');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/code/Subset_SSCT_for_Ernie.xlsx');
groupIdx = demographicData.Dx;
tdcIdx = find(strcmp(groupIdx, 'TDC'));
fxsIdx = find(strcmp(groupIdx, 'FXS'));

% Generate frequency bins.
addpath('/srv/Makoto/Tools/siyisCodeFromRamesh')
freqRange   = [1 100];
numFreqBins = 100;
wtFreqBins = logspace(log10(1), log10(100), numFreqBins);
[~,freqIdx1] = min(abs(wtFreqBins-1));
[~,freqIdx2] = min(abs(wtFreqBins-2));
[~,freqIdx3] = min(abs(wtFreqBins-4));
[~,freqIdx4] = min(abs(wtFreqBins-8));
[~,freqIdx5] = min(abs(wtFreqBins-13));
[~,freqIdx6] = min(abs(wtFreqBins-20));
[~,freqIdx7] = min(abs(wtFreqBins-40));
[~,freqIdx8] = min(abs(wtFreqBins-80));
freqIdx = [freqIdx1 freqIdx2 freqIdx3 freqIdx4 freqIdx5 freqIdx6 freqIdx7 freqIdx8];
freqLabels = [1 2 4 8 13 20 40 80];
timeBins = -1000:3999;


%% Cz Median ERSP.

allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMedian;
end
itcTDC = stackData(:,:,:,tdcIdx);
itcFXS = stackData(:,:,:,fxsIdx);

% Prepare the plotting data.
thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
gammaIdx = freqIdx7;

onsetIdx  = find(timeBins>=0 & timeBins<=500);
assrIdx   = find(timeBins>=0 & timeBins<=3000);
offsetIdx = find(timeBins>=3000 & timeBins<=3500);

fxsThetaOnset = squeeze(mean(mean(itcFXS(:, thetaIdx, onsetIdx, :),2),3));
tdcThetaOnset = squeeze(mean(mean(itcTDC(:, thetaIdx, onsetIdx, :),2),3));

fxs40Hz = squeeze(mean(mean(itcFXS(:, gammaIdx, assrIdx, :),2),3));
tdc40Hz = squeeze(mean(mean(itcTDC(:, gammaIdx, assrIdx, :),2),3));

fxsThetaOffset = squeeze(mean(mean(itcFXS(:, thetaIdx, offsetIdx, :),2),3));
tdcThetaOffset = squeeze(mean(mean(itcTDC(:, thetaIdx, offsetIdx, :),2),3));

%%
figure('position', [800   100   700   800])
subplot(3,3,1)
topoplot(mean(fxsThetaOnset,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('FXS Onset Theta')

subplot(3,3,2)
topoplot(mean(fxs40Hz,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('FXS 40Hz ASSR')

subplot(3,3,3)
topoplot(mean(fxsThetaOffset,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('FXS Offset Theta')
currentPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'ITC[0-1]');
set(gca,'position', currentPosition);

subplot(3,3,4)
topoplot(mean(tdcThetaOnset,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('TDC Onset Theta')

subplot(3,3,5)
topoplot(mean(tdc40Hz,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('TDC 40Hz ASSR')

subplot(3,3,6)
topoplot(mean(tdcThetaOffset,2), EEG.chanlocs, 'maplimits', [-0.6 0.6])
title('TDC Offset Theta')

subplot(3,3,7)
fxsData = fxsThetaOnset(55,:);
tdcData = tdcThetaOnset(55,:);
[H,P,CI,STATS] = ttest2(fxsData, tdcData);
h = daboxplot({fxsData', tdcData'},'linkline', 1, 'xtlabels', {'FXS' 'TDC'}, 'whiskers',0,'outliers',1,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0 0.4470 0.7410];
ylabel('ITC [0-1] at Cz')
xlim([0.5 2.5])
ylim([0 0.7])
title(sprintf('t(%d)=%.1f, p=%.3f', STATS.df, STATS.tstat, P))

subplot(3,3,8)
fxsData = fxs40Hz(55,:);
tdcData = tdc40Hz(55,:);
[H,P,CI,STATS] = ttest2(fxsData, tdcData);
h = daboxplot({fxsData', tdcData'},'linkline', 1, 'xtlabels', {'FXS' 'TDC'}, 'whiskers',0,'outliers',1,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0 0.4470 0.7410];
ylabel('ITC [0-1] at Cz')
xlim([0.5 2.5])
ylim([0 0.7])
title(sprintf('t(%d)=%.1f, p=%.3f', STATS.df, STATS.tstat, P))

subplot(3,3,9)
fxsData = fxsThetaOffset(55,:);
tdcData = tdcThetaOffset(55,:);
[H,P,CI,STATS] = ttest2(fxsData, tdcData);
h = daboxplot({fxsData', tdcData'},'linkline', 1, 'xtlabels', {'FXS' 'TDC'}, 'whiskers',0,'outliers',1,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0 0.4470 0.7410];
ylabel('ITC [0-1] at Cz')
xlim([0.5 2.5])
ylim([0 0.7])
title(sprintf('t(%d)=%.1f, p=%.3f', STATS.df, STATS.tstat, P))

sgtitle('ITC')

exportgraphics(gcf, '/srv/Makoto/ASSR/p0340_2x2TopoComparison_ITC/ITC.pdf', 'ContentType', 'vector')
% print('/srv/Makoto/ASSR/p0340_2x2TopoComparison_ITC/ITC', '-r200', '-djpeg95')
% print('/srv/Makoto/ASSR/p0340_2x2TopoComparison_ITC/ITC', '-r200', '-dsvg')