% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

% Load dummy data.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');

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


%% Cz Mean ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMean.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecErspMean;
end
stackData      = 10*log10(stackData);
stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,1:1000,:),3));
erspTDC        = mean(stackData_base(:,:,:,tdcIdx),4);
erspFXS        = mean(stackData_base(:,:,:,fxsIdx),4);

figure('position', [582 471 1300 400])
subplot(1,2,1)
imagesc(timeBins, 1:100, squeeze(erspFXS(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at FCz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(erspTDC(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at FCz', length(tdcIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'ERSP (dB)')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0301_visErspItc/fczMeanErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/fczMeanErsp', '-r200', '-dsvg')


% Scalp topos.
figure('position', [582 471 1300 400])
subplot(1,2,1)
topoplot(mean(erspFXS(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.8 0.8]);
title('FXS grand-mean 40-Hz ERSP')
subplot(1,2,2)
topoplot(mean(erspTDC(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.8 0.8]);
title('TDC grand-mean 40-Hz ERSP')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

print('/srv/Makoto/ASSR/p0301_visErspItc/topoMeanErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/topoMeanErsp', '-r200', '-dsvg')


% Individual FXS scalp topos.
subjNames = cellfun(@(x) x(1:4), {allMats.name}, 'UniformOutput', false);
figure
for withinFxsIdx = 1:length(fxsIdx)
    plottingData = mean(stackData_base(:,freqIdx7,1000:4000,fxsIdx(withinFxsIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
        subjIdx = fxsIdx(withinFxsIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (FXS) 40-Hz mean ERSP', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMeanErsp_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end

% Individual TDC scalp topos.
for withinTdcIdx = 1:length(tdcIdx)
    plottingData = mean(stackData_base(:,freqIdx7,1000:4000,tdcIdx(withinTdcIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
        subjIdx = tdcIdx(withinTdcIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (TDC) 40-Hz mean ERSP', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMeanErsp_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end


%% Cz Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecErspMedian;
end
stackData      = 10*log10(stackData);
stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,1:1000,:),3));
erspTDC        = mean(stackData_base(:,:,:,tdcIdx),4);
erspFXS        = mean(stackData_base(:,:,:,fxsIdx),4);

figure('position', [582 471 1300 400])
subplot(1,2,1)
imagesc(timeBins, 1:100, squeeze(erspFXS(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at FCz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(erspTDC(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at FCz', length(tdcIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'ERSP (dB)')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0301_visErspItc/fczMedianErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/fczMedianErsp', '-r200', '-dsvg')


% Scalp topos.
figure('position', [582 471 1300 400])
subplot(1,2,1)
topoplot(mean(erspFXS(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-1.0 1.0]);
title('FXS grand-mean 40-Hz ERSP')
subplot(1,2,2)
topoplot(mean(erspTDC(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-1.0 1.0]);
title('TDC grand-mean 40-Hz ERSP')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

print('/srv/Makoto/ASSR/p0301_visErspItc/topoMedianErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/topoMedianErsp', '-r200', '-dsvg')


% Individual FXS scalp topos.
subjNames = cellfun(@(x) x(1:4), {allMats.name}, 'UniformOutput', false);
figure
for withinFxsIdx = 1:length(fxsIdx)
    plottingData = mean(stackData_base(:,freqIdx7,1000:4000,fxsIdx(withinFxsIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
        subjIdx = fxsIdx(withinFxsIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (FXS) 40-Hz median ERSP', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMedianErsp_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end

% Individual TDC scalp topos.
for withinTdcIdx = 1:length(tdcIdx)
    plottingData = mean(stackData_base(:,freqIdx7,1000:4000,tdcIdx(withinTdcIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
        subjIdx = tdcIdx(withinTdcIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (TDC) 40-Hz median ERSP', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMedianErsp_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end


%% Cz Mean ITC.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMean.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMean;
end
itcTDC = mean(stackData(:,:,:,tdcIdx),4);
itcFXS = mean(stackData(:,:,:,fxsIdx),4);

figure('position', [582 471 1300 400])
subplot(1,2,1)
imagesc(timeBins, 1:100, squeeze(itcFXS(55,:,:)), [0 0.35]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at FCz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(itcTDC(55,:,:)), [0 0.35]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at FCz', length(tdcIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'ITC [0-1]')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0301_visErspItc/fczMeanItc', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/fczMeanItc', '-r200', '-dsvg')


% Scalp topos.
figure('position', [582 471 1300 400])
subplot(1,2,1)
topoplot(mean(itcFXS(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.35 0.35]);
title('FXS grand-mean 40-Hz ITC')
subplot(1,2,2)
topoplot(mean(itcTDC(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.35 0.35]);
title('TDC grand-mean 40-Hz ITC')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '[0-1]')

print('/srv/Makoto/ASSR/p0301_visErspItc/topoMeanItc', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/topoMeanItc', '-r200', '-dsvg')


% Individual FXS scalp topos.
subjNames = cellfun(@(x) x(1:4), {allMats.name}, 'UniformOutput', false);
figure
for withinFxsIdx = 1:length(fxsIdx)
    plottingData = mean(stackData(:,freqIdx7,1000:4000,fxsIdx(withinFxsIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '[0-1]')
        subjIdx = fxsIdx(withinFxsIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (FXS) 40-Hz mean ITC', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMeanItc_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end

% Individual TDC scalp topos.
for withinTdcIdx = 1:length(tdcIdx)
    plottingData = mean(stackData(:,freqIdx7,1000:4000,tdcIdx(withinTdcIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '[0-1]')
        subjIdx = tdcIdx(withinTdcIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (TDC) 40-Hz mean ITC', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMeanItc_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end


%% Cz Median ITC.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMedian;
end
itcTDC = mean(stackData(:,:,:,tdcIdx),4);
itcFXS = mean(stackData(:,:,:,fxsIdx),4);

figure('position', [582 471 1300 400])
subplot(1,2,1)
imagesc(timeBins, 1:100, squeeze(itcFXS(55,:,:)), [0 0.55]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at FCz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(itcTDC(55,:,:)), [0 0.55]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at FCz', length(tdcIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', set(get(cbarHandle, 'title'), 'string', 'ITC [0-1]'))
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0301_visErspItc/fczMedianErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/fczMedianErsp', '-r200', '-dsvg')


% Scalp topos.
figure('position', [582 471 1300 400])
subplot(1,2,1)
topoplot(mean(itcFXS(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.55 0.55]);
title('FXS grand-mean 40-Hz ITC')
subplot(1,2,2)
topoplot(mean(itcTDC(:,freqIdx7,1000:4000),3), EEG.chanlocs, 'maplimits', [-0.55 0.55]);
title('TDC grand-mean 40-Hz ITC')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '[0-1]')

print('/srv/Makoto/ASSR/p0301_visErspItc/topoMedianItc', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0301_visErspItc/topoMedianItc', '-r200', '-dsvg')


% Individual FXS scalp topos.
subjNames = cellfun(@(x) x(1:4), {allMats.name}, 'UniformOutput', false);
figure
for withinFxsIdx = 1:length(fxsIdx)
    plottingData = mean(stackData(:,freqIdx7,1000:4000,fxsIdx(withinFxsIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '[0-1]')
        subjIdx = fxsIdx(withinFxsIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (FXS) 40-Hz median ITC', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMedianItc_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end

% Individual TDC scalp topos.
for withinTdcIdx = 1:length(tdcIdx)
    plottingData = mean(stackData(:,freqIdx7,1000:4000,tdcIdx(withinTdcIdx)),3);
    colorRange = [-max(abs(plottingData)) max(abs(plottingData))];
    topoplot(plottingData, EEG.chanlocs, 'maplimits', colorRange);
    currentHandle = gca;
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '[0-1]')
        subjIdx = tdcIdx(withinTdcIdx);
    subjName = subjNames{subjIdx};
    title(sprintf('%s (TDC) 40-Hz median ITC', subjName))
    print(sprintf('/srv/Makoto/ASSR/p0301_visErspItc/individualTopos/topoMedianItc_%s', subjName), '-r200', '-djpeg95')
    cla(currentHandle)
end


%% Visualize 2287's ITC at FCz.
targetIdx = find(contains(subjNames, '2287'));
itcTensor = stackData(:,:,:,targetIdx);

figure
imagesc(squeeze(itcTensor(55,:,:)), [0 0.6])
axis xy
colormap('jet')