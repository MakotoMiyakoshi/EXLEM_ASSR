% 07/12/2024 Makoto. Created.


%% Original.
clear
clc
close all

% Compute grand-mean and -median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*elecErspMean.mat');
meanERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    meanERSP(  :,:,:,matIdx) = elecErspMean;
end
meanERSP_log = 10*log10(meanERSP);
grandMeanERSP = mean(meanERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP', 'grandMeanERSP', '-nocompression')


clear
clc
close all

% Compute grand-mean and -median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*elecErspMedian.mat');
medianERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    medianERSP(  :,:,:,matIdx) = elecErspMedian;
end
medianERSP_log = 10*log10(medianERSP);
grandMedianERSP = mean(medianERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP', 'grandMedianERSP', '-nocompression')
disp('Done.')

%% Brain-class.
clear
clc
close all

allMats = dir('/srv/Makoto/ASSR/p0420_epochBrain/*elecErspMean.mat');
meanERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    meanERSP(  :,:,:,matIdx) = elecErspMean;
end
meanERSP_log = 10*log10(meanERSP);
grandMeanERSP = mean(meanERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP_brain', 'grandMeanERSP', '-nocompression')


clear
clc
close all

% Compute grand-mean and -median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0420_epochBrain/*elecErspMedian.mat');
medianERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    medianERSP(  :,:,:,matIdx) = elecErspMedian;
end
medianERSP_log = 10*log10(medianERSP);
grandMedianERSP = mean(medianERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP_brain', 'grandMedianERSP', '-nocompression')
disp('Done.')

%% Non-Brain-class.
clear
clc
close all

allMats = dir('/srv/Makoto/ASSR/p0430_epochNonBrain/*elecErspMean.mat');
meanERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    meanERSP(  :,:,:,matIdx) = elecErspMean;
end
meanERSP_log = 10*log10(meanERSP);
grandMeanERSP = mean(meanERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP_nonBrain', 'grandMeanERSP', '-nocompression')


clear
clc
close all

% Compute grand-mean and -median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0430_epochNonBrain/*elecErspMedian.mat');
medianERSP   = single(zeros(128, 100, 5000, length(allMats)));
for matIdx = 1:length(allMats)

    loadName = allMats(matIdx).name;
    loadPath = allMats(matIdx).folder;
    load([loadPath filesep loadName]);
    medianERSP(  :,:,:,matIdx) = elecErspMedian;
end
medianERSP_log = 10*log10(medianERSP);
grandMedianERSP = mean(medianERSP_log,4);
save('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP_nonBrain', 'grandMedianERSP', '-nocompression')
disp('Done.')


%% Visualize Mean ERSP.
clear
clc
close all

originalData      = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP.mat');
brainClassData    = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP_brain.mat');
nonBrainClassData = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMeanERSP_nonBrain.mat');

brainPlusNonBrain = brainClassData.grandMeanERSP + nonBrainClassData.grandMeanERSP;

residual = originalData.grandMeanERSP - brainPlusNonBrain;


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


figure('position', [100         600        2350         435])

subplot(1,4,1)
plotData = squeeze(originalData.grandMeanERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Original at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,2)
plotData = squeeze(brainClassData.grandMeanERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Brain-class at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,3)
plotData = squeeze(nonBrainClassData.grandMeanERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('NonBrain-class at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,4)
plotData = squeeze(residual(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Original - (Brain + Non-Brain)')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0470_validation_ERSP/meanERSP', '-r200', '-djpeg95')


%% Visualize median ERSP.
clear
clc
close all

originalData      = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP.mat');
brainClassData    = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP_brain.mat');
nonBrainClassData = load('/srv/Makoto/ASSR/p0470_validation_ERSP/grandMedianERSP_nonBrain.mat');

brainPlusNonBrain = brainClassData.grandMedianERSP + nonBrainClassData.grandMedianERSP;

residual = originalData.grandMedianERSP - brainPlusNonBrain;


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


figure('position', [100         600        2350         435])

subplot(1,4,1)
plotData = squeeze(originalData.grandMedianERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Original at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,2)
plotData = squeeze(brainClassData.grandMedianERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Brain-class at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,3)
plotData = squeeze(nonBrainClassData.grandMedianERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('NonBrain-class at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,4,4)
plotData = squeeze(residual(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Original - (Brain + Non-Brain)')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0470_validation_ERSP/medianERSP', '-r200', '-djpeg95')


%%
figure('position', [100         600        1200         435])

subplot(1,2,1)
plotData = squeeze(originalData.grandMedianERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Original at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
plotData = squeeze(brainClassData.grandMedianERSP(55,:,:)+nonBrainClassData.grandMedianERSP(55,:,:));
plotData = bsxfun(@minus, plotData, mean(plotData(:,1:1000),2));
imagesc(-1000:3999, 1:100, plotData, [-1 1]);
colormap jet
axis xy
title('Brain + Non-brain at Cz')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

originalPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')
set(gca,'position', originalPosition)

print('/srv/Makoto/ASSR/p0470_validation_ERSP/splitAndCombineTest', '-r200', '-djpeg95')