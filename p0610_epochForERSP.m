% 08/05/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all


%% Prepare the common data.
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


%% Perform 40Hz-focused epoch selection.

% Register all the precomputed ERSP data. 
allMats     = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMean.mat');
stackData   = zeros(128, 100, 5000, length(allMats));
allMatNames = cellfun(@(x) x(1:4), {allMats.name}, 'UniformOutput',false)';

% Loop across the set files.
allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    % 40-Hz band-pass filter.
    EEG = pop_firws(EEG, 'fcutoff', [37.5 42.5], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', 660, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Epoch the data.
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [0.5  3], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Nakayoshi special.
    data2D = permute(EEG.data, [3 1 2]);
    data2D = double(data2D(:,:));
    [V,D] = eigs(data2D*data2D', size(data2D, 1));
    firstEigenVector = V(:,1);
    numPositiveEpochs = sum(firstEigenVector>0);
    if numPositiveEpochs>=50
        [~, sortingIdx] = sort(firstEigenVector, 'ascend');
    else
        [~, sortingIdx] = sort(firstEigenVector, 'descend');
    end
    exclusionIdx = sortingIdx(1:round(length(sortingIdx)*0.2));

    % Extract good epochs.
    disp(sprintf('%d/%d', setIdx, length(allSets)))
    currentMatIdx = find(strcmp(allMatNames, dataName));
    currentMat    = allMats(currentMatIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    epochRejectedData = elecErspMean;
    epochRejectedData(:,exclusionIdx,:) = [];
    stackData(:,:,:,matIdx) = elecErspMean;








    % Save the data.
    pop_saveset(EEG, 'filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0600_epochForERP_trialRej');
end















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
imagesc(timeBins, 1:100, squeeze(erspFXS(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at Cz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(erspTDC(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at Cz', length(tdcIdx)))
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

print('/srv/Makoto/ASSR/p0300_ersp_itc/czMeanErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0300_ersp_itc/czMeanErsp', '-r200', '-dsvg')


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
imagesc(timeBins, 1:100, squeeze(erspFXS(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at Cz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(erspTDC(55,:,:)), [-1.5 1.5]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at Cz', length(tdcIdx)))
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

print('/srv/Makoto/ASSR/p0300_ersp_itc/czMedianErsp', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0300_ersp_itc/czMedianErsp', '-r200', '-dsvg')


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
imagesc(timeBins, 1:100, squeeze(itcFXS(55,:,:)), [0 0.35]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at Cz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(itcTDC(55,:,:)), [0 0.35]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at Cz', length(tdcIdx)))
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

print('/srv/Makoto/ASSR/p0300_ersp_itc/czMeanItc', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0300_ersp_itc/czMeanItc', '-r200', '-dsvg')


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
imagesc(timeBins, 1:100, squeeze(itcFXS(55,:,:)), [0 0.55]); colormap jet; axis xy; title(sprintf('FXS (n=%d) at Cz', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

subplot(1,2,2)
imagesc(timeBins, 1:100, squeeze(itcTDC(55,:,:)), [0 0.55]); colormap jet; axis xy; title(sprintf('TDC (n=%d) at Cz', length(tdcIdx)))
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

print('/srv/Makoto/ASSR/p0300_ersp_itc/czMedianItc', '-r200', '-djpeg95')
print('/srv/Makoto/ASSR/p0300_ersp_itc/czMedianItc', '-r200', '-dsvg')
