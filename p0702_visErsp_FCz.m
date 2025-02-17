% 12/03/2024 Makoto. Used again.
% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools/groupSIFT-master')
addpath('/srv/Makoto/Tools/frank-pk-DataViz-3.2.3.0/daboxplot')

% Preselect datasets recorded from 13+ yrs old subjects.
allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
cleaningMatrix = zeros(length(allSets),5); % numCh, meandB, stddB, mindB, maxdB.
ageList        = cell(length(allSets),1);
sexGroupList = cell(length(allSets),2);
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath, 'loadmode', 'info');

    cleaningMatrix(setIdx,:) = [length(EEG.etc.rejectedChannelIdx) mean(EEG.etc.asrPowerReductionDb) std(EEG.etc.asrPowerReductionDb) min(EEG.etc.asrPowerReductionDb) max(EEG.etc.asrPowerReductionDb)];

    % EEG.etc.Sex
    % EEG.etc.Dx
    sexGroupList{setIdx,1} = EEG.etc.Sex{1,1};
    sexGroupList{setIdx,2} = EEG.etc.Dx{1,1};
    ageList{setIdx,1}      = EEG.etc.AgeAtVisit{1,1};
end

subjNames  = cellfun(@(x) x(1:4), {allSets.name},    'UniformOutput', false);
groupNames = cellfun(@(x) x(),    sexGroupList(:,2), 'UniformOutput', false);

ageVec = str2num(cell2mat(cellfun(@(x) x(1:2), ageList, 'UniformOutput', false)));
    %{
    figure
    bar(sort(ageVec))
    %}

% Load dummy data.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/p0000_raw/Final_SSCT_Adults_Adolescents_editedLAD_11152024_modified.xlsx');
groupIdx = demographicData.Var6;

% Sort the subjIDs.
[subjID_sorted, sortingIdx] = sort(demographicData.Var1);
groupIdx = groupIdx(sortingIdx);

hcsIdx = find(strcmp(groupIdx, 'TDC'));
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
freqIdx      = [freqIdx1 freqIdx2 freqIdx3 freqIdx4 freqIdx5 freqIdx6 freqIdx7 freqIdx8];
freqLabels   = [1 2 4 8 13 20 40 80];

thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
gammaIdx = freqIdx7;

timeBins    = -1000:3999;
baselineIdx = find(timeBins>=-1000 & timeBins<=0);
onsetIdx    = find(timeBins>0 & timeBins<=500);
assrIdx     = find(timeBins>0 & timeBins<=3000);
offsetIdx   = find(timeBins>3000 & timeBins<=3500);
afterZeroIdx = find(timeBins>0);

%% Load Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');

stackData = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecErspMedian;
end
stackData      = 10*log10(stackData);
stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,baselineIdx,:),3));
erspFXS        = mean(stackData_base(:,:,:,fxsIdx),4);
erspHCS        = mean(stackData_base(:,:,:,hcsIdx),4);

% Load p-value tensors.
load('/srv/Makoto/ASSR/p0700_visErssp_Ersp/ersp_pTensor_FXS.mat')
ersp_pTensor_FXS = pTensor_rTale;
load('/srv/Makoto/ASSR/p0700_visErssp_Ersp/ersp_pTensor_HCS.mat')
ersp_pTensor_HCS = pTensor_rTale;
clear pTensor_rTale


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot ERSP at FCz. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the time-freq domain ROI for FXS.
elecIdx = 6; % 6, FCz; 55, Cz.
[~, index40Hz] = min(abs(wtFreqBins-40));
FCzERSP_FXS = squeeze(erspFXS(elecIdx,:,:));
FCzERSP_HCS = squeeze(erspHCS(elecIdx,:,:));

% Plot.
figure('position', [50 50 1800 600])
subplot(1,3,1)
imagesc(timeBins, 1:100, FCzERSP_FXS, [-2 2])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

subplot(1,3,2)
imagesc(timeBins, 1:100, FCzERSP_HCS, [-2 2])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

subplot(1,3,3)

erspFXS_group = permute(stackData_base(elecIdx,:,:,fxsIdx), [4 2 3 1]);
erspHCS_group = permute(stackData_base(elecIdx,:,:,hcsIdx), [4 2 3 1]);
[H,P,CI,STATS] = ttest2(erspFXS_group, erspHCS_group);
significanceCutoff = tinv([0.025 0.975], median(STATS.df(:)));
plottingMaps = squeeze(STATS.tstat);
% plottingMaps(plottingMaps>significanceCutoff(1) & plottingMaps<significanceCutoff(2)) = 0;
% imagesc(timeBins, 1:100, plottingMaps, [-4.5 4.5])

% Cluster-level correction.
[mask, tScore, pValue, surroMassOfCluster] = clusterLevelPermutationTest(permute(erspFXS_group, [2 3 1]), permute(erspHCS_group, [2 3 1]), 0, 0.05, 200);

massOfClusterVec = zeros(max(mask(:)),1);
for clsIdx = 1:max(mask(:))
    massOfClusterVec(clsIdx) = sum(tScore(mask==clsIdx)); % Output is 1D.
end
threshold = [prctile(surroMassOfCluster(:,1),5) prctile(surroMassOfCluster(:,2),95)];
survivedClusterIdx = find(massOfClusterVec<threshold(1) | massOfClusterVec>threshold(2));
finalMask = logical(zeros(size(mask)));
for clsIdx = 1:length(survivedClusterIdx)
    currentClsIdx = survivedClusterIdx(clsIdx);
    survivedMask = mask==currentClsIdx;
    finalMask = finalMask+survivedMask;
end

grandMeanFXS = squeeze(mean(erspFXS_group,1));
grandMeanHCS = squeeze(mean(erspHCS_group,1));
axesHandle = imagesc(timeBins, 1:100, plottingMaps.*finalMask, [-4.5 4.5]);
axis xy;
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels);
customColormap = jet(64);
customColormap(33,:) = [0.7 0.7 0.7];
colormap(axesHandle.Parent, customColormap)
axis xy
axis square
title('FXS-HCS')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 't')

sgtitle('Grand-mean ERSP at FCz')

print('/srv/Makoto/ASSR/p0702_visErsp_FCz/ERSP', '-djpeg', '-r200')
print('/srv/Makoto/ASSR/p0702_visErsp_FCz/ERSP', '-dsvg', '-r200')