% 12/03/2024 Makoto. Used again.
% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

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
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');

stackData = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMedian;
end
% stackData      = 10*log10(stackData);
% stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,baselineIdx,:),3));
itcFXS         = mean(stackData(:,:,:,fxsIdx),4);
itcHCS         = mean(stackData(:,:,:,hcsIdx),4);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot ERSP at FCz. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the time-freq domain ROI for FXS.
elecIdx = 6; % 6, FCz; 55, Cz.
[~, index40Hz] = min(abs(wtFreqBins-40));
FCzITC_FXS = squeeze(itcFXS(elecIdx,:,:));
FCzITC_HCS = squeeze(itcHCS(elecIdx,:,:));

% Plot.
figure('position', [50 50 1800 600])
subplot(1,3,1)
imagesc(timeBins, 1:100, FCzITC_FXS, [0 0.7])
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
set(get(cbarHandle, 'title'), 'string', 'ITC')

subplot(1,3,2)
imagesc(timeBins, 1:100, FCzITC_HCS, [0 0.7])
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
set(get(cbarHandle, 'title'), 'string', 'ITC')

sgtitle('Grand-mean ITC at FCz')

print('/srv/Makoto/ASSR/p0802_visItc_FCz/Itc', '-djpeg', '-r200')
print('/srv/Makoto/ASSR/p0802_visItc_FCz/Itc', '-dsvg', '-r200')