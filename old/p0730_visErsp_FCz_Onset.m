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
load('/srv/Makoto/ASSR/p0700_visErsp/ersp_pTensor_FXS.mat')
ersp_pTensor_FXS = pTensor_rTale;
load('/srv/Makoto/ASSR/p0700_visErsp/ersp_pTensor_HCS.mat')
ersp_pTensor_HCS = pTensor_rTale;
clear pTensor_rTale


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FCz ERSP: Gamma ASSR. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the time-freq domain ROI for FXS.
elecIdx = 6; % 6, FCz; 55, Cz.
[~, index40Hz] = min(abs(wtFreqBins-40));
FCzERSP_FXS = squeeze(erspFXS(elecIdx,:,:));
FCzERSP_HCS = squeeze(erspHCS(elecIdx,:,:));
erspROI_FXS = squeeze(mean(stackData_base(elecIdx,index40Hz,assrIdx,fxsIdx),3));
erspROI_HCS = squeeze(mean(stackData_base(elecIdx,index40Hz,assrIdx,hcsIdx),3));

[H,P,CI,STATS] = ttest2(erspROI_FXS, erspROI_HCS);

% % LME.
% allAbove13Idx     = find(ageVec>12);
% subjNames_13plus  = subjNames( allAbove13Idx)';
% groupNames_13plus = groupNames(allAbove13Idx)';
% ageVec_13plus     = ageVec(    allAbove13Idx)';
% fxsIdx_preselected = find(contains(groupNames_13plus, 'FXS'));
% hcsIdx_preselected = find(contains(groupNames_13plus, 'TDC'));
% fxsAge_preselected = ageVec_13plus(fxsIdx_preselected);
% tdcAge_preselected = ageVec_13plus(hcsIdx_preselected);
% 
% % Combine data to a table.
% SubjectIdx    = [1:length(erspROI_FXS)+length(erspROI_HCS)]';
% Group         = [repmat({'FXS'}, [length(erspROI_FXS), 1]); repmat({'HCS'}, [length(erspROI_HCS), 1])];
% Age           = [fxsAge_preselected'; tdcAge_preselected'];
% meanERSP      = [erspROI_FXS; erspROI_HCS];
% gammaTable = table(SubjectIdx, Group, Age, meanERSP, 'variableNames', {'SubjectIdx' 'Group', 'Age', 'meanERSP'});
% gammaTable.SubjectIdx = categorical(gammaTable.SubjectIdx);
% gammaTable.Group      = categorical(gammaTable.Group);
% 
% % Start LME.
% lmeModel1 = fitlme(gammaTable, 'meanERSP ~ 1 + Group + Age + (1|SubjectIdx)');
% lmeModel2 = fitlme(gammaTable, 'meanERSP ~ 1 + Group * Age + (1|SubjectIdx)');
% % lmeModel3 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + Group * Age + (1|SubjectIdx)'); % meaningless.
% compare(lmeModel1, lmeModel2)
% lmeModel1.Rsquared % R^2 0.7525
% sqrt(lmeModel1.Rsquared.Adjusted) % r=0.8675.


%% Visualize the 1x3 plot.
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
h = daboxplot({erspROI_FXS, erspROI_HCS},'linkline', 1, 'xtlabels', {'FXS' 'TDC'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0 0.4470 0.7410];
ylabel('Power (dB)')
xlim([0.5 2.5])
ylim([0 3])
title(sprintf('40Hz power, t(%d)=%.2f, p=%.4f', STATS.df, STATS.tstat, P))
%title(sprintf('40Hz power, R^2=%.3f, t(%d)=%.2f, p=%.4f', lmeModel1.Rsquared.Adjusted, lmeModel1.Coefficients.DF(2), lmeModel1.Coefficients.tStat(2), lmeModel1.Coefficients.pValue(2)))
axis square

sgtitle('Grand-mean ERSP at FCz')

print('/srv/Makoto/ASSR/p0720_visErsp_FCz_40Hz/ERSP_40Hz', '-djpeg', '-r200')
print('/srv/Makoto/ASSR/p0720_visErsp_FCz_40Hz/ERSP_40Hz', '-dsvg', '-r200')