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

% tooYoungIdx = find(ageVec<6);
% subjNames( tooYoungIdx)'
% groupNames(tooYoungIdx)'
% 
% earlyAdolescenceIdx = find(ageVec>=6 & ageVec<=12);
% subjNames( earlyAdolescenceIdx)'
% groupNames(earlyAdolescenceIdx)'

allAbove13Idx = find(ageVec>12);
subjNames_13plus  = subjNames( allAbove13Idx)';
groupNames_13plus = groupNames(allAbove13Idx)';

% Load dummy data.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/code/Subset_SSCT_for_Ernie.xlsx');
groupIdx = demographicData.Dx;

% Apply preselection for 13+ yrs olds.
groupIdx = groupIdx(allAbove13Idx);

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

% Apply preselection for 13+ yrs olds.
allMats = allMats(allAbove13Idx);

stackData = zeros(128, 100, 5000, length(allMats));
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
load('/srv/Makoto/ASSR/p0700_visErspItc_13yrsPlus_stats/ersp_pTensor_FXS.mat')
ersp_pTensor_FXS = pTensor_rTale;
load('/srv/Makoto/ASSR/p0700_visErspItc_13yrsPlus_stats/ersp_pTensor_HCS.mat')
ersp_pTensor_HCS = pTensor_rTale;
clear pTensor_rTale


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FCz ERSP: Gamma ASSR. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the time-freq domain ROI for FXS.
elecIdx = 11; % 6, FCz; 55, Cz; 11, Fz.
CzERSP_FXS = squeeze(erspFXS(elecIdx,:,:));
CzERSP_HCS = squeeze(erspHCS(elecIdx,:,:));

% Define functional ROI for FXS ERSP at '40Hz ROI' (i.e. 'Lemniscal ROI').
CzERSP_sigMask_FXS = squeeze(ersp_pTensor_FXS(elecIdx,:,:))<0.05;
    %{
    figure
    imagesc(CzERSP_sigMask_FXS)
    %}
CzERSP_sigMask_FXS(freqIdx2:end, :) = 0;

% Define functional ROI for HCS ERSP at '40Hz ROI' (i.e. 'Lemniscal ROI').
CzERSP_sigMask_HCS = squeeze(ersp_pTensor_HCS(elecIdx,:,:))<0.05;
    %{
    figure
    imagesc(CzERSP_sigMask_HCS)
    %}
CzERSP_sigMask_HCS(freqIdx2:end, :) = 0;

% Extract ERSP values.
CzERSP_sigMask_FXS_3D = permute(CzERSP_sigMask_FXS, [3 1 2]);
erspROI_FXS = squeeze(bsxfun(@times, stackData_base(elecIdx,:,:,fxsIdx), CzERSP_sigMask_FXS_3D));
erspROI_FXS = squeeze(sum(sum(erspROI_FXS,1),2)/sum(sum(CzERSP_sigMask_FXS)));

CzERSP_sigMask_HCS_3D = permute(CzERSP_sigMask_HCS, [3 1 2]);
erspROI_HCS = squeeze(bsxfun(@times, stackData_base(elecIdx,:,:,hcsIdx), CzERSP_sigMask_HCS_3D));
erspROI_HCS = squeeze(sum(sum(erspROI_HCS,1),2)/sum(sum(CzERSP_sigMask_HCS)));


%%%%%%%%%%%%
%%% LME. %%%
%%%%%%%%%%%%
allAbove13Idx     = find(ageVec>12);
subjNames_13plus  = subjNames( allAbove13Idx)';
groupNames_13plus = groupNames(allAbove13Idx)';
ageVec_13plus     = ageVec(    allAbove13Idx)';

fxsIdx_preselected = find(contains(groupNames_13plus, 'FXS'));
hcsIdx_preselected = find(contains(groupNames_13plus, 'TDC'));

fxsAge_preselected = ageVec_13plus(fxsIdx_preselected);
tdcAge_preselected = ageVec_13plus(hcsIdx_preselected);

% Combine data to a table.
SubjectIdx    = [1:length(erspROI_FXS)+length(erspROI_HCS)]';
Group         = [repmat({'FXS'}, [length(erspROI_FXS), 1]); repmat({'HCS'}, [length(erspROI_HCS), 1])];
Age           = [fxsAge_preselected'; tdcAge_preselected'];
meanERSP      = [erspROI_FXS; erspROI_HCS];

dataTable = table(SubjectIdx, Group, Age, meanERSP, 'variableNames', {'SubjectIdx' 'Group', 'Age', 'meanERSP'});
dataTable.SubjectIdx = categorical(dataTable.SubjectIdx);
dataTable.Group      = categorical(dataTable.Group);

% Start LME.
lmeModel1 = fitlme(dataTable, 'meanERSP ~ 1 + Group + Age + (1|SubjectIdx)');
lmeModel2 = fitlme(dataTable, 'meanERSP ~ 1 + Group * Age + (1|SubjectIdx)');
% lmeModel3 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + Group * Age + (1|SubjectIdx)'); % meaningless.
compare(lmeModel1, lmeModel2)
lmeModel1.Rsquared % R^2 0.7525
sqrt(lmeModel1.Rsquared.Adjusted) % r=0.8675.


%% Visualize the 2x2 plot.
figure('position', [50 50 1800 600])

subplot(1,3,1)
imagesc(timeBins, 1:100, CzERSP_FXS, [-2 2])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS, Cz grand mean (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

subplot(1,3,2)
imagesc(timeBins, 1:100, CzERSP_HCS, [-2 2])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS, Cz grand mean (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

subplot(1,3,3)
h = daboxplot({erspROI_FXS, erspROI_HCS},'linkline', 1, 'xtlabels', {'FXS' 'TDC'}, 'whiskers',0,'outliers',1,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0 0.4470 0.7410];
ylabel('Theta power at Cz')
xlim([0.5 2.5])
ylim([0 3])
title(sprintf('Continuous delta, R^2=%.3f, t(%d)=%.2f, p=%.4f', lmeModel1.Rsquared.Adjusted, lmeModel1.Coefficients.DF(2), lmeModel1.Coefficients.tStat(2), lmeModel1.Coefficients.pValue(2)))

sgtitle('Grand-mean ERSP at Fz')

print('/srv/Makoto/ASSR/p0750_visErspItc_Fz_Delta/ERSP_delta', '-djpeg', '-r200')
print('/srv/Makoto/ASSR/p0750_visErspItc_Fz_Delta/ERSP_delta', '-dsvg',  '-r200')