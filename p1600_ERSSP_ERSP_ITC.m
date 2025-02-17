% 01/08/2025 Makoto. Used.
% 12/10/2024 Makoto. Integrating p1021 and p1100.
% 12/04/2024 Makoto. Modified.
% 12/02/2024 Makoto. Used again.
% 11/29/2024 Makoto. Used again.
% 10/25/2024 Makoto. Used.
% 10/04/2024 Makoto. Modified.
% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools')

% Preselect datasets recorded from 13+ yrs old subjects.
allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
cleaningMatrix = zeros(length(allSets),5); % numCh, meandB, stddB, mindB, maxdB.
sexGroupList = cell(length(allSets),2);
for setIdx = 1:length(allSets)
    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath, 'loadmode', 'info');
    cleaningMatrix(setIdx,:) = [length(EEG.etc.rejectedChannelIdx) mean(EEG.etc.asrPowerReductionDb) std(EEG.etc.asrPowerReductionDb) min(EEG.etc.asrPowerReductionDb) max(EEG.etc.asrPowerReductionDb)];
    sexGroupList{setIdx,1} = EEG.etc.Sex{1,1};
    sexGroupList{setIdx,2} = EEG.etc.Dx{1,1};
end

subjNames  = cellfun(@(x) x(1:4), {allSets.name},    'UniformOutput', false);
groupNames = cellfun(@(x) x(),    sexGroupList(:,2), 'UniformOutput', false);
    %{
    figure
    bar(sort(ageVec))
    %}


%%
% Load the latest demographic data.
demoTable = readtable('/srv/Makoto/ASSR/p0000_raw/Final_List_SSCT_Adults_Adolescents_01062025.xlsx');
eegIdVec = demoTable.eegid;

% Exclusion list by Craig.
exclusionList = [44 532 554 879 2537];
[C, preselectingIdx] = setdiff(eegIdVec, exclusionList);
demoTable = demoTable(preselectingIdx, :);
eegIdVec2 = demoTable.eegid;

% Obtain valid subject IDs.
validEegId   = demoTable.eegid;
subjNamesVec = cellfun(@str2num, subjNames');
[C, preselectingIdx] = intersect(subjNamesVec, validEegId); % Choosing 62 out of 68 is correct.

% Apply the trimming.
subjNames  = subjNames(preselectingIdx);
groupNames = groupNames(preselectingIdx);
sexGroupList = sexGroupList(preselectingIdx,:);

% Validate.
if sum(subjNamesVec(preselectingIdx)-eegIdVec2) ~= 0
    error('Two lists do not match.')
end

% Obtain FMRP.
fmrpVec = demoTable.fmrp_replicate_trim_mean;




%%
subjNames  = cellfun(@(x) x(1:4), subjNames,         'UniformOutput', false);
groupNames = cellfun(@(x) x(),    sexGroupList(:,2), 'UniformOutput', false);
    %{
    figure
    bar(sort(ageVec))
    %}

fxsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'FXS'));
fxsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'FXS'));
hcsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'TDC'));
hcsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'TDC'));
fxsIdx   = sort([fxsM_Idx; fxsF_Idx]);
hcsIdx   = sort([hcsM_Idx; hcsF_Idx]);

% Obtain EEG IDs.
fxsM_Names = subjNames(fxsM_Idx);
fxsF_Names = subjNames(fxsF_Idx);
hcsM_Names = subjNames(hcsM_Idx);
hcsF_Names = subjNames(hcsF_Idx);
fxsNames   = subjNames(fxsIdx);
hcsNames   = subjNames(hcsIdx);

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

thetaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
gammaIdx = freqIdx7;

timeBins    = -1000:4000;
baselineIdx = find(timeBins>=-1000 & timeBins<=0);
onsetIdx    = find(timeBins>0 & timeBins<=500);
assrIdx     = find(timeBins>0 & timeBins<=3000);
offsetIdx   = find(timeBins>3000 & timeBins<=3500);
afterZeroIdx = find(timeBins>0);


%% Load Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');
allMats = allMats(preselectingIdx);
stackData_ERSP = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData_ERSP(:,:,:,matIdx) = elecErspMedian;
end
stackData_ERSP = 10*log10(stackData_ERSP);
stackData_base = bsxfun(@minus, stackData_ERSP, mean(stackData_ERSP(:,:,1:1000,:),3));
erspFXS_M      = stackData_base(:,:,:,fxsM_Idx);
erspFXS_F      = stackData_base(:,:,:,fxsF_Idx);
erspHCS_M      = stackData_base(:,:,:,hcsM_Idx);
erspHCS_F      = stackData_base(:,:,:,hcsF_Idx);
erspFXS        = stackData_base(:,:,:,fxsIdx);
erspHCS        = stackData_base(:,:,:,hcsIdx);


%% Calculate ERSSP. This takes 10 min.

% % FXS.
% pTensor     = single(zeros(size(erspFXS,1), size(erspFXS,2), size(erspFXS,3)));
% ersspTensor = single(zeros(size(erspFXS,2), size(erspFXS,3), size(erspFXS,4)));
% for subjIdx = 1:size(erspFXS, 4)
%     currentData = erspFXS(:,:,:,subjIdx);
%     for chIdx = 1:size(currentData, 1)
%         for freqIdxIdx = 1:size(currentData,2)
%             currentFreqP = squeeze(normcdf(currentData(chIdx,freqIdxIdx,:), squeeze(mean(currentData(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(currentData(chIdx,freqIdxIdx,baselineIdx),0,3))));
%             pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
%         end
%     end
%     pTensor_rTale   = 1-pTensor;
%     pTensor_sigMask = pTensor_rTale<0.05;
%     pTensorSumMasks = squeeze(sum(pTensor_sigMask,1));
%     ersspTensor(:,:,subjIdx) = pTensorSumMasks;
% end
% save('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_ersp_FXS', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')
% 
% 
% % HCS.
% pTensor     = single(zeros(size(erspHCS,1), size(erspHCS,2), size(erspHCS,3)));
% ersspTensor = single(zeros(size(erspHCS,2), size(erspHCS,3), size(erspHCS,4)));
% for subjIdx = 1:size(erspHCS, 4)
%     currentData = erspHCS(:,:,:,subjIdx);
%     for chIdx = 1:size(currentData, 1)
%         for freqIdxIdx = 1:size(currentData,2)
%             currentFreqP = squeeze(normcdf(currentData(chIdx,freqIdxIdx,:), squeeze(mean(currentData(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(currentData(chIdx,freqIdxIdx,baselineIdx),0,3))));
%             pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
%         end
%     end
%     pTensor_rTale   = 1-pTensor;
%     pTensor_sigMask = pTensor_rTale<0.05;
%     pTensorSumMasks = squeeze(sum(pTensor_sigMask,1));
%     ersspTensor(:,:,subjIdx) = pTensorSumMasks;
% end
% save('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_ersp_HCS', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')


%% Load Median ITC.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
allMats = allMats(preselectingIdx);
stackData_ITC = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData_ITC(:,:,:,matIdx) = elecItcMedian;
end
itcFXS_M      = stackData_ITC(:,:,:,fxsM_Idx);
itcFXS_F      = stackData_ITC(:,:,:,fxsF_Idx);
itcHCS_M      = stackData_ITC(:,:,:,hcsM_Idx);
itcHCS_F      = stackData_ITC(:,:,:,hcsF_Idx);
itcFXS        = stackData_ITC(:,:,:,fxsIdx);
itcHCS        = stackData_ITC(:,:,:,hcsIdx);


%% Calculate ERSSP. This takes 10 min.

% % FXS.
% pTensor     = single(zeros(size(itcFXS,1), size(itcFXS,2), size(itcFXS,3)));
% ersspTensor = single(zeros(size(itcFXS,2), size(itcFXS,3), size(itcFXS,4)));
% for subjIdx = 1:size(itcFXS, 4)
%     currentData = itcFXS(:,:,:,subjIdx);
%     for chIdx = 1:size(currentData, 1)
%         for freqIdxIdx = 1:size(currentData,2)
%             currentFreqP = squeeze(normcdf(currentData(chIdx,freqIdxIdx,:), squeeze(mean(currentData(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(currentData(chIdx,freqIdxIdx,baselineIdx),0,3))));
%             pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
%         end
%     end
%     pTensor_rTale   = 1-pTensor;
%     pTensor_sigMask = pTensor_rTale<0.05;
%     pTensorSumMasks = squeeze(sum(pTensor_sigMask,1));
%     ersspTensor(:,:,subjIdx) = pTensorSumMasks;
% end
% save('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_itc_FXS', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')
% 
% 
% % HCS.
% pTensor     = single(zeros(size(itcHCS,1), size(itcHCS,2), size(itcHCS,3)));
% ersspTensor = single(zeros(size(itcHCS,2), size(itcHCS,3), size(itcHCS,4)));
% for subjIdx = 1:size(itcHCS, 4)
%     currentData = itcHCS(:,:,:,subjIdx);
%     for chIdx = 1:size(currentData, 1)
%         for freqIdxIdx = 1:size(currentData,2)
%             currentFreqP = squeeze(normcdf(currentData(chIdx,freqIdxIdx,:), squeeze(mean(currentData(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(currentData(chIdx,freqIdxIdx,baselineIdx),0,3))));
%             pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
%         end
%     end
%     pTensor_rTale   = 1-pTensor;
%     pTensor_sigMask = pTensor_rTale<0.05;
%     pTensorSumMasks = squeeze(sum(pTensor_sigMask,1));
%     ersspTensor(:,:,subjIdx) = pTensorSumMasks;
% end
% save('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_itc_HCS', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')


%%
FXS_ersp = load('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_ersp_FXS.mat');
HCS_ersp = load('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_ersp_HCS.mat');
FXS_ersp = FXS_ersp.ersspTensor;
HCS_ersp = HCS_ersp.ersspTensor;

FXS_itc = load('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_itc_FXS.mat');
HCS_itc = load('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_itc_HCS.mat');
FXS_itc = FXS_itc.ersspTensor;
HCS_itc = HCS_itc.ersspTensor;

fxs_demographics = [subjNames(fxsIdx)' sexGroupList(fxsIdx,:)];
hcs_demographics = [subjNames(hcsIdx)' sexGroupList(hcsIdx,:)];

%% Bonferroni correction mask?

% % Apply multiple comparison.
% FXS_M(FXS_M<=6) = 0;
% FXS_F(FXS_F<=6) = 0;
% HCS_M(HCS_M<=6) = 0;
% HCS_F(HCS_F<=6) = 0;


%% Plot grand-average ERSP-ERSSP.
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw ERSSP plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position', [50 50 1000 1200])
mainLayout = tiledlayout(3,2, 'TileSpacing', 'tight', 'Padding', 'tight');

% Plot main ERSSP plots.
nexttile(mainLayout,1)
imagesc(timeBins, 1:100, mean(FXS_ersp,3), [0 55])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')

colormap jet
axis xy
axis square
title(sprintf('FXS Power (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,2)
imagesc(timeBins, 1:100, mean(HCS_ersp,3), [0 55])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
colormap jet
axis xy
axis square
title(sprintf('HCS Power (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,3)
imagesc(timeBins, 1:100, mean(FXS_itc,3), [0 85])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
colormap jet
axis xy
axis square
title(sprintf('FXS ITC (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,4)
imagesc(timeBins, 1:100, mean(HCS_itc,3), [0 85])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
colormap jet
axis xy
axis square
title(sprintf('HCS ITC (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')


%% Transplant p1100 here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for Onset's num cells. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
[~,gammaIdx]    = min(abs(wtFreqBins-40));
assrTimeIdx     = find(timeBins>0 & timeBins<=3000);
onsetVpTimeIdx  = find(timeBins>0 & timeBins<=500);
offsetVpTimeIdx = find(timeBins>3000 & timeBins<=3500);

fxsErsp = zeros(size(FXS_ersp,3),3); % Onset, ASSR, Offset.
for dataIdx = 1:size(FXS_ersp,3)
    currentData = FXS_ersp(:,:,dataIdx);
    currentMeanOnset  = mean(mean(currentData(deltaThetaAlphaIdx, onsetVpTimeIdx)));
    currentMeanAssr   = mean(mean(currentData(gammaIdx,           assrTimeIdx)));
    currentMeanOffset = mean(mean(currentData(deltaThetaAlphaIdx, offsetVpTimeIdx)));
    fxsErsp(dataIdx,:) = [currentMeanOnset currentMeanAssr currentMeanOffset];
end

hcsErsp = zeros(size(HCS_ersp,3),3); % Onset, ASSR, Offset.
for dataIdx = 1:size(HCS_ersp,3)
    currentData = HCS_ersp(:,:,dataIdx);
    currentMeanOnset  = mean(mean(currentData(deltaThetaAlphaIdx, onsetVpTimeIdx)));
    currentMeanAssr   = mean(mean(currentData(gammaIdx,           assrTimeIdx)));
    currentMeanOffset = mean(mean(currentData(deltaThetaAlphaIdx, offsetVpTimeIdx)));
    hcsErsp(dataIdx,:) = [currentMeanOnset currentMeanAssr currentMeanOffset];
end

fxsItc = zeros(size(FXS_itc,3),3); % Onset, ASSR, Offset.
for dataIdx = 1:size(FXS_itc,3)
    currentData = FXS_itc(:,:,dataIdx);
    currentMeanOnset  = mean(mean(currentData(deltaThetaAlphaIdx, onsetVpTimeIdx)));
    currentMeanAssr   = mean(mean(currentData(gammaIdx,           assrTimeIdx)));
    currentMeanOffset = mean(mean(currentData(deltaThetaAlphaIdx, offsetVpTimeIdx)));
    fxsItc(dataIdx,:) = [currentMeanOnset currentMeanAssr currentMeanOffset];
end

hcsItc = zeros(size(HCS_itc,3),3); % Onset, ASSR, Offset.
for dataIdx = 1:size(HCS_itc,3)
    currentData = HCS_itc(:,:,dataIdx);
    currentMeanOnset  = mean(mean(currentData(deltaThetaAlphaIdx, onsetVpTimeIdx)));
    currentMeanAssr   = mean(mean(currentData(gammaIdx,           assrTimeIdx)));
    currentMeanOffset = mean(mean(currentData(deltaThetaAlphaIdx, offsetVpTimeIdx)));
    hcsItc(dataIdx,:) = [currentMeanOnset currentMeanAssr currentMeanOffset];
end


%% Stats.

[H1,P1,CI1,STATS1] = ttest2(fxsErsp, hcsErsp);
d1 = cohensD(fxsErsp, hcsErsp);
disp(sprintf('Onset ERSP:  t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS1.df(1), STATS1.tstat(1), P1(1), d1(1)))
disp(sprintf('ASSR ERSP:   t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS1.df(2), STATS1.tstat(2), P1(2), d1(2)))
disp(sprintf('Offset ERSP: t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS1.df(3), STATS1.tstat(3), P1(3), d1(3)))
    % Onset ERSP:  t(60)=3.573, p=0.001, Cohen's d = 0.909.
    % ASSR ERSP:   t(60)=-0.889, p=0.378, Cohen's d = -0.226.
    % Offset ERSP: t(60)=1.754, p=0.085, Cohen's d = 0.446.

[H2,P2,CI2,STATS2] = ttest2(fxsItc, hcsItc);
d2 = cohensD(fxsItc, hcsItc);
disp(sprintf('Onset ITC:  t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS2.df(1), STATS2.tstat(1), P2(1), d2(1)))
disp(sprintf('ASSR ITC:   t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS2.df(2), STATS2.tstat(2), P2(2), d2(2)))
disp(sprintf('Offset ITC: t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS2.df(3), STATS2.tstat(3), P2(3), d2(3)))
    % Onset ITC:  t(60)=1.550, p=0.126, Cohen's d = 0.395.
    % ASSR ITC:   t(60)=-2.895, p=0.005, Cohen's d = -0.737.
    % Offset ITC: t(60)=1.803, p=0.076, Cohen's d = 0.459.


%% Plot bargraphs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for ERSP. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot bar graphs at 5 of (2,6), 1 of (1,3)
nestedLayout1 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 5;
nexttile(nestedLayout1, 1)
barHandle = bar([mean(fxsErsp(:,1)) mean(hcsErsp(:,1))]);
hold on
errorbar(1:2, [mean(fxsErsp(:,1)) mean(hcsErsp(:,1))], [0 0], [std(fxsErsp(:,1))/sqrt(length(fxsErsp)) std(hcsErsp(:,1))/sqrt(length(fxsErsp))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 40])
title('OnsetER')
set(gca, 'box', 'off')

% Plot bar graphs at 5 of (2,6), 2 of (1,3)
nexttile(nestedLayout1, 2)
barHandle = bar([mean(fxsErsp(:,2)) mean(hcsErsp(:,2))]);
hold on
errorbar(1:2, [mean(fxsErsp(:,2)) mean(hcsErsp(:,2))], [0 0], [std(fxsErsp(:,2))/sqrt(length(fxsErsp)) std(hcsErsp(:,2))/sqrt(length(fxsErsp))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 40])
title('ASSR')
set(gca, 'box', 'off')

% Plot bar graphs at 5 of (2,6), 3 of (1,3)
nexttile(nestedLayout1, 3)
barHandle = bar([mean(fxsErsp(:,3)) mean(hcsErsp(:,3))]);
hold on
errorbar(1:2, [mean(fxsErsp(:,3)) mean(hcsErsp(:,3))], [0 0], [std(fxsErsp(:,3))/sqrt(length(fxsErsp)) std(hcsErsp(:,3))/sqrt(length(fxsErsp))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 40])
title('OffsetER')
set(gca, 'box', 'off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for ITC. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot bar graphs at 6 of (2,6), 1 of (1,3)
nestedLayout2 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 6;
nexttile(nestedLayout2, 1)
barHandle = bar([mean(fxsItc(:,1)) mean(hcsItc(:,1))]);
hold on
errorbar(1:2, [mean(fxsItc(:,1)) mean(hcsItc(:,1))], [0 0], [std(fxsItc(:,1))/sqrt(length(fxsItc)) std(hcsItc(:,1))/sqrt(length(fxsItc))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 80])
title('OnsetER')
set(gca, 'box', 'off')

% Plot bar graphs at 6 of (2,6), 2 of (1,3)
nestedLayout2 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 6;
nexttile(nestedLayout2, 2)
barHandle = bar([mean(fxsItc(:,2)) mean(hcsItc(:,2))]);
hold on
errorbar(1:2, [mean(fxsItc(:,2)) mean(hcsItc(:,2))], [0 0], [std(fxsItc(:,2))/sqrt(length(fxsItc)) std(hcsItc(:,2))/sqrt(length(fxsItc))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 80])
title('ASSR')
set(gca, 'box', 'off')

% Plot bar graphs at 6 of (2,6), 3 of (1,3)
nestedLayout2 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 6;
nexttile(nestedLayout2, 3)
barHandle = bar([mean(fxsItc(:,3)) mean(hcsItc(:,3))]);
hold on
errorbar(1:2, [mean(fxsItc(:,3)) mean(hcsItc(:,3))], [0 0], [std(fxsItc(:,3))/sqrt(length(fxsItc)) std(hcsItc(:,3))/sqrt(length(fxsItc))], 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS' 'HCS'})
ylim([0 80])
title('OffsetER')
set(gca, 'box', 'off')


% Export the results.
combinedDemographics = [fxs_demographics; hcs_demographics];
fmrp     = [fmrpVec(fxsIdx); fmrpVec(hcsIdx)];
subjID   = combinedDemographics(:,1);
sex      = combinedDemographics(:,2);
group    = combinedDemographics(:,3);
powerOnsetER  = [fxsErsp(:,1); hcsErsp(:,1)];
powerASSR     = [fxsErsp(:,2); hcsErsp(:,2)];
powerOffsetER = [fxsErsp(:,3); hcsErsp(:,3)];
itcOnsetER    = [fxsItc(:,1);  hcsItc(:,1)];
itcASSR       = [fxsItc(:,2);  hcsItc(:,2)];
itcOffsetER   = [fxsItc(:,3);  hcsItc(:,3)];
resultTableErssp = table(subjID, sex, group, powerOnsetER, powerASSR, powerOffsetER, itcOnsetER, itcASSR, itcOffsetER);

exportgraphics(gcf, '/srv/Makoto/ASSR/p1600_ERSSP_erspItc/ERSSP_erspItc.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')
save('/srv/Makoto/ASSR/p1600_ERSSP_erspItc/resultTableErssp.mat', 'resultTableErssp')