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

% Cleaning stats.
mean(cleaningMatrix,1)  % 0.0441   -2.5726    1.8670   -9.9652   -0.5785
std(cleaningMatrix,0,1) % 0.2069    1.8581    0.9892    5.9418    0.6645
min(cleaningMatrix)     %      0   -8.0536    0.1790  -29.8985   -2.7451   
max(cleaningMatrix)     % 1.0000   -0.1294    5.0535   -0.9597    0.1880


%%
subjNames  = cellfun(@(x) x(1:4), {allSets.name},    'UniformOutput', false);
groupNames = cellfun(@(x) x(),    sexGroupList(:,2), 'UniformOutput', false);
    %{
    figure
    bar(sort(ageVec))
    %}

fxsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'FXS'));
fxsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'FXS'));
hcsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'TDC'));
hcsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'TDC'));
fxsIdx = sort([fxsM_Idx; fxsF_Idx]);
hcsIdx = sort([hcsM_Idx; hcsF_Idx]);

% Obtain EEG IDs.
fxsM_Names = subjNames(fxsM_Idx);
fxsF_Names = subjNames(fxsF_Idx);
hcsM_Names = subjNames(hcsM_Idx);
hcsF_Names = subjNames(hcsF_Idx);
fxsNames = subjNames(fxsIdx);
hcsNames = subjNames(hcsIdx);

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


%%
% Load Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');
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

%% Calculate ERSSP. This takes 10 min.

% % FXS_M.
% pTensor     = single(zeros(size(erspFXS_M,1), size(erspFXS_M,2), size(erspFXS_M,3)));
% ersspTensor = single(zeros(size(erspFXS_M,2), size(erspFXS_M,3), size(erspFXS_M,4)));
% for subjIdx = 1:size(erspFXS_M, 4)
%     currentData = erspFXS_M(:,:,:,subjIdx);
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
% save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_FXS_M', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')
% 
% % FXS_F.
% pTensor     = single(zeros(size(erspFXS_F,1), size(erspFXS_F,2), size(erspFXS_F,3)));
% ersspTensor = single(zeros(size(erspFXS_F,2), size(erspFXS_F,3), size(erspFXS_F,4)));
% for subjIdx = 1:size(erspFXS_F, 4)
%     currentData = erspFXS_F(:,:,:,subjIdx);
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
% save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_FXS_F', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')
% 
% % HCS_M.
% pTensor     = single(zeros(size(erspHCS_M,1), size(erspHCS_M,2), size(erspHCS_M,3)));
% ersspTensor = single(zeros(size(erspHCS_M,2), size(erspHCS_M,3), size(erspHCS_M,4)));
% for subjIdx = 1:size(erspHCS_M, 4)
%     currentData = erspHCS_M(:,:,:,subjIdx);
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
% save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_HCS_M', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')
% 
% % HCS_F.
% pTensor     = single(zeros(size(erspHCS_F,1), size(erspHCS_F,2), size(erspHCS_F,3)));
% ersspTensor = single(zeros(size(erspHCS_F,2), size(erspHCS_F,3), size(erspHCS_F,4)));
% for subjIdx = 1:size(erspHCS_F, 4)
%     currentData = erspHCS_F(:,:,:,subjIdx);
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
% save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_HCS_F', 'ersspTensor', 'timeBins', 'wtFreqBins', '-nocompression')


%%
FXS_M = load('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_FXS_M.mat');
FXS_F = load('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_FXS_F.mat');
HCS_M = load('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_HCS_M.mat');
HCS_F = load('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ersp_HCS_F.mat');
FXS_M = FXS_M.ersspTensor;
FXS_F = FXS_F.ersspTensor;
HCS_M = HCS_M.ersspTensor;
HCS_F = HCS_F.ersspTensor;

% clear erspFXS_* erspHCS_* stackData_ERSP stackData_base
% save dumpedFiles

%%

% % Apply multiple comparison.
% FXS_M(FXS_M<=6) = 0;
% FXS_F(FXS_F<=6) = 0;
% HCS_M(HCS_M<=6) = 0;
% HCS_F(HCS_F<=6) = 0;

%% ERSSP plot with peak masks.

gmFXS_M  = mean(FXS_M,3);
gmFXS_F  = mean(FXS_F,3);
%gmFXS_MF = gmFXS_M-gmFXS_F;
gmHCS_M  = mean(HCS_M,3);
gmHCS_F  = mean(HCS_F,3);
%gmHCS_MF = gmHCS_M-gmHCS_F;

gmFXS_M(gmFXS_M<prctile(gmFXS_M(:), 95))= 0;
gmFXS_F(gmFXS_F<prctile(gmFXS_F(:), 95))= 0;
%gmFXS_MF(gmFXS_MF<prctile(gmFXS_MF(:), 95))= 0;
gmHCS_M(gmHCS_M<prctile(gmHCS_M(:), 95))= 0;
gmHCS_F(gmHCS_F<prctile(gmHCS_F(:), 95))= 0;
%gmHCS_MF(gmHCS_MF<prctile(gmHCS_MF(:), 95))= 0;

gmFXS_M(gmFXS_M>0)= 1;
gmFXS_F(gmFXS_F>0)= 1;
%gmFXS_MF(gmFXS_MF>0)= 1;
gmHCS_M(gmHCS_M>0)= 1;
gmHCS_F(gmHCS_F>0)= 1;
%gmHCS_MF(gmHCS_MF>0)= 1;

lowDeltaIdx = find(wtFreqBins<=2);
gmFXS_M(lowDeltaIdx,:)=0;
gmFXS_F(lowDeltaIdx,:)=0;
%gmFXS_MF(lowDeltaIdx,:)=0;
gmHCS_M(lowDeltaIdx,:)=0;
gmHCS_F(lowDeltaIdx,:)=0;
%gmHCS_MF(lowDeltaIdx,:)=0;


%% Plot grand-average ERSP-ERSSP.
%%%%%%%%%%%%%%%%%%%%%%%
%%% Tailor the mask.%%%
%%%%%%%%%%%%%%%%%%%%%%%
deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
nullMask = zeros(size(gmFXS_M));

% FXS_M.
currentData = gmFXS_M;
currentData(currentData<prctile(currentData(:), 95))= 0;
inclusionMask = nullMask;
inclusionMask(deltaThetaAlphaIdx, [onsetIdx offsetIdx]) = 1;
%inclusionMask(freqIdx7, assrIdx) = 1;
gmFXS_M_masked = logical(currentData.*inclusionMask);

% FXS_F.
currentData = gmFXS_F;
currentData(currentData<prctile(currentData(:), 95))= 0;
inclusionMask = nullMask;
inclusionMask(deltaThetaAlphaIdx, [onsetIdx offsetIdx]) = 1;
%inclusionMask(freqIdx7, assrIdx) = 1;
gmFXS_F_masked = logical(currentData.*inclusionMask);

% HCS_M.
currentData = gmHCS_M;
currentData(currentData<prctile(currentData(:), 95))= 0;
inclusionMask = nullMask;
inclusionMask(deltaThetaAlphaIdx, [onsetIdx offsetIdx]) = 1;
%inclusionMask(freqIdx7, assrIdx) = 1;
gmHCS_M_masked = logical(currentData.*inclusionMask);

% HCS_F.
currentData = gmHCS_F;
currentData(currentData<prctile(currentData(:), 95))= 0;
inclusionMask = nullMask;
inclusionMask(deltaThetaAlphaIdx, [onsetIdx offsetIdx]) = 1;
%inclusionMask(freqIdx7, assrIdx) = 1;
gmHCS_F_masked = logical(currentData.*inclusionMask);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw ERSSP plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position', [50 50 1000 1200])
mainLayout = tiledlayout(3,2, 'TileSpacing', 'tight', 'Padding', 'tight');

% Plot main ERSSP plots.
nexttile(mainLayout,1)
imagesc(timeBins, 1:100, mean(FXS_M,3), [0 55])
hold on
contour(timeBins, 1:100, logical(gmFXS_M_masked),  '-k')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS Male (n=%d)', length(fxsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,2)
imagesc(timeBins, 1:100, mean(FXS_F,3), [0 55])
hold on
contour(timeBins, 1:100, logical(gmFXS_F_masked),  '-k')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS Female (n=%d)', length(fxsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,3)
imagesc(timeBins, 1:100, mean(HCS_M,3), [0 55])
hold on
contour(timeBins, 1:100, logical(gmHCS_M_masked),  '-k')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS Male (n=%d)', length(hcsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')

nexttile(mainLayout,4)
imagesc(timeBins, 1:100, mean(HCS_F,3), [0 55])
hold on
contour(timeBins, 1:100, logical(gmHCS_F_masked),  '-k')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS Female (n=%d)', length(hcsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Num ch')


% Transplant p1100 here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for Onset's num cells. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);

fxsmSumCellOnsetVec = zeros(size(FXS_M,3),1);
fxsmSumChanOnsetVec = zeros(size(FXS_M,3),1);
for dataIdx = 1:size(FXS_M,3)
    currentData = FXS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, onsetIdx) = 1;
    gmFXS_M_masked = currentData.*inclusionMask;
    fxsmSumCellOnsetVec(dataIdx) = sum(sum(logical(gmFXS_M_masked)));
    fxsmSumChanOnsetVec(dataIdx) = sum(sum(gmFXS_M_masked))/sum(sum(logical(gmFXS_M_masked)));
end

fxsfSumCellOnsetVec = zeros(size(FXS_F,3),1);
fxsfSumChanOnsetVec = zeros(size(FXS_F,3),1);
for dataIdx = 1:size(FXS_F,3)
    currentData = FXS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, onsetIdx) = 1;
    gmFXS_F_masked = currentData.*inclusionMask;
    fxsfSumCellOnsetVec(dataIdx) = sum(sum(logical(gmFXS_F_masked)));
    fxsfSumChanOnsetVec(dataIdx) = sum(sum(gmFXS_F_masked))/sum(sum(logical(gmFXS_F_masked)));
end

hcsmSumCellOnsetVec = zeros(size(HCS_M,3),1);
hcsmSumChanOnsetVec = zeros(size(HCS_M,3),1);
for dataIdx = 1:size(HCS_M,3)
    currentData = HCS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, onsetIdx) = 1;
    gmHCS_M_masked = currentData.*inclusionMask;
    hcsmSumCellOnsetVec(dataIdx) = sum(sum(logical(gmHCS_M_masked)));
    hcsmSumChanOnsetVec(dataIdx) = sum(sum(gmHCS_M_masked))/sum(sum(logical(gmHCS_M_masked)));
end

hcsfSumCellOnsetVec = zeros(size(HCS_F,3),1);
hcsfSumChanOnsetVec = zeros(size(HCS_F,3),1);
for dataIdx = 1:size(HCS_F,3)
    currentData = HCS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, onsetIdx) = 1;
    gmHCS_F_masked = currentData.*inclusionMask;
    hcsfSumCellOnsetVec(dataIdx) = sum(sum(logical(gmHCS_F_masked)));
    hcsfSumChanOnsetVec(dataIdx) = sum(sum(gmHCS_F_masked))/sum(sum(logical(gmHCS_F_masked)));
end

% Convert NaN to 0.
fxsmSumChanOnsetVec(isnan(fxsmSumChanOnsetVec)) = 0;
fxsfSumChanOnsetVec(isnan(fxsfSumChanOnsetVec)) = 0;
hcsmSumChanOnsetVec(isnan(hcsmSumChanOnsetVec)) = 0;
hcsfSumChanOnsetVec(isnan(hcsfSumChanOnsetVec)) = 0;

% Perform stats.
[numCellStats_onset, numCellDf_onset, numCellPvals_onset] = statcond({fxsmSumCellOnsetVec' fxsfSumCellOnsetVec'; hcsmSumCellOnsetVec' hcsfSumCellOnsetVec'});
numCellOnsetMean = [mean(fxsmSumCellOnsetVec) mean(fxsfSumCellOnsetVec) mean(hcsmSumCellOnsetVec) mean(hcsfSumCellOnsetVec)];
numCellOnsetSe   = [std(fxsmSumCellOnsetVec)/sqrt(length(fxsmSumCellOnsetVec)) std(fxsfSumCellOnsetVec)/sqrt(length(fxsfSumCellOnsetVec)) std(hcsmSumCellOnsetVec)/sqrt(length(hcsmSumCellOnsetVec)) std(hcsfSumCellOnsetVec)/sqrt(length(hcsfSumCellOnsetVec))];
numChanOnsetMean = [mean(fxsmSumChanOnsetVec) mean(fxsfSumChanOnsetVec) mean(hcsmSumChanOnsetVec) mean(hcsfSumChanOnsetVec)];
numChanOnsetSe   = [std(fxsmSumChanOnsetVec)/sqrt(length(fxsmSumChanOnsetVec)) std(fxsfSumChanOnsetVec)/sqrt(length(fxsfSumChanOnsetVec)) std(hcsmSumChanOnsetVec)/sqrt(length(hcsmSumChanOnsetVec)) std(hcsfSumChanOnsetVec)/sqrt(length(hcsfSumChanOnsetVec))];
    % T  {[4.7606]}    {[10.5698]}    {[4.1548]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P {[0.0318]}    {[0.0016]}    {[0.0445]}

[numCellOnset_d_MaleFemale, numCellOnset_d_FxsHcs] = cohensD_mainEffects2x2(fxsmSumCellOnsetVec', fxsfSumCellOnsetVec', hcsmSumCellOnsetVec', hcsfSumCellOnsetVec');
numCellOnset_F2                                    = cohensF2(numCellStats_onset{3}, numCellDf_onset{3}(1), numCellDf_onset{3}(2));
disp(sprintf('%.3f %.3f %.3f',numCellOnset_d_MaleFemale, numCellOnset_d_FxsHcs, numCellOnset_F2))
    % Cohen 0.370 0.457 0.047

% Plot bar graphs at 5 of (2,6), 1 of (1,3)
nestedLayout1 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 5;
nexttile(nestedLayout1, 1)
barHandle = bar(numCellOnsetMean);
hold on
errorbar(1:4, numCellOnsetMean, [0 0 0 0], numCellOnsetSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
ylim([1000 7000])
title('OnsetER')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for ASSR's num cells. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot bar graphs at 5 of (2,6), 2 of (1,3)
nexttile(nestedLayout1, 2)
barHandle = bar([3000 3000 3000 3000]);

barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
ylim([1000 7000])
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('ASSR')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for Offset's num cells. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fxsmSumCellOffsetVec = zeros(size(FXS_M,3),1);
fxsmSumChanOffsetVec = zeros(size(FXS_M,3),1);
for dataIdx = 1:size(FXS_M,3)
    currentData = FXS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, offsetIdx) = 1;
    gmFXS_M_masked = currentData.*inclusionMask;
    fxsmSumCellOffsetVec(dataIdx) = sum(sum(logical(gmFXS_M_masked)));
    fxsmSumChanOffsetVec(dataIdx) = sum(sum(gmFXS_M_masked))/sum(sum(logical(gmFXS_M_masked)));
end

fxsfSumCellOffsetVec = zeros(size(FXS_F,3),1);
fxsfSumChanOffsetVec = zeros(size(FXS_F,3),1);
for dataIdx = 1:size(FXS_F,3)
    currentData = FXS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, offsetIdx) = 1;
    gmFXS_F_masked = currentData.*inclusionMask;
    fxsfSumCellOffsetVec(dataIdx) = sum(sum(logical(gmFXS_F_masked)));
    fxsfSumChanOffsetVec(dataIdx) = sum(sum(gmFXS_F_masked))/sum(sum(logical(gmFXS_F_masked)));
end

hcsmSumCellOffsetVec = zeros(size(HCS_M,3),1);
hcsmSumChanOffsetVec = zeros(size(HCS_M,3),1);
for dataIdx = 1:size(HCS_M,3)
    currentData = HCS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, offsetIdx) = 1;
    gmHCS_M_masked = currentData.*inclusionMask;
    hcsmSumCellOffsetVec(dataIdx) = sum(sum(logical(gmHCS_M_masked)));
    hcsmSumChanOffsetVec(dataIdx) = sum(sum(gmHCS_M_masked))/sum(sum(logical(gmHCS_M_masked)));
end

hcsfSumCellOffsetVec = zeros(size(HCS_F,3),1);
hcsfSumChanOffsetVec = zeros(size(HCS_F,3),1);
for dataIdx = 1:size(HCS_F,3)
    currentData = HCS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;
    inclusionMask = nullMask;
    inclusionMask(deltaThetaAlphaIdx, offsetIdx) = 1;
    gmHCS_F_masked = currentData.*inclusionMask;
    hcsfSumCellOffsetVec(dataIdx) = sum(sum(logical(gmHCS_F_masked)));
    hcsfSumChanOffsetVec(dataIdx) = sum(sum(gmHCS_F_masked))/sum(sum(logical(gmHCS_F_masked)));
end

% Convert NaN to 0.
fxsmSumChanOffsetVec(isnan(fxsmSumChanOffsetVec)) = 0;
fxsfSumChanOffsetVec(isnan(fxsfSumChanOffsetVec)) = 0;
hcsmSumChanOffsetVec(isnan(hcsmSumChanOffsetVec)) = 0;
hcsfSumChanOffsetVec(isnan(hcsfSumChanOffsetVec)) = 0;

% Perform stats.
[numCellStats_offset, numCellDf_offset, numCellPvals_offset] = statcond({fxsmSumCellOffsetVec' fxsfSumCellOffsetVec'; hcsmSumCellOffsetVec' hcsfSumCellOffsetVec'}); % FXS>HCS t(88)=19.9, p<0.001
numCellOffsetMean = [mean(fxsmSumCellOffsetVec) mean(fxsfSumCellOffsetVec) mean(hcsmSumCellOffsetVec) mean(hcsfSumCellOffsetVec)];
numCellOffsetSe   = [std(fxsmSumCellOffsetVec)/sqrt(length(fxsmSumCellOffsetVec)) std(fxsfSumCellOffsetVec)/sqrt(length(fxsfSumCellOffsetVec)) std(hcsmSumCellOffsetVec)/sqrt(length(hcsmSumCellOffsetVec)) std(hcsfSumCellOffsetVec)/sqrt(length(hcsfSumCellOffsetVec))];
numChanOffsetMean = [mean(fxsmSumChanOffsetVec) mean(fxsfSumChanOffsetVec) mean(hcsmSumChanOffsetVec) mean(hcsfSumChanOffsetVec)];
numChanOffsetSe   = [std(fxsmSumChanOffsetVec)/sqrt(length(fxsmSumChanOffsetVec)) std(fxsfSumChanOffsetVec)/sqrt(length(fxsfSumChanOffsetVec)) std(hcsmSumChanOffsetVec)/sqrt(length(hcsmSumChanOffsetVec)) std(hcsfSumChanOffsetVec)/sqrt(length(hcsfSumChanOffsetVec))];
    % T  {[1.4094]}    {[1.2350]}    {[5.7993e-04]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.2384]}    {[0.2695]}    {[0.9808]}

[numCellOffset_d_MaleFemale, numCellOffset_d_FxsHcs] = cohensD_mainEffects2x2(fxsmSumCellOffsetVec', fxsfSumCellOffsetVec', hcsmSumCellOffsetVec', hcsfSumCellOffsetVec');
numCellOffset_F2                                    = cohensF2(numCellStats_offset{3}, numCellDf_offset{3}(1), numCellDf_offset{3}(2));
disp(sprintf('%.3f %.3f %.3f',numCellOffset_d_MaleFemale, numCellOffset_d_FxsHcs, numCellOffset_F2))
    % Cohen 0.209 0.197 0.000

% Plot bar graphs at 5 of (2,6), 3 of (1,3)
nexttile(nestedLayout1, 3)
barHandle = bar(numCellOffsetMean);
hold on
errorbar(1:4, numCellOffsetMean, [0 0 0 0], numCellOffsetSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of t-f cells')
ylim([1000 7000])
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('OffsetER')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for Onset's num channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform stats.
[numChanStats_onset, numChanDf_onset, numChanPvals_onset] = statcond({fxsmSumChanOnsetVec' fxsfSumChanOnsetVec'; hcsmSumChanOnsetVec' hcsfSumChanOnsetVec'}); % FXS>HCS t(88)=19.9, p<0.001
    % T  {[7.2456]}    {[9.9459]}    {[1.5782]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.0085]}    {[0.0022]}    {[0.2123]}

[numChanOnset_d_MaleFemale, numChanOnset_d_FxsHcs] = cohensD_mainEffects2x2(fxsmSumChanOnsetVec', fxsfSumChanOnsetVec', hcsmSumChanOnsetVec', hcsfSumChanOnsetVec');
numChanOnset_F2                                    = cohensF2(numChanStats_onset{3}, numChanDf_onset{3}(1), numChanDf_onset{3}(2));
disp(sprintf('%.3f %.3f %.3f',numChanOnset_d_MaleFemale, numChanOnset_d_FxsHcs, numChanOnset_F2))
    % Cohen 0.480 0.629 0.018

% Plot bar graphs at 6 of (2,6), 3 of (1,3)
nestedLayout2 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 6;
nexttile(nestedLayout2, 1)
barHandle = bar(numChanOnsetMean);
hold on
errorbar(1:4, numChanOnsetMean, [0 0 0 0], numChanOnsetSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('OnsetER')
ylim([10 80])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for ASSR's num channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform stats.
fxsmSumChanAssrVec = squeeze(mean(FXS_M(gammaIdx,assrIdx,:),2));
fxsfSumChanAssrVec = squeeze(mean(FXS_F(gammaIdx,assrIdx,:),2));
hcsmSumChanAssrVec = squeeze(mean(HCS_M(gammaIdx,assrIdx,:),2));
hcsfSumChanAssrVec = squeeze(mean(HCS_F(gammaIdx,assrIdx,:),2));
[assrStats, assrDf, assrPvals] = statcond({fxsmSumChanAssrVec' fxsfSumChanAssrVec'; hcsmSumChanAssrVec' hcsfSumChanAssrVec'});
    % stats2 = {[0.0151]}    {[3.2381]}    {[0.3414]}
    % df2    = {[1 88]}    {[1 88]}    {[1 88]}
    % pvals2 = {[0.9023]}    {[0.0754]}    {[0.5605]}

[numChanAssr_d_MaleFemale, numChanAssr_d_FxsHcs] = cohensD_mainEffects2x2(fxsmSumChanAssrVec', fxsfSumChanAssrVec', hcsmSumChanAssrVec', hcsfSumChanAssrVec');
numChanAssr_F2                                   = cohensF2(assrStats{3}, assrDf{3}(1), assrDf{3}(2));
disp(sprintf('%.3f %.3f %.3f',numChanAssr_d_MaleFemale, numChanAssr_d_FxsHcs, numChanAssr_F2))
    % Cohen 0.029 0.287 0.004

nexttile(nestedLayout2, 2)
numChanAssrMean = [mean(fxsmSumChanAssrVec) mean(fxsfSumChanAssrVec) mean(hcsmSumChanAssrVec) mean(hcsfSumChanAssrVec)];
numChanAssrSe   = [std(fxsmSumChanAssrVec)/sqrt(length(fxsmSumChanAssrVec)) std(fxsfSumChanAssrVec)/sqrt(length(fxsfSumChanAssrVec))  std(hcsmSumChanAssrVec)/sqrt(length(hcsmSumChanAssrVec))  std(hcsfSumChanAssrVec)/sqrt(length(hcsfSumChanAssrVec))];
barHandle = bar(numChanAssrMean);
hold on
errorbar(1:4, numChanAssrMean, [0 0 0 0], numChanAssrSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('ASSR')
ylim([10 80])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw bargraphs for Offset's num channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform stats.
[numChanStats_offset, numChanDf_offset, numChanPvals_offset] = statcond({fxsmSumChanOffsetVec' fxsfSumChanOffsetVec'; hcsmSumChanOffsetVec' hcsfSumChanOffsetVec'}); % FXS>HCS t(88)=19.9, p<0.001
    % stats2 = {[0.6185]}    {[5.2262]}    {[7.9363]}
    % df2    = {[1 88]}    {[1 88]}    {[1 88]}
    % pvals2 = {[0.4337]}    {[0.0247]}    {[0.0060]}

[numChanOffset_d_MaleFemale, numChanOffset_d_FxsHcs] = cohensD_mainEffects2x2(fxsmSumChanOffsetVec', fxsfSumChanOffsetVec', hcsmSumChanOffsetVec', hcsfSumChanOffsetVec');
numChanOffset_F2                                     = cohensF2(numChanStats_offset{3}, numChanDf_offset{3}(1), numChanDf_offset{3}(2));
disp(sprintf('%.3f %.3f %.3f',numChanOffset_d_MaleFemale, numChanOffset_d_FxsHcs, numChanOffset_F2))
    % Cohen 0.152 0.567 0.090

nexttile(nestedLayout2, 3)
barHandle = bar(numChanOffsetMean);
hold on
errorbar(1:4, numChanOffsetMean, [0 0 0 0], numChanOffsetSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('OffsetER')
ylim([10 80])

exportgraphics(gcf, '/srv/Makoto/ASSR/p1200_maleFemale_ERSP_individualMask/ERSSP_ERSP_maleFemale.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')