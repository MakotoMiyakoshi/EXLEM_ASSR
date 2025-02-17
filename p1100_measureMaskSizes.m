% 12/04/2024 Makoto. Modified.
% 12/02/2024 Makoto. Used again.
% 11/29/2024 Makoto. Used again.
% 10/25/2024 Makoto. Used.
% 10/04/2024 Makoto. Modified.
% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all


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
figure('visible', 'off')
for dataIdx = 1:size(FXS_M,3)

    currentData = FXS_M(:,:,dataIdx);

    imagesc(timeBins, 1:100, currentData, [0 max(currentData(:))])
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    colormap jet
    axis xy
    axis square
    title(sprintf('FXS M %03d', dataIdx))
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
    line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    xlim([-800 3800])
    xlabel('Latency (ms)')
    ylabel('Frequency (Hz)')
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'NumCh')
    print(sprintf('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/FXS_M%03d', dataIdx), '-r200', '-djpeg95')
    cla
end

for dataIdx = 1:size(FXS_F,3)

    currentData = FXS_F(:,:,dataIdx);

    imagesc(timeBins, 1:100, currentData, [0 max(currentData(:))])
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    colormap jet
    axis xy
    axis square
    title(sprintf('FXS F %03d', dataIdx))
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
    line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    xlim([-800 3800])
    xlabel('Latency (ms)')
    ylabel('Frequency (Hz)')
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'NumCh')
    print(sprintf('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/FXS_F%03d', dataIdx), '-r200', '-djpeg95')
    cla
end

for dataIdx = 1:size(HCS_M,3)

    currentData = HCS_M(:,:,dataIdx);

    imagesc(timeBins, 1:100, currentData, [0 max(currentData(:))])
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    colormap jet
    axis xy
    axis square
    title(sprintf('HCS M %03d', dataIdx))
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
    line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    xlim([-800 3800])
    xlabel('Latency (ms)')
    ylabel('Frequency (Hz)')
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'NumCh')
    print(sprintf('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/HCS_M%03d', dataIdx), '-r200', '-djpeg95')
    cla
end

for dataIdx = 1:size(HCS_F,3)

    currentData = HCS_F(:,:,dataIdx);

    imagesc(timeBins, 1:100, currentData, [0 max(currentData(:))])
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    colormap jet
    axis xy
    axis square
    title(sprintf('HCS F %03d', dataIdx))
    set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
    line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
    line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    xlim([-800 3800])
    xlabel('Latency (ms)')
    ylabel('Frequency (Hz)')
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'NumCh')
    print(sprintf('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/HCS_F%03d', dataIdx), '-r200', '-djpeg95')
    cla
end


%%

%%%%%%%%%%%%%%
%%% Onset. %%%
%%%%%%%%%%%%%%
deltaThetaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);

fxsmVec = zeros(size(FXS_M,3),1);
for dataIdx = 1:size(FXS_M,3)

    currentData = FXS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, onsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    fxsmVec(dataIdx) = sum(sum(currentData));
end

fxsfVec = zeros(size(FXS_F,3),1);
for dataIdx = 1:size(FXS_F,3)

    currentData = FXS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, onsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    fxsfVec(dataIdx) = sum(sum(currentData));
end

hcsmVec = zeros(size(HCS_M,3),1);
for dataIdx = 1:size(HCS_M,3)

    currentData = HCS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, onsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    hcsmVec(dataIdx) = sum(sum(currentData));
end

hcsfVec = zeros(size(HCS_F,3),1);
for dataIdx = 1:size(HCS_F,3)

    currentData = HCS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, onsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    hcsfVec(dataIdx) = sum(sum(currentData));
end


figure
bar([mean(fxsmVec) mean(fxsfVec) mean(hcsmVec) mean(hcsfVec)])
[stats_onset, df_onset, pvals_onset] = statcond({fxsmVec' fxsfVec'; hcsmVec' hcsfVec'});

onsetMaskSizeMean = [mean(fxsmVec) mean(fxsfVec) mean(hcsmVec) mean(hcsfVec)];
onsetMaskSizeSe   = [std(fxsmVec)/sqrt(length(fxsmVec)) std(fxsfVec)/sqrt(length(fxsfVec)) std(hcsmVec)/sqrt(length(hcsmVec)) std(hcsfVec)/sqrt(length(hcsfVec))];

barHandle = bar(onsetMaskSizeMean);
hold on
errorbar(1:4, onsetMaskSizeMean, [0 0 0 0], onsetMaskSizeSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Onset numCh, FXS>HCS t(88)=19.9, p<0.001')

print('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/onset', '-r200', '-djpeg95')

%%
%%%%%%%%%%%%%%%
%%% Offset. %%%
%%%%%%%%%%%%%%%
fxsmVec = zeros(size(FXS_M,3),1);
for dataIdx = 1:size(FXS_M,3)

    currentData = FXS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, offsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    fxsmVec(dataIdx) = sum(sum(currentData));
end

fxsfVec = zeros(size(FXS_F,3),1);
for dataIdx = 1:size(FXS_F,3)

    currentData = FXS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, offsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    fxsfVec(dataIdx) = sum(sum(currentData));
end

hcsmVec = zeros(size(HCS_M,3),1);
for dataIdx = 1:size(HCS_M,3)

    currentData = HCS_M(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, offsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    hcsmVec(dataIdx) = sum(sum(currentData));
end

hcsfVec = zeros(size(HCS_F,3),1);
for dataIdx = 1:size(HCS_F,3)

    currentData = HCS_F(:,:,dataIdx);
    currentData(currentData<prctile(currentData(:), 95))= 0;

    % Cut out everything outside the time-freq of interest. 
    nonThetaIdx = setdiff(1:100, deltaThetaIdx);
    nonOnsetIdx = setdiff(1:5001, offsetIdx);
    currentData(nonThetaIdx, nonOnsetIdx) = 0;
    hcsfVec(dataIdx) = sum(sum(currentData));
end


figure
bar([mean(fxsmVec) mean(fxsfVec) mean(hcsmVec) mean(hcsfVec)])
[stats_offset, df_offset, pvals_offset] = statcond({fxsmVec' fxsfVec'; hcsmVec' hcsfVec'});

offsetMaskSizeMean = [mean(fxsmVec) mean(fxsfVec) mean(hcsmVec) mean(hcsfVec)];
offsetMaskSizeSe   = [std(fxsmVec)/sqrt(length(fxsmVec)) std(fxsfVec)/sqrt(length(fxsfVec)) std(hcsmVec)/sqrt(length(hcsmVec)) std(hcsfVec)/sqrt(length(hcsfVec))];

barHandle = bar(offsetMaskSizeMean);
hold on
errorbar(1:4, offsetMaskSizeMean, [0 0 0 0], offsetMaskSizeSe, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Offset numCh, FXS>HCS t(88)=12.7, p<0.005')

print('/srv/Makoto/ASSR/p1100_measureMaskSizes/images/offset', '-r200', '-djpeg95')
