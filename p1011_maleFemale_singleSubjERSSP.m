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


%%%%%%%%%%%%%%%%%%%%%%%
%%% Tailor the mask.%%%
%%%%%%%%%%%%%%%%%%%%%%%
% FXS_M.
[L,NUM] = bwlabel(gmFXS_M, 8);
sumPixelVec = 1:NUM;
for maskIdx = 1:NUM
    sumPixels = sum(sum(L==maskIdx));
    sumPixelVec(maskIdx) = sumPixels;
end
gmFXS_M_onsetMask  = L==1;
gmFXS_M_offsetMask = L==2;

% FXS_F.
[L,NUM] = bwlabel(gmFXS_F, 8);
sumPixelVec = 1:NUM;
for maskIdx = 1:NUM
    sumPixels = sum(sum(L==maskIdx));
    sumPixelVec(maskIdx) = sumPixels;
end
gmFXS_F_onsetMask  = L==1;
gmFXS_F_offsetMask = L==17;

% % FXS_MF.
% [L,NUM] = bwlabel(gmFXS_MF, 8);
% sumPixelVec = 1:NUM;
% for maskIdx = 1:NUM
%     sumPixels = sum(sum(L==maskIdx));
%     sumPixelVec(maskIdx) = sumPixels;
% end
% gmFXS_MF_onsetMask  = L==2;
% gmFXS_MF_offsetMask = L==322 | L==335;

% HCS_M.
[L,NUM] = bwlabel(gmHCS_M, 8);
sumPixelVec = 1:NUM;
for maskIdx = 1:NUM
    sumPixels = sum(sum(L==maskIdx));
    sumPixelVec(maskIdx) = sumPixels;
end
gmHCS_M_onsetMask  = L==2;
gmHCS_M_offsetMask = L==13;

% HCS_M.
[L,NUM] = bwlabel(gmHCS_F, 8);
sumPixelVec = 1:NUM;
for maskIdx = 1:NUM
    sumPixels = sum(sum(L==maskIdx));
    sumPixelVec(maskIdx) = sumPixels;
end
gmHCS_F_onsetMask  = L==2;
gmHCS_F_offsetMask = L==19 | L==23 | L==24 | L==25 | L==26 | L==27 | L==28 | L==29 | L==30 | L==31; 

% % HCS_MF.
% [L,NUM] = bwlabel(gmHCS_MF, 8);
% sumPixelVec = 1:NUM;
% for maskIdx = 1:NUM
%     sumPixels = sum(sum(L==maskIdx));
%     sumPixelVec(maskIdx) = sumPixels;
% end
% gmHCS_MF_onsetMask  = L==28;
% gmHCS_MF_offsetMask = L==369 | L==384;


%%
%%%%%%%%%%%%%%%%%%%
%%% Draw plots. %%%
%%%%%%%%%%%%%%%%%%%
figure('position', [50 50 1000 1200])
mainLayout = tiledlayout(3,2, 'TileSpacing', 'tight', 'Padding', 'tight');


% Plot main ERSSP plots.
nexttile(mainLayout,1)
imagesc(timeBins, 1:100, mean(FXS_M,3), [0 55])
hold on
contour(timeBins, 1:100, gmFXS_M_onsetMask,  '-k')
contour(timeBins, 1:100, gmFXS_M_offsetMask, '-k')
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
set(get(cbarHandle, 'title'), 'string', 'NumCh')

nexttile(mainLayout,2)
imagesc(timeBins, 1:100, mean(FXS_F,3), [0 55])
hold on
contour(timeBins, 1:100, gmFXS_F_onsetMask,  '-k')
contour(timeBins, 1:100, gmFXS_F_offsetMask, '-k')
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
set(get(cbarHandle, 'title'), 'string', 'NumCh')

nexttile(mainLayout,3)
imagesc(timeBins, 1:100, mean(HCS_M,3), [0 55])
hold on
contour(timeBins, 1:100, gmHCS_M_onsetMask,  '-k')
contour(timeBins, 1:100, gmHCS_M_offsetMask, '-k')
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
set(get(cbarHandle, 'title'), 'string', 'NumCh')

nexttile(mainLayout,4)
imagesc(timeBins, 1:100, mean(HCS_F,3), [0 55])
hold on
contour(timeBins, 1:100, gmHCS_F_onsetMask,  '-k')
contour(timeBins, 1:100, gmHCS_F_offsetMask, '-k')
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
set(get(cbarHandle, 'title'), 'string', 'NumCh')


%%
% Plot bar graphs 1.
nestedLayout1 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 5;

nexttile(nestedLayout1, 1)
onsetMaskSize  = [sum(sum(gmFXS_M_onsetMask))  sum(sum(gmFXS_F_onsetMask))  sum(sum(gmHCS_M_onsetMask))  sum(sum(gmHCS_F_onsetMask))];
barHandle = bar(onsetMaskSize);
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of time-freq cells')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Onset mask area')
ylim([0 12000])

nexttile(nestedLayout1, 3)
offsetMaskSize  = [sum(sum(gmFXS_M_offsetMask))  sum(sum(gmFXS_F_offsetMask))  sum(sum(gmHCS_M_offsetMask))  sum(sum(gmHCS_F_offsetMask))];
barHandle = bar(offsetMaskSize);
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of time-freq cells')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Offset mask area')
ylim([0 12000])














%%
% Plot bar graphs 2.
nestedLayout2 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 6;

nexttile(nestedLayout2, 1)
tmp = bsxfun(@times, FXS_M, gmFXS_M_onsetMask);
tmp(tmp==0)=NaN;
gmFXS_M_onset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, FXS_F, gmFXS_F_onsetMask);
tmp(tmp==0)=NaN;
gmFXS_F_onset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, HCS_M, gmHCS_M_onsetMask);
tmp(tmp==0)=NaN;
gmHCS_M_onset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, HCS_F, gmHCS_F_onsetMask);
tmp(tmp==0)=NaN;
gmHCS_F_onset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
individualErssp_onset_mean = [mean(gmFXS_M_onset_individuals) mean(gmFXS_F_onset_individuals) mean(gmHCS_M_onset_individuals) mean(gmHCS_F_onset_individuals)];
individualErssp_onset_se   = [std(gmFXS_M_onset_individuals)/sqrt(length(gmFXS_M_onset_individuals)) std(gmFXS_F_onset_individuals)/sqrt(length(gmFXS_F_onset_individuals))  std(gmHCS_M_onset_individuals)/sqrt(length(gmHCS_M_onset_individuals))  std(gmHCS_F_onset_individuals)/sqrt(length(gmHCS_F_onset_individuals))];
barHandle = bar(individualErssp_onset_mean);
hold on
errorbar(1:4, individualErssp_onset_mean, [0 0 0 0], individualErssp_onset_se, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Onset numCh')
ylim([20 50])

[stats1, df1, pvals1] = statcond({gmFXS_M_onset_individuals' gmFXS_F_onset_individuals'; gmHCS_M_onset_individuals' gmHCS_F_onset_individuals'});
% stats1 = {[0.9926]}    {[0.9913]}    {[0.1317]}
% df1    = {[1 88]}    {[1 88]}    {[1 88]}
% pvals1 = {[0.3218]}    {[0.3222]}    {[0.7176]}

[H,P,CI,STATS] = ttest2(gmFXS_M_onset_individuals, gmFXS_F_onset_individuals);
[H,P,CI,STATS] = ttest2(gmHCS_M_onset_individuals, gmHCS_F_onset_individuals);



nexttile(nestedLayout2, 2)
gmFXS_M_assr_individuals = squeeze(mean(FXS_M(gammaIdx,assrIdx,:),2));
gmFXS_F_assr_individuals = squeeze(mean(FXS_F(gammaIdx,assrIdx,:),2));
gmHCS_M_assr_individuals = squeeze(mean(HCS_M(gammaIdx,assrIdx,:),2));
gmHCS_F_assr_individuals = squeeze(mean(HCS_F(gammaIdx,assrIdx,:),2));
individualErssp_assr_mean = [mean(gmFXS_M_assr_individuals) mean(gmFXS_F_assr_individuals) mean(gmHCS_M_assr_individuals) mean(gmHCS_F_assr_individuals)];
individualErssp_assr_se  = [std(gmFXS_M_assr_individuals)/sqrt(length(gmFXS_M_assr_individuals)) std(gmFXS_F_assr_individuals)/sqrt(length(gmFXS_F_assr_individuals))  std(gmHCS_M_assr_individuals)/sqrt(length(gmHCS_M_assr_individuals))  std(gmHCS_F_assr_individuals)/sqrt(length(gmHCS_F_assr_individuals))];
barHandle = bar(individualErssp_assr_mean);
hold on
errorbar(1:4, individualErssp_assr_mean, [0 0 0 0], individualErssp_assr_se, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('ASSR numCh')
ylim([20 50])

[stats2, df2, pvals2] = statcond({gmFXS_M_assr_individuals' gmFXS_F_assr_individuals'; gmHCS_M_assr_individuals' gmHCS_F_assr_individuals'});
% stats2 = {[0.0151]}    {[3.2381]}    {[0.3414]}
% df2    = {[1 88]}    {[1 88]}    {[1 88]}
% pvals2 = {[0.9023]}    {[0.0754]}    {[0.5605]}


nexttile(nestedLayout2, 3)
tmp = bsxfun(@times, FXS_M, gmFXS_M_offsetMask);
tmp(tmp==0)=NaN;
gmFXS_M_offset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, FXS_F, gmFXS_F_offsetMask);
tmp(tmp==0)=NaN;
gmFXS_F_offset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, HCS_M, gmHCS_M_offsetMask);
tmp(tmp==0)=NaN;
gmHCS_M_offset_individuals = squeeze(nanmean(nanmean(tmp,1),2));
tmp = bsxfun(@times, HCS_F, gmHCS_F_offsetMask);
tmp(tmp==0)=NaN;
gmHCS_F_offset_individuals = squeeze(nanmean(nanmean(tmp,1),2));

individualErssp_offset_mean = [mean(gmFXS_M_offset_individuals) mean(gmFXS_F_offset_individuals) mean(gmHCS_M_offset_individuals) mean(gmHCS_F_offset_individuals)];
individualErssp_offset_se   = [std(gmFXS_M_offset_individuals)/sqrt(length(gmFXS_M_offset_individuals)) std(gmFXS_F_offset_individuals)/sqrt(length(gmFXS_F_offset_individuals))  std(gmHCS_M_offset_individuals)/sqrt(length(gmHCS_M_offset_individuals))  std(gmHCS_F_offset_individuals)/sqrt(length(gmHCS_F_offset_individuals))];

barHandle = bar(individualErssp_offset_mean);
hold on
errorbar(1:4, individualErssp_offset_mean, [0 0 0 0], individualErssp_offset_se, 'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Counts of channels')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('Offset numCh')
ylim([20 50])

[stats3, df3, pvals3] = statcond({gmFXS_M_offset_individuals' gmFXS_F_offset_individuals'; gmHCS_M_offset_individuals' gmHCS_F_offset_individuals'});
% stats3 = {[0.4231]}    {[0.0874]}    {[2.2727]}
% df3    = {[1 88]}    {[1 88]}    {[1 88]}
% pvals3 = {[0.5171]}    {[0.7682]}    {[0.1353]}

print('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ERSP_maleFemale', '-dsvg')
print('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_ERSP_maleFemale', '-djpeg95', '-r200')