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

%%

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

% allAbove13Idx = find(ageVec>12);
% subjNames_13plus  = subjNames( allAbove13Idx)';
% groupNames_13plus = groupNames(allAbove13Idx)';
% 
% % Load dummy data.
% EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');
% 
% % Load the demographic data.
% demographicData = readtable('/srv/Makoto/ASSR/code/Subset_SSCT_for_Ernie.xlsx');
% groupIdx = demographicData.Dx;
% 
% % Apply preselection for 13+ yrs olds.
% groupIdx = groupIdx(allAbove13Idx);

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

thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
gammaIdx = freqIdx7;

timeBins    = -1000:3999;
baselineIdx = find(timeBins>=-1000 & timeBins<=0);
onsetIdx    = find(timeBins>0 & timeBins<=500);
assrIdx     = find(timeBins>0 & timeBins<=3000);
offsetIdx   = find(timeBins>3000 & timeBins<=3500);
afterZeroIdx = find(timeBins>0);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Plot ERSSP on power. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');

stackData = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecErspMedian;
end
stackData      = 10*log10(stackData);
stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,1:1000,:),3));
erspFXS_M      = mean(stackData_base(:,:,:,fxsM_Idx),4);
erspFXS_F      = mean(stackData_base(:,:,:,fxsF_Idx),4);
erspHCS_M      = mean(stackData_base(:,:,:,hcsM_Idx),4);
erspHCS_F      = mean(stackData_base(:,:,:,hcsF_Idx),4);
% erspFXS_M      = median(stackData_base(:,:,:,fxsM_Idx),4);
% erspFXS_F      = median(stackData_base(:,:,:,fxsF_Idx),4);
% erspHCS_M      = median(stackData_base(:,:,:,hcsM_Idx),4);
% erspHCS_F      = median(stackData_base(:,:,:,hcsF_Idx),4);


%%%%%%%%%%%%%%
%%% FXS_M. %%%
%%%%%%%%%%%%%%
pTensor = zeros(size(erspFXS_M));
for chIdx = 1:size(erspFXS_M, 1)
    for freqIdxIdx = 1:size(erspFXS_M,2)
        currentFreqP = normcdf(squeeze(erspFXS_M(chIdx,freqIdxIdx,:)), squeeze(mean(erspFXS_M(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspFXS_M(chIdx,freqIdxIdx,baselineIdx),0,3)));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_FXS_M = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/pTensor_FXS_M', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')


%%%%%%%%%%%%%%
%%% FXS_F. %%%
%%%%%%%%%%%%%%
pTensor = zeros(size(erspFXS_F));
for chIdx = 1:size(erspFXS_F, 1)
    for freqIdxIdx = 1:size(erspFXS_F,2)
        currentFreqP = normcdf(squeeze(erspFXS_F(chIdx,freqIdxIdx,:)), squeeze(mean(erspFXS_F(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspFXS_F(chIdx,freqIdxIdx,baselineIdx),0,3)));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_FXS_F = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/pTensor_FXS_F', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')


%%%%%%%%%%%%%%
%%% HCS_M. %%%
%%%%%%%%%%%%%%
pTensor = zeros(size(erspHCS_M));
for chIdx = 1:size(erspHCS_M, 1)
    for freqIdxIdx = 1:size(erspHCS_M,2)
        currentFreqP = normcdf(squeeze(erspHCS_M(chIdx,freqIdxIdx,:)), squeeze(mean(erspHCS_M(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspHCS_M(chIdx,freqIdxIdx,baselineIdx),0,3)));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_HCS_M = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/pTensor_HCS_M', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')


%%%%%%%%%%%%%%
%%% HCS_F. %%%
%%%%%%%%%%%%%%
pTensor = zeros(size(erspHCS_F));
for chIdx = 1:size(erspHCS_F, 1)
    for freqIdxIdx = 1:size(erspHCS_F,2)
        currentFreqP = normcdf(erspHCS_F(chIdx,freqIdxIdx,:), squeeze(mean(erspHCS_F(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspHCS_F(chIdx,freqIdxIdx,baselineIdx),0,3)));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_HCS_F = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/pTensor_HCS_F', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot FXS-M FXS-F difference. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position', [150 150 2000 1000])
subplot(2,3,1)
imagesc(timeBins, 1:100, pTensorSum_FXS_M, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS Male grand mean (n=%d)', length(fxsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh')

subplot(2,3,2)
imagesc(timeBins, 1:100, pTensorSum_FXS_F, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS Female grand mean (n=%d)', length(fxsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');

subplot(2,3,3)
diffTensor = pTensorSum_FXS_M-pTensorSum_FXS_F;
imagesc(timeBins, 1:100, diffTensor, [-128 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title('FXS Male-Female')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot HCS-M HCS-F difference. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4)
imagesc(timeBins, 1:100, pTensorSum_HCS_M, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS Male grand mean (n=%d)', length(hcsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh')

subplot(2,3,5)
imagesc(timeBins, 1:100, pTensorSum_HCS_F, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS Female grand mean (n=%d)', length(hcsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');

subplot(2,3,6)
diffTensor = pTensorSum_HCS_M-pTensorSum_HCS_F;
imagesc(timeBins, 1:100, diffTensor, [-128 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title('HCS Male-Female')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');

print('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/ersspPower_M_F', '-dsvg')
print('/srv/Makoto/ASSR/p0701_visErssp_Ersp_maleFemale/ersspPower_M_F', '-djpeg95', '-r200')