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

%thetaIdx = find(wtFreqBins>=2 & wtFreqBins<=4);
thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
gammaIdx = freqIdx7;

timeBins    = -1000:4000;
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
erspFXS_M      = stackData_base(:,:,:,fxsM_Idx);
erspFXS_F      = stackData_base(:,:,:,fxsF_Idx);
erspHCS_M      = stackData_base(:,:,:,hcsM_Idx);
erspHCS_F      = stackData_base(:,:,:,hcsF_Idx);


%%%%%%%%%%%%%%
%%% FXS_M. %%%
%%%%%%%%%%%%%%
pTensor = zeros(size(elecErspMedian));
ERSSP_FXS_M = single(zeros(100, 5001, length(fxsM_Idx)));
for subjIdx = 1:length(fxsM_Idx)
    currentSubjData = erspFXS_M(:,:,:,subjIdx);
    for chIdx = 1:size(currentSubjData, 1)
        for freqIdx = 1:size(currentSubjData,2)
            currentFreqP = normcdf(squeeze(currentSubjData(chIdx,freqIdx,:)), squeeze(mean(currentSubjData(chIdx,freqIdx,baselineIdx),3)), squeeze(std(currentSubjData(chIdx,freqIdx,baselineIdx),0,3)));
            pTensor(chIdx,freqIdx,:) = currentFreqP;
        end
    end

    % Determine functional ROI.
    pTensor_rTale = 1-pTensor;
    pTensor_sigMask = pTensor_rTale<0.05;
    ERSSP = squeeze(sum(pTensor_sigMask));
    ERSSP_FXS_M(:,:,subjIdx) = ERSSP;
end


%%%%%%%%%%%%%%
%%% FXS_F. %%%
%%%%%%%%%%%%%%
ERSSP_FXS_F = single(zeros(100, 5001, length(fxsF_Idx)));
for subjIdx = 1:length(fxsF_Idx)
    currentSubjData = erspFXS_F(:,:,:,subjIdx);
    for chIdx = 1:size(currentSubjData, 1)
        for freqIdx = 1:size(currentSubjData,2)
            currentFreqP = normcdf(squeeze(currentSubjData(chIdx,freqIdx,:)), squeeze(mean(currentSubjData(chIdx,freqIdx,baselineIdx),3)), squeeze(std(currentSubjData(chIdx,freqIdx,baselineIdx),0,3)));
            pTensor(chIdx,freqIdx,:) = currentFreqP;
        end
    end

    % Determine functional ROI.
    pTensor_rTale = 1-pTensor;
    pTensor_sigMask = pTensor_rTale<0.05;
    ERSSP = squeeze(sum(pTensor_sigMask));
    ERSSP_FXS_F(:,:,subjIdx) = ERSSP;
end


%%%%%%%%%%%%%%
%%% HCS_M. %%%
%%%%%%%%%%%%%%
ERSSP_HCS_M = single(zeros(100, 5001, length(hcsM_Idx)));
for subjIdx = 1:length(hcsM_Idx)
    currentSubjData = erspHCS_M(:,:,:,subjIdx);
    for chIdx = 1:size(currentSubjData, 1)
        for freqIdx = 1:size(currentSubjData,2)
            currentFreqP = normcdf(squeeze(currentSubjData(chIdx,freqIdx,:)), squeeze(mean(currentSubjData(chIdx,freqIdx,baselineIdx),3)), squeeze(std(currentSubjData(chIdx,freqIdx,baselineIdx),0,3)));
            pTensor(chIdx,freqIdx,:) = currentFreqP;
        end
    end

    % Determine functional ROI.
    pTensor_rTale = 1-pTensor;
    pTensor_sigMask = pTensor_rTale<0.05;
    ERSSP = squeeze(sum(pTensor_sigMask));
    ERSSP_HCS_M(:,:,subjIdx) = ERSSP;
end


%%%%%%%%%%%%%%
%%% HCS_F. %%%
%%%%%%%%%%%%%%
ERSSP_HCS_F = single(zeros(100, 5001, length(hcsF_Idx)));
for subjIdx = 1:length(hcsF_Idx)
    currentSubjData = erspHCS_F(:,:,:,subjIdx);
    for chIdx = 1:size(currentSubjData, 1)
        for freqIdx = 1:size(currentSubjData,2)
            currentFreqP = normcdf(squeeze(currentSubjData(chIdx,freqIdx,:)), squeeze(mean(currentSubjData(chIdx,freqIdx,baselineIdx),3)), squeeze(std(currentSubjData(chIdx,freqIdx,baselineIdx),0,3)));
            pTensor(chIdx,freqIdx,:) = currentFreqP;
        end
    end

    % Determine functional ROI.
    pTensor_rTale = 1-pTensor;
    pTensor_sigMask = pTensor_rTale<0.05;
    ERSSP = squeeze(sum(pTensor_sigMask));
    ERSSP_HCS_F(:,:,subjIdx) = ERSSP;
end

save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_FXS_M', 'ERSSP_FXS_M')
save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_FXS_F', 'ERSSP_FXS_F')
save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_HCS_M', 'ERSSP_HCS_M')
save('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ERSSP_HCS_F', 'ERSSP_HCS_F')


%% Visualize the results.
ERSSP_FXS_M(ERSSP_FXS_M<=6) = 0;
ERSSP_FXS_F(ERSSP_FXS_F<=6) = 0;
ERSSP_HCS_M(ERSSP_HCS_M<=6) = 0;
ERSSP_HCS_F(ERSSP_HCS_F<=6) = 0;

    %{
    trueDiff = (mean(onset_FXS_M)-mean(onset_FXS_F))-(mean(onset_HCS_M)-mean(onset_HCS_F));
    concatenated = [onset_FXS_M; onset_FXS_F; onset_HCS_M; onset_HCS_F];
    permute_concatenated = randi(length(concatenated), 68, 8000);
    surroData = concatenated(permute_concatenated);
    surro_FXS_M = surroData(1:23, 1   :2000);
    surro_FXS_F = surroData(24:35,2001:4000);
    surro_HCS_M = surroData(36:57,4001:6000);
    surro_HCS_F = surroData(58:68,6001:8000);
    surroDiff = (mean(surro_FXS_M)-mean(surro_FXS_F))-(mean(surro_HCS_M)-mean(surro_HCS_F));
    figure; hist(surroDiff, 50)
    %}


%% Massive ANOVA.
resultMatrix = zeros(size(ERSSP_FXS_M, 1), size(ERSSP_FXS_M, 2));
for freqIdx = 1:size(ERSSP_FXS_M, 1)
    for timeIdx = 1:size(ERSSP_FXS_M, 2)
        [stats, df, pvals] = statcond({squeeze(ERSSP_FXS_M(thetaIdx(freqIdx), onsetIdx(timeIdx), :))' ...
                                        squeeze(ERSSP_FXS_F(thetaIdx(freqIdx), onsetIdx(timeIdx), :))'; ...
                                        squeeze(ERSSP_HCS_M(thetaIdx(freqIdx), onsetIdx(timeIdx), :))' ...
                                        squeeze(ERSSP_HCS_F(thetaIdx(freqIdx), onsetIdx(timeIdx), :))});
        resultMatrix(freqIdx,timeIdx) = pvals{3};
    end
end

figure
imagesc(resultMatrix)

significanceMask = resultMatrix<0.05;

figure
imagesc(significanceMask)

onset_FXS_M = squeeze(mean(mean(ERSSP_FXS_M(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
onset_FXS_F = squeeze(mean(mean(ERSSP_FXS_F(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
onset_HCS_M = squeeze(mean(mean(ERSSP_HCS_M(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
onset_HCS_F = squeeze(mean(mean(ERSSP_HCS_F(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';

figure
bar([mean(onset_FXS_M) mean(onset_FXS_F) mean(onset_HCS_M) mean(onset_HCS_F)]);
    % resultMatrix = zeros(length(thetaIdx), length(onsetIdx));
    % for freqIdx = 1:length(thetaIdx)
    %     for timeIdx = 1:length(onsetIdx)
    % 
    %         [stats, df, pvals] = statcond({squeeze(ERSSP_FXS_M(thetaIdx(freqIdx), onsetIdx(timeIdx), :))' ...
    %             squeeze(ERSSP_FXS_F(thetaIdx(freqIdx), onsetIdx(timeIdx), :))'; ...
    %             squeeze(ERSSP_HCS_M(thetaIdx(freqIdx), onsetIdx(timeIdx), :))' ...
    %             squeeze(ERSSP_HCS_F(thetaIdx(freqIdx), onsetIdx(timeIdx), :))});
    %         resultMatrix(freqIdx,timeIdx) = pvals{3};
    %     end
    % end
    % 
    % figure
    % imagesc(resultMatrix)
    % 
    % significanceMask = resultMatrix<0.05;
    % 
    % figure
    % imagesc(significanceMask)
    % 
    % onset_FXS_M = squeeze(mean(mean(ERSSP_FXS_M(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
    % onset_FXS_F = squeeze(mean(mean(ERSSP_FXS_F(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
    % onset_HCS_M = squeeze(mean(mean(ERSSP_HCS_M(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
    % onset_HCS_F = squeeze(mean(mean(ERSSP_HCS_F(thetaIdx(1:10), onsetIdx(1:170), :),1),2))';
    % 
    % figure
    % bar([mean(onset_FXS_M) mean(onset_FXS_F) mean(onset_HCS_M) mean(onset_HCS_F)]);





onset_FXS_M = squeeze(sum(sum(ERSSP_FXS_M(thetaIdx, onsetIdx, :),1),2)/(length(thetaIdx)*length(onsetIdx)))';
onset_FXS_F = squeeze(sum(sum(ERSSP_FXS_F(thetaIdx, onsetIdx, :),1),2)/(length(thetaIdx)*length(onsetIdx)))';
onset_HCS_M = squeeze(sum(sum(ERSSP_HCS_M(thetaIdx, onsetIdx, :),1),2)/(length(thetaIdx)*length(onsetIdx)))';
onset_HCS_F = squeeze(sum(sum(ERSSP_HCS_F(thetaIdx, onsetIdx, :),1),2)/(length(thetaIdx)*length(onsetIdx)))';
    % onset_FXS_M = median(reshape(ERSSP_FXS_M(thetaIdx, onsetIdx, :), [length(thetaIdx)*length(onsetIdx) size(ERSSP_FXS_M,3)]),1);
    % onset_FXS_F = median(reshape(ERSSP_FXS_F(thetaIdx, onsetIdx, :), [length(thetaIdx)*length(onsetIdx) size(ERSSP_FXS_F,3)]),1);
    % onset_HCS_M = median(reshape(ERSSP_HCS_M(thetaIdx, onsetIdx, :), [length(thetaIdx)*length(onsetIdx) size(ERSSP_HCS_M,3)]),1);
    % onset_HCS_F = median(reshape(ERSSP_HCS_F(thetaIdx, onsetIdx, :), [length(thetaIdx)*length(onsetIdx) size(ERSSP_HCS_F,3)]),1);

[stats1, df1, pvals1] = statcond({onset_FXS_M onset_FXS_F; onset_HCS_M onset_HCS_F});

assr_FXS_M = squeeze(sum(ERSSP_FXS_M(gammaIdx, assrIdx, :),2)/length(assrIdx));
assr_FXS_F = squeeze(sum(ERSSP_FXS_F(gammaIdx, assrIdx, :),2)/length(assrIdx));
assr_HCS_M = squeeze(sum(ERSSP_HCS_M(gammaIdx, assrIdx, :),2)/length(assrIdx));
assr_HCS_F = squeeze(sum(ERSSP_HCS_F(gammaIdx, assrIdx, :),2)/length(assrIdx));
[stats2, df2, pvals2] = statcond({assr_FXS_M' assr_FXS_F'; assr_HCS_M' assr_HCS_F'});

offset_FXS_M = squeeze(sum(sum(ERSSP_FXS_M(thetaIdx, offsetIdx, :),1),2)/(length(thetaIdx)*length(offsetIdx)));
offset_FXS_F = squeeze(sum(sum(ERSSP_FXS_F(thetaIdx, offsetIdx, :),1),2)/(length(thetaIdx)*length(offsetIdx)));
offset_HCS_M = squeeze(sum(sum(ERSSP_HCS_M(thetaIdx, offsetIdx, :),1),2)/(length(thetaIdx)*length(offsetIdx)));
offset_HCS_F = squeeze(sum(sum(ERSSP_HCS_F(thetaIdx, offsetIdx, :),1),2)/(length(thetaIdx)*length(offsetIdx)));
[stats3, df3, pvals3] = statcond({offset_FXS_M' offset_FXS_F'; offset_HCS_M' offset_HCS_F'});

% figure
% subplot(1,3,1)
% bar([mean(onset_FXS_M)  mean(onset_FXS_F)  mean(onset_HCS_M)  mean(onset_HCS_F)])
% title('Onset Response')
% subplot(1,3,2)
% bar([mean(assr_FXS_M)   mean(assr_FXS_F)   mean(assr_HCS_M)   mean(assr_HCS_F)])
% title('ASSR')
% subplot(1,3,3)
% bar([mean(offset_FXS_M) mean(offset_FXS_F) mean(offset_HCS_M) mean(offset_HCS_F)])
% title('Offset Response')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot FXS-M FXS-F difference. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position', [150 150 2000 1000])
subplot(2,3,1)
imagesc(timeBins, 1:100, mean(ERSSP_FXS_M,3), [0 55])
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
imagesc(timeBins, 1:100, mean(ERSSP_FXS_F,3), [0 55])
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
diffTensor = mean(ERSSP_FXS_M,3)-mean(ERSSP_FXS_F,3);
imagesc(timeBins, 1:100, diffTensor, [-25 25])
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
imagesc(timeBins, 1:100, mean(ERSSP_HCS_M,3), [0 55])
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
imagesc(timeBins, 1:100, mean(ERSSP_HCS_F,3), [0 55])
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
diffTensor = mean(ERSSP_HCS_M,3)-mean(ERSSP_HCS_F,3);
imagesc(timeBins, 1:100, diffTensor, [-25 25])
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

% print('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ersspPower_M_F', '-dsvg')
% print('/srv/Makoto/ASSR/p1010_maleFemale_singleSubjERSSP/ersspPower_M_F', '-djpeg95', '-r200')


%%
meanERSSP = mean(ERSSP_FXS_M,3);
pMatrix   = zeros(size(meanERSSP));
for freqIdx = 1:size(meanERSSP,1)
    currentFreqP = normcdf(squeeze(meanERSSP(freqIdx,:)), mean(vec(meanERSSP(:,baselineIdx))), std(vec(meanERSSP(:,baselineIdx))));
    pMatrix(freqIdx,:) = currentFreqP';
end

% Determine functional ROI.
pMatrix_rTale = 1-pMatrix;
pMatrix_sigMask = pMatrix_rTale<0.0000000000000001;
figure
imagesc(pMatrix_sigMask)