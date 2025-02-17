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

thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
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


% Load Median ITC.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
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

% Prepare data for topoplot.
erspFXS_M_onsetResp  = squeeze(mean(mean(erspFXS_M(:,thetaIdx,onsetIdx,:),2),3));
erspFXS_M_assr       = squeeze(mean(mean(erspFXS_M(:,gammaIdx,assrIdx,:),2),3));
erspFXS_M_offsetResp = squeeze(mean(mean(erspFXS_M(:,thetaIdx,offsetIdx,:),2),3));
erspFXS_F_onsetResp  = squeeze(mean(mean(erspFXS_F(:,thetaIdx,onsetIdx,:),2),3));
erspFXS_F_assr       = squeeze(mean(mean(erspFXS_F(:,gammaIdx,assrIdx,:),2),3));
erspFXS_F_offsetResp = squeeze(mean(mean(erspFXS_F(:,thetaIdx,offsetIdx,:),2),3));
erspHCS_M_onsetResp  = squeeze(mean(mean(erspHCS_M(:,thetaIdx,onsetIdx,:),2),3));
erspHCS_M_assr       = squeeze(mean(mean(erspHCS_M(:,gammaIdx,assrIdx,:),2),3));
erspHCS_M_offsetResp = squeeze(mean(mean(erspHCS_M(:,thetaIdx,offsetIdx,:),2),3));
erspHCS_F_onsetResp  = squeeze(mean(mean(erspHCS_F(:,thetaIdx,onsetIdx,:),2),3));
erspHCS_F_assr       = squeeze(mean(mean(erspHCS_F(:,gammaIdx,assrIdx,:),2),3));
erspHCS_F_offsetResp = squeeze(mean(mean(erspHCS_F(:,thetaIdx,offsetIdx,:),2),3));
itcFXS_M_onsetResp   = squeeze(mean(mean(itcFXS_M(:,thetaIdx,onsetIdx,:),2),3));
itcFXS_M_assr        = squeeze(mean(mean(itcFXS_M(:,gammaIdx,assrIdx,:),2),3));
itcFXS_M_offsetResp  = squeeze(mean(mean(itcFXS_M(:,thetaIdx,offsetIdx,:),2),3));
itcFXS_F_onsetResp   = squeeze(mean(mean(itcFXS_F(:,thetaIdx,onsetIdx,:),2),3));
itcFXS_F_assr        = squeeze(mean(mean(itcFXS_F(:,gammaIdx,assrIdx,:),2),3));
itcFXS_F_offsetResp  = squeeze(mean(mean(itcFXS_F(:,thetaIdx,offsetIdx,:),2),3));
itcHCS_M_onsetResp   = squeeze(mean(mean(itcHCS_M(:,thetaIdx,onsetIdx,:),2),3));
itcHCS_M_assr        = squeeze(mean(mean(itcHCS_M(:,gammaIdx,assrIdx,:),2),3));
itcHCS_M_offsetResp  = squeeze(mean(mean(itcHCS_M(:,thetaIdx,offsetIdx,:),2),3));
itcHCS_F_onsetResp   = squeeze(mean(mean(itcHCS_F(:,thetaIdx,onsetIdx,:),2),3));
itcHCS_F_assr        = squeeze(mean(mean(itcHCS_F(:,gammaIdx,assrIdx,:),2),3));
itcHCS_F_offsetResp  = squeeze(mean(mean(itcHCS_F(:,thetaIdx,offsetIdx,:),2),3));

erspOnsetFxsCombined     = [mean(erspFXS_M_onsetResp,2); mean(erspFXS_F_onsetResp,2)];
erspOnsetFxsCombinedMax  = max(erspOnsetFxsCombined);
erspAssrFxsCombined      = [mean(erspFXS_M_assr,2); mean(erspFXS_F_assr,2)];
erspAssrFxsCombinedMax   = max(erspAssrFxsCombined);
erspOffsetFxsCombined    = [mean(erspFXS_M_offsetResp,2); mean(erspFXS_F_offsetResp,2)];
erspOffsetFxsCombinedMax = max(erspOffsetFxsCombined);
erspOnsetHcsCombined     = [mean(erspHCS_M_onsetResp,2); mean(erspHCS_F_onsetResp,2)];
erspOnsetHcsCombinedMax  = max(erspOnsetHcsCombined);
erspAssrHcsCombined      = [mean(erspHCS_M_assr,2); mean(erspHCS_F_assr,2)];
erspAssrHcsCombinedMax   = max(erspAssrHcsCombined);
erspOffsetHcsCombined    = [mean(erspHCS_M_offsetResp,2); mean(erspHCS_F_offsetResp,2)];
erspOffsetHcsCombinedMax = max(erspOffsetHcsCombined);
itcOnsetFxsCombined      = [mean(itcFXS_M_onsetResp,2); mean(itcFXS_F_onsetResp,2)];
itcOnsetFxsCombinedMax   = max(itcOnsetFxsCombined);
itcAssrFxsCombined       = [mean(itcFXS_M_assr,2); mean(itcFXS_F_assr,2)];
itcAssrFxsCombinedMax    = max(itcAssrFxsCombined);
itcOffsetFxsCombined     = [mean(itcFXS_M_offsetResp,2); mean(itcFXS_F_offsetResp,2)];
itcOffsetFxsCombinedMax  = max(itcOffsetFxsCombined);
itcOnsetHcsCombined      = [mean(itcHCS_M_onsetResp,2); mean(itcHCS_F_onsetResp,2)];
itcOnsetHcsCombinedMax   = max(itcOnsetHcsCombined);
itcAssrHcsCombined       = [mean(itcHCS_M_assr,2); mean(itcHCS_F_assr,2)];
itcAssrHcsCombinedMax    = max(itcAssrHcsCombined);
itcOffsetHcsCombined     = [mean(itcHCS_M_offsetResp,2); mean(itcHCS_F_offsetResp,2)];
itcOffsetHcsCombinedMax  = max(itcOffsetHcsCombined);

%%
figure

% ERSP Onset, FXS.
subplot(3,6,1)
handle1 = topoplot(mean(erspFXS_M_onsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOnsetFxsCombinedMax]);
title('Onset, FXS M')

subplot(3,6,6*1+1)
handle2 = topoplot(mean(erspFXS_F_onsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOnsetFxsCombinedMax]);
title('Onset, FXS F')

subplot(3,6,6*2+1)
[H,P,CI,STATS] = ttest2(erspFXS_M_onsetResp', erspFXS_F_onsetResp');
signicanceMap = P<0.05;
handle3 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('Onset, FXS M-F')

% ERSP ASSR, FXS
subplot(3,6,2)
handle4 = topoplot(mean(erspFXS_M_assr,2), EEG.chanlocs, 'maplimits', [0 erspAssrFxsCombinedMax]);
title('ASSR, FXS M')

subplot(3,6,6*1+2)
handle5 = topoplot(mean(erspFXS_F_assr,2), EEG.chanlocs, 'maplimits', [0 erspAssrFxsCombinedMax]);
title('ASSR, FXS F')

subplot(3,6,6*2+2)
[H,P,CI,STATS] = ttest2(erspFXS_M_assr', erspFXS_F_assr');
signicanceMap = P<0.05;
handle6 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('ASSR, FXS M-F')

% ERSP Offset, FXS.
subplot(3,6,3)
handle7 = topoplot(mean(erspFXS_M_offsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOffsetFxsCombinedMax]);
title('Offset, FXS M')

subplot(3,6,6*1+3)
handle8 = topoplot(mean(erspFXS_F_offsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOffsetFxsCombinedMax]);
title('Offset, FXS F')

subplot(3,6,6*2+3)
[H,P,CI,STATS] = ttest2(erspFXS_M_offsetResp', erspFXS_F_offsetResp');
signicanceMap = P<0.05;
handle9 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('Offset, FXS M-F')




% ERSP Onset, HCS.
subplot(3,6,4)
handle10 = topoplot(mean(erspHCS_M_onsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOnsetFxsCombinedMax]);
title('Onset, HCS M')

subplot(3,6,6*1+4)
handle11 = topoplot(mean(erspHCS_F_onsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOnsetFxsCombinedMax]);
title('Onset, HCS F')

subplot(3,6,6*2+4)
[H,P,CI,STATS] = ttest2(erspHCS_M_onsetResp', erspHCS_F_onsetResp');
signicanceMap = P<0.05;
handle12 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('Onset, HCS M-F')

% ERSP ASSR, HCS
subplot(3,6,5)
handle13 = topoplot(mean(erspHCS_M_assr,2), EEG.chanlocs, 'maplimits', [0 erspAssrFxsCombinedMax]);
title('ASSR, HCS M')

subplot(3,6,6*1+5)
handle14 = topoplot(mean(erspHCS_F_assr,2), EEG.chanlocs, 'maplimits', [0 erspAssrFxsCombinedMax]);
title('ASSR, HCS F')

subplot(3,6,6*2+5)
[H,P,CI,STATS] = ttest2(erspHCS_M_assr', erspHCS_F_assr');
signicanceMap = P<0.05;
handle15 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('ASSR, HCS M-F')

% ERSP Offset, HCS.
subplot(3,6,6)
handle16 = topoplot(mean(erspHCS_M_offsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOffsetFxsCombinedMax]);
title('Offset, HCS M')

subplot(3,6,6*1+6)
handle17 = topoplot(mean(erspHCS_F_offsetResp,2), EEG.chanlocs, 'maplimits', [0 erspOffsetFxsCombinedMax]);
title('Offset, HCS F')

subplot(3,6,6*2+6)
[H,P,CI,STATS] = ttest2(erspHCS_M_offsetResp', erspHCS_F_offsetResp');
signicanceMap = P<0.05;
handle18 = topoplot(signicanceMap.*STATS.tstat, EEG.chanlocs, 'maplimits', [max(abs(signicanceMap.*STATS.tstat))*-1 max(abs(signicanceMap.*STATS.tstat))]);
title('Offset, HCS M-F')

colormap(handle1.Parent, 'hot')
colormap(handle2.Parent, 'hot')
colormap(handle3.Parent, 'jet')
colormap(handle4.Parent, 'hot')
colormap(handle5.Parent, 'hot')
colormap(handle6.Parent, 'jet')
colormap(handle7.Parent, 'hot')
colormap(handle8.Parent, 'hot')
colormap(handle9.Parent, 'jet')
colormap(handle10.Parent, 'hot')
colormap(handle11.Parent, 'hot')
colormap(handle12.Parent, 'jet')
colormap(handle13.Parent, 'hot')
colormap(handle14.Parent, 'hot')
colormap(handle15.Parent, 'jet')
colormap(handle16.Parent, 'hot')
colormap(handle17.Parent, 'hot')
colormap(handle18.Parent, 'jet')

exportgraphics(gcf, '/srv/Makoto/ASSR/p1000_maleFemaleComparison_topos/ERSP_topos_3x6.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')


%% Test all-channel ERSP values.

[~,P1,~,STATS1] = ttest(mean(erspFXS_M_onsetResp,2),  mean(erspFXS_F_onsetResp,2));
[~,P2,~,STATS2] = ttest(mean(erspFXS_M_assr,2),       mean(erspFXS_F_assr,2));
[~,P3,~,STATS3] = ttest(mean(erspFXS_M_offsetResp,2), mean(erspFXS_F_offsetResp,2));
[~,P4,~,STATS4] = ttest(mean(erspHCS_M_onsetResp,2),  mean(erspHCS_F_onsetResp,2));
[~,P5,~,STATS5] = ttest(mean(erspHCS_M_assr,2),       mean(erspHCS_F_assr,2));
[~,P6,~,STATS6] = ttest(mean(erspHCS_M_offsetResp,2), mean(erspHCS_F_offsetResp,2));

[~,P1,~,STATS1] = ttest(mean(erspFXS_M_onsetResp,2),  mean(erspHCS_M_onsetResp,2));
[~,P2,~,STATS2] = ttest(mean(erspFXS_M_assr,2),       mean(erspHCS_M_assr,2));
[~,P3,~,STATS3] = ttest(mean(erspFXS_M_offsetResp,2), mean(erspHCS_M_offsetResp,2));
[~,P4,~,STATS4] = ttest(mean(erspFXS_F_onsetResp,2),  mean(erspHCS_F_onsetResp,2));
[~,P5,~,STATS5] = ttest(mean(erspFXS_F_assr,2),       mean(erspHCS_F_assr,2));
[~,P6,~,STATS6] = ttest(mean(erspFXS_F_offsetResp,2), mean(erspHCS_F_offsetResp,2));

[stats1, df1, pvals1] = statcond({mean(erspFXS_M_onsetResp,2)'  mean(erspFXS_F_onsetResp,2)';  mean(erspHCS_M_onsetResp,2)',  mean(erspHCS_F_onsetResp,2)'});
[stats2, df2, pvals2] = statcond({mean(erspFXS_M_assr,2)'       mean(erspFXS_F_assr,2)';       mean(erspHCS_M_assr,2)',       mean(erspHCS_F_assr,2)'});
[stats3, df3, pvals3] = statcond({mean(erspFXS_M_offsetResp,2)' mean(erspFXS_F_offsetResp,2)'; mean(erspHCS_M_offsetResp,2)', mean(erspHCS_F_offsetResp,2)'});

figure
subplot(1,3,1)
bar([mean(mean(erspFXS_M_onsetResp,2))  mean(mean(erspFXS_F_onsetResp,2))  mean(mean(erspHCS_M_onsetResp,2)),  mean(mean(erspHCS_F_onsetResp,2))])
title('Onset Response')
subplot(1,3,2)
bar([mean(mean(erspFXS_M_assr,2))  mean(mean(erspFXS_F_assr,2))  mean(mean(erspHCS_M_assr,2)),  mean(mean(erspHCS_F_assr,2))])
title('ASSR')
subplot(1,3,3)
bar([mean(mean(erspFXS_M_offsetResp,2))  mean(mean(erspFXS_F_offsetResp,2))  mean(mean(erspHCS_M_offsetResp,2)),  mean(mean(erspHCS_F_offsetResp,2))])
title('Offset Response')


%% Repeat ANOVA for 128 times.

figure

% Onset-evoked responses.
pvalsMatrix = ones(128,3);
meanMatrix  = zeros(128,4);
stdMatrix   = zeros(128,4);
for chIdx = 1:128
    FXS_M = erspFXS_M_onsetResp(chIdx,:);
    FXS_F = erspFXS_F_onsetResp(chIdx,:);
    HCS_M = erspHCS_M_onsetResp(chIdx,:);
    HCS_F = erspHCS_F_onsetResp(chIdx,:);
    [stats, df, pvals] = statcond({FXS_M FXS_F; HCS_M HCS_F});
    meanMatrix(chIdx,:)  = [mean(FXS_M) mean(FXS_F) mean(HCS_M) mean(HCS_F)];
    stdMatrix(chIdx,:)   = [std(FXS_M)  std(FXS_F)  std(HCS_M)  std(HCS_F)];
    pvalsMatrix(chIdx,:) = cell2mat(pvals);
end

subplot(2,3,1)
topoplot(pvalsMatrix(:,3)<0.05, EEG.chanlocs)
significantInteractionIdx = find(pvalsMatrix(:,3)<0.05);
sigIntValues = meanMatrix(significantInteractionIdx,:);
subplot(2,3,4)
bar(mean(sigIntValues))
title('ERSP, Onset-evoked resp.')

% ASSR.
for chIdx = 1:128
    FXS_M = erspFXS_M_assr(chIdx,:);
    FXS_F = erspFXS_F_assr(chIdx,:);
    HCS_M = erspHCS_M_assr(chIdx,:);
    HCS_F = erspHCS_F_assr(chIdx,:);
    [stats, df, pvals] = statcond({FXS_M FXS_F; HCS_M HCS_F});
    meanMatrix(chIdx,:)  = [mean(FXS_M) mean(FXS_F) mean(HCS_M) mean(HCS_F)];
    stdMatrix(chIdx,:)   = [std(FXS_M)  std(FXS_F)  std(HCS_M)  std(HCS_F)];
    pvalsMatrix(chIdx,:) = cell2mat(pvals);
end

subplot(2,3,2)
topoplot(pvalsMatrix(:,3)<0.05, EEG.chanlocs)
significantInteractionIdx = find(pvalsMatrix(:,3)<0.05);
sigIntValues = meanMatrix(significantInteractionIdx,:);
subplot(2,3,5)
bar(mean(sigIntValues))
title('ERSP, ASSR.')

% Offset-evoked responses.
pvalsMatrix = ones(128,3);
meanMatrix  = zeros(128,4);
stdMatrix   = zeros(128,4);
for chIdx = 1:128
    FXS_M = erspFXS_M_offsetResp(chIdx,:);
    FXS_F = erspFXS_F_offsetResp(chIdx,:);
    HCS_M = erspHCS_M_offsetResp(chIdx,:);
    HCS_F = erspHCS_F_offsetResp(chIdx,:);
    [stats, df, pvals] = statcond({FXS_M FXS_F; HCS_M HCS_F});
    meanMatrix(chIdx,:)  = [mean(FXS_M) mean(FXS_F) mean(HCS_M) mean(HCS_F)];
    stdMatrix(chIdx,:)   = [std(FXS_M)  std(FXS_F)  std(HCS_M)  std(HCS_F)];
    pvalsMatrix(chIdx,:) = cell2mat(pvals);
end

subplot(2,3,3)
topoplot(pvalsMatrix(:,3)<0.05, EEG.chanlocs)
significantInteractionIdx = find(pvalsMatrix(:,3)<0.05);
sigIntValues = meanMatrix(significantInteractionIdx,:);
subplot(2,3,6)
bar(mean(sigIntValues))
title('ERSP, Offset-evoked resp.')

exportgraphics(gcf, '/srv/Makoto/ASSR/p1000_maleFemaleComparison_topos/ERSP_1x3_allChanANOVA.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')

