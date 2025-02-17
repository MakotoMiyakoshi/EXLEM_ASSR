% 01/08/2025 Makoto. Updated.
% 12/11/2024 Makoto. Created.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools')
addpath('/srv/Makoto/Tools/frank-pk-DataViz-3.2.3.0/daboxplot')

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


%%
% Load the latest demographic data.
demoTable = readtable('/srv/Makoto/ASSR/p0000_raw/Final_List_SSCT_Adults_Adolescents_01062025.xlsx');
eegIdVec = demoTable.eegid;

% Exclusion list by Craig.
exclusionList = [44 532 554 879 2537];
[~, preselectingIdx1] = setdiff(eegIdVec, exclusionList);
demoTable = demoTable(preselectingIdx1, :);

% Obtain valid subject IDs.
validEegId   = demoTable.eegid;
subjNamesVec = cellfun(@str2num, subjNames');
[~, preselectingIdx2] = intersect(subjNamesVec, validEegId);

% Apply the trimming.
subjNames  = subjNames(preselectingIdx2);
groupNames = groupNames(preselectingIdx2);
sexGroupList = sexGroupList(preselectingIdx2,:);

% Validate the subject order.
if sum(subjNamesVec(preselectingIdx2)-validEegId) ~= 0
    error('Two lists do not match.')
end

% Obtain FMRP.
fmrpVec = demoTable.fmrp_replicate_trim_mean;


%%
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
freqIdx = [freqIdx1 freqIdx2 freqIdx3 freqIdx4 freqIdx5 freqIdx6 freqIdx7 freqIdx8];
freqLabels = [1 2 4 8 13 20 40 80];

deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
gammaIdx = freqIdx7;

timeBins    = -1000:4000;
baselineIdx = find(timeBins>=-1000 & timeBins<=0);
onsetIdx    = find(timeBins>0 & timeBins<=500);
assrIdx     = find(timeBins>0 & timeBins<=3000);
offsetIdx   = find(timeBins>3000 & timeBins<=3500);
afterZeroIdx = find(timeBins>0);

%% Cz Median itc.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
allMats = allMats(preselectingIdx2);
stackData = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMedian;
end
itcFXS_M      = stackData(:,:,:,fxsM_Idx);
itcFXS_F      = stackData(:,:,:,fxsF_Idx);
itcHCS_M      = stackData(:,:,:,hcsM_Idx);
itcHCS_F      = stackData(:,:,:,hcsF_Idx);
itcFXS        = stackData(:,:,:,fxsIdx);
itcHCS        = stackData(:,:,:,hcsIdx);

fxs_demographics = [subjNames(fxsIdx)' sexGroupList(fxsIdx,:)];
hcs_demographics = [subjNames(hcsIdx)' sexGroupList(hcsIdx,:)];

fxsmGm_FCz = squeeze(mean(itcFXS_M(6,:,:,:),4));
fxsfGm_FCz = squeeze(mean(itcFXS_F(6,:,:,:),4));
hcsmGm_FCz = squeeze(mean(itcHCS_M(6,:,:,:),4));
hcsfGm_FCz = squeeze(mean(itcHCS_F(6,:,:,:),4));
fxsGm_FCz  = squeeze(mean(itcFXS(6,:,:,:),4));
hcsGm_FCz  = squeeze(mean(itcHCS(6,:,:,:),4));

fxsmGm_Cz = squeeze(mean(itcFXS_M(55,:,:,:),4));
fxsfGm_Cz = squeeze(mean(itcFXS_F(55,:,:,:),4));
hcsmGm_Cz = squeeze(mean(itcHCS_M(55,:,:,:),4));
hcsfGm_Cz = squeeze(mean(itcHCS_F(55,:,:,:),4));
fxsGm_Cz  = squeeze(mean(itcFXS(55,:,:,:),4));
hcsGm_Cz  = squeeze(mean(itcHCS(55,:,:,:),4));


% Scalp topos.
onsetFXS_topo = squeeze(mean(mean(mean(itcFXS(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));
onsetHCS_topo = squeeze(mean(mean(mean(itcHCS(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));

assrFXS_topo = squeeze(mean(mean(mean(itcFXS(:, gammaIdx, assrIdx, :),2),3),4));
assrHCS_topo = squeeze(mean(mean(mean(itcHCS(:, gammaIdx, assrIdx, :),2),3),4));

offsetFXS_topo = squeeze(mean(mean(mean(itcFXS(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));
offsetHCS_topo = squeeze(mean(mean(mean(itcHCS(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));


% ROI mean.
onsetFXS_roiMean = squeeze(mean(mean(itcFXS(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));
onsetHCS_roiMean = squeeze(mean(mean(itcHCS(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));

assrFXS_roiMean = squeeze(mean(mean(itcFXS(6, gammaIdx, assrIdx, :),2),3));
assrHCS_roiMean = squeeze(mean(mean(itcHCS(6, gammaIdx, assrIdx, :),2),3));

offsetFXS_roiMean = squeeze(mean(mean(itcFXS(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));
offsetHCS_roiMean = squeeze(mean(mean(itcHCS(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));


%%
close all

figure('position', [50 50 1000 1200])
mainLayout = tiledlayout(3,1,"TileSpacing", "tight", "Padding", "tight");

% itc.
nestedLayout1 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 1;
nexttile(nestedLayout1, 1)
imagesc(timeBins, 1:100, fxsGm_Cz, [0 0.6])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
axis xy
title(sprintf('FXS (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
ylabel('Frequency (Hz)')

nexttile(nestedLayout1, 2)
imagesc(timeBins, 1:100, hcsGm_Cz, [0 0.6])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels, 'box', 'off')
axis xy
title(sprintf('HCS (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
ylabel('Frequency (Hz)')


% Topos.
nestedLayout2 = tiledlayout(mainLayout, 1, 6, 'TileSpacing', 'none', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 2;
nexttile(nestedLayout2, 1)
topoplot(onsetFXS_topo, EEG.chanlocs, 'maplimits', [0 0.4]);
nexttile(nestedLayout2, 2)
topoplot(onsetHCS_topo, EEG.chanlocs, 'maplimits', [0 0.4]);
nexttile(nestedLayout2, 3)
topoplot(assrFXS_topo, EEG.chanlocs, 'maplimits', [0 0.6]);
nexttile(nestedLayout2, 4)
topoplot(assrHCS_topo, EEG.chanlocs, 'maplimits', [0 0.6]);
nexttile(nestedLayout2, 5)
topoplot(offsetFXS_topo, EEG.chanlocs, 'maplimits', [0 0.4]);
nexttile(nestedLayout2, 6)
topoplot(offsetHCS_topo, EEG.chanlocs, 'maplimits', [0 0.4]);
colormap jet


% Bar graphs.
nestedLayout3 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout3.Layout.Tile = 3;

nexttile(nestedLayout3, 1)

h = daboxplot({onsetFXS_roiMean, onsetHCS_roiMean},'linkline', 1, 'xtlabels', {'FXS' 'HCS'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [1 0 0];
h.bx(2).FaceColor = [0 0 1];
h.bx(1).FaceAlpha = 0.2;
h.bx(2).FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
title('OnsetER at Cz')

[H1,P1,CI1,STATS1] = ttest2(onsetFXS_roiMean, onsetHCS_roiMean);
d1 = cohensD(onsetFXS_roiMean, onsetHCS_roiMean);
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS1.df, STATS1.tstat, P1, d1))
    % t(60)=0.874, p=0.385, Cohen's d = 0.223.


nexttile(nestedLayout3, 2)

h = daboxplot({assrFXS_roiMean, assrHCS_roiMean},'linkline', 1, 'xtlabels', {'FXS' 'HCS'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [1 0 0];
h.bx(2).FaceColor = [0 0 1];
h.bx(1).FaceAlpha = 0.2;
h.bx(2).FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
title('ASSR at FCz')

[H2,P2,CI2,STATS2] = ttest2(assrFXS_roiMean, assrHCS_roiMean);
d2 = cohensD(assrFXS_roiMean, assrHCS_roiMean);
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS2.df, STATS2.tstat, P2, d2))
    % t(60)=-1.372, p=0.175, Cohen's d = -0.349.


nexttile(nestedLayout3, 3)

h = daboxplot({offsetFXS_roiMean, offsetHCS_roiMean},'linkline', 1, 'xtlabels', {'FXS' 'HCS'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [1 0 0];
h.bx(2).FaceColor = [0 0 1];
h.bx(1).FaceAlpha = 0.2;
h.bx(2).FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
title('OffsetER at Cz')

[H3,P3,CI3,STATS3] = ttest2(offsetFXS_roiMean, offsetHCS_roiMean);
d3 = cohensD(offsetFXS_roiMean, offsetHCS_roiMean);
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS3.df, STATS3.tstat, P3, d3))
    % t(60)=1.572, p=0.121, Cohen's d = 0.400.

% Export the results.
combinedDemographics = [fxs_demographics; hcs_demographics];
fmrp     = [fmrpVec(fxsIdx); fmrpVec(hcsIdx)];
subjID   = combinedDemographics(:,1);
sex      = combinedDemographics(:,2);
group    = combinedDemographics(:,3);
onsetER  = [onsetFXS_roiMean; onsetHCS_roiMean];
ASSR     = [assrFXS_roiMean; assrHCS_roiMean];
offsetER = [offsetFXS_roiMean; offsetHCS_roiMean];
resultTableItc = table(subjID, sex, group, onsetER, ASSR, offsetER);

exportgraphics(gcf, '/srv/Makoto/ASSR/p1312_itc/itcs.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')
save('/srv/Makoto/ASSR/p1312_itc/resultTableItc.mat', 'resultTableItc')