% 12/11/2024 Makoto. Created.
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

%% Cz Median ERSP.
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

fxsmGm_FCz = squeeze(mean(erspFXS_M(6,:,:,:),4));
fxsfGm_FCz = squeeze(mean(erspFXS_F(6,:,:,:),4));
hcsmGm_FCz = squeeze(mean(erspHCS_M(6,:,:,:),4));
hcsfGm_FCz = squeeze(mean(erspHCS_F(6,:,:,:),4));

fxsmGm_Cz = squeeze(mean(erspFXS_M(55,:,:,:),4));
fxsfGm_Cz = squeeze(mean(erspFXS_F(55,:,:,:),4));
hcsmGm_Cz = squeeze(mean(erspHCS_M(55,:,:,:),4));
hcsfGm_Cz = squeeze(mean(erspHCS_F(55,:,:,:),4));



% Scalp topos.
onsetFXS_M_topo = squeeze(mean(mean(mean(erspFXS_M(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));
onsetFXS_F_topo = squeeze(mean(mean(mean(erspFXS_F(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));
onsetHCS_M_topo = squeeze(mean(mean(mean(erspHCS_M(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));
onsetHCS_F_topo = squeeze(mean(mean(mean(erspHCS_F(:, deltaThetaAlphaIdx, onsetIdx, :),2),3),4));
assrFXS_M_topo = squeeze(mean(mean(mean(erspFXS_M(:, gammaIdx, assrIdx, :),2),3),4));
assrFXS_F_topo = squeeze(mean(mean(mean(erspFXS_F(:, gammaIdx, assrIdx, :),2),3),4));
assrHCS_M_topo = squeeze(mean(mean(mean(erspHCS_M(:, gammaIdx, assrIdx, :),2),3),4));
assrHCS_F_topo = squeeze(mean(mean(mean(erspHCS_F(:, gammaIdx, assrIdx, :),2),3),4));
offsetFXS_M_topo = squeeze(mean(mean(mean(erspFXS_M(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));
offsetFXS_F_topo = squeeze(mean(mean(mean(erspFXS_F(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));
offsetHCS_M_topo = squeeze(mean(mean(mean(erspHCS_M(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));
offsetHCS_F_topo = squeeze(mean(mean(mean(erspHCS_F(:, deltaThetaAlphaIdx, offsetIdx, :),2),3),4));

% ROI mean.
onsetFXS_M_roiMean = squeeze(mean(mean(erspFXS_M(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));
onsetFXS_F_roiMean = squeeze(mean(mean(erspFXS_F(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));
onsetHCS_M_roiMean = squeeze(mean(mean(erspHCS_M(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));
onsetHCS_F_roiMean = squeeze(mean(mean(erspHCS_F(55, deltaThetaAlphaIdx, onsetIdx, :),2),3));
assrFXS_M_roiMean = squeeze(mean(mean(erspFXS_M(6, gammaIdx, assrIdx, :),2),3));
assrFXS_F_roiMean = squeeze(mean(mean(erspFXS_F(6, gammaIdx, assrIdx, :),2),3));
assrHCS_M_roiMean = squeeze(mean(mean(erspHCS_M(6, gammaIdx, assrIdx, :),2),3));
assrHCS_F_roiMean = squeeze(mean(mean(erspHCS_F(6, gammaIdx, assrIdx, :),2),3));
offsetFXS_M_roiMean = squeeze(mean(mean(erspFXS_M(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));
offsetFXS_F_roiMean = squeeze(mean(mean(erspFXS_F(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));
offsetHCS_M_roiMean = squeeze(mean(mean(erspHCS_M(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));
offsetHCS_F_roiMean = squeeze(mean(mean(erspHCS_F(55, deltaThetaAlphaIdx, offsetIdx, :),2),3));


%%
close all

figure('position', [50 50 1000 1200])
mainLayout = tiledlayout(4,1,"TileSpacing", "tight", "Padding", "tight");

% ERSP.
nestedLayout1 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 1;
nexttile(nestedLayout1, 1)
imagesc(timeBins, 1:100, fxsmGm_Cz, [-1.5 1.5])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
axis xy
title(sprintf('FXS Male (n=%d)', length(fxsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
ylabel('Frequency (Hz)')

nexttile(nestedLayout1, 2)
imagesc(timeBins, 1:100, fxsfGm_Cz, [-1.5 1.5])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
axis xy
title(sprintf('FXS Female (n=%d)', length(fxsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
ylabel('Frequency (Hz)')

nestedLayout2 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 2;
nexttile(nestedLayout2, 1)
imagesc(timeBins, 1:100, hcsmGm_Cz, [-1.5 1.5])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
axis xy
title(sprintf('HCS Male (n=%d)', length(hcsM_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

nexttile(nestedLayout2, 2)
imagesc(timeBins, 1:100, hcsfGm_Cz, [-1.5 1.5])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
axis xy
title(sprintf('HCS Female (n=%d)', length(hcsF_Idx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 1)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 1, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')

% Topos.
nestedLayout3 = tiledlayout(mainLayout, 2, 6, 'TileSpacing', 'none', 'Padding', 'tight');
nestedLayout3.Layout.Tile = 3;
nexttile(nestedLayout3, 1)
topoplot(onsetFXS_M_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 2)
topoplot(onsetFXS_F_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 3)
topoplot(assrFXS_M_topo, EEG.chanlocs, 'maplimits', [-1.3 1.3]);
nexttile(nestedLayout3, 4)
topoplot(assrFXS_F_topo, EEG.chanlocs, 'maplimits', [-1.3 1.3]);
nexttile(nestedLayout3, 5)
topoplot(offsetFXS_M_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 6)
topoplot(offsetFXS_F_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 7)
topoplot(onsetHCS_M_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 8)
topoplot(onsetHCS_F_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 9)
topoplot(assrHCS_M_topo, EEG.chanlocs, 'maplimits', [-1.3 1.3]);
nexttile(nestedLayout3, 10)
topoplot(assrHCS_F_topo, EEG.chanlocs, 'maplimits', [-1.3 1.3]);
nexttile(nestedLayout3, 11)
topoplot(offsetHCS_M_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);
nexttile(nestedLayout3, 12)
topoplot(offsetHCS_F_topo, EEG.chanlocs, 'maplimits', [-0.8 0.8]);

colormap jet

% Bar graphs
nestedLayout4 = tiledlayout(mainLayout, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout4.Layout.Tile = 4;
nexttile(nestedLayout4, 1)
barHandle = bar([mean(onsetFXS_M_roiMean) mean(onsetFXS_F_roiMean) mean(onsetHCS_M_roiMean) mean(onsetHCS_F_roiMean)])
hold on
errorbar(1:4, [mean(onsetFXS_M_roiMean) mean(onsetFXS_F_roiMean) mean(onsetHCS_M_roiMean) mean(onsetHCS_F_roiMean)], ...
              [0 0 0 0], ...
              [ std(onsetFXS_M_roiMean)  std(onsetFXS_F_roiMean)  std(onsetHCS_M_roiMean)  std(onsetHCS_F_roiMean)], ...
              'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('OnsetER at Cz')
ylim([0 1.5])
[tStats_onset, df_onset, pvals_onset] = statcond({onsetFXS_M_roiMean' onsetFXS_F_roiMean'; onsetHCS_M_roiMean' onsetHCS_F_roiMean'});
    % t: {[10.5023]}    {[10.3215]}    {[2.5985]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p: {[0.0017]}    {[0.0018]}    {[0.1105]}

[onsetD_MaleFemale, onsetD_FxsHcs] = cohensD_mainEffects2x2(onsetFXS_M_roiMean', onsetFXS_F_roiMean', onsetHCS_M_roiMean', onsetHCS_F_roiMean');
onsetD_F2                          = cohensF2(tStats_onset{3}, df_onset{3}(1), df_onset{3}(2));
disp(sprintf('%.3f %.3f %.3f',onsetD_MaleFemale, onsetD_FxsHcs, onsetD_F2))
    % Cohen 0.580 0.659 0.030

nexttile(nestedLayout4, 2)
barHandle = bar([mean(assrFXS_M_roiMean) mean(assrFXS_F_roiMean) mean(assrHCS_M_roiMean) mean(assrHCS_F_roiMean)])
hold on
errorbar(1:4, [mean(assrFXS_M_roiMean) mean(assrFXS_F_roiMean) mean(assrHCS_M_roiMean) mean(assrHCS_F_roiMean)], ...
              [0 0 0 0], ...
              [ std(assrFXS_M_roiMean)  std(assrFXS_F_roiMean)  std(assrHCS_M_roiMean)  std(assrHCS_F_roiMean)], ...
              'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('ASSR at FCz')
ylim([0 2.5])
[tStats_assr, df_assr, pvals_assr] = statcond({assrFXS_M_roiMean' assrFXS_F_roiMean'; assrHCS_M_roiMean' assrHCS_F_roiMean'});
    % t: {[0.0531]}    {[5.9497]}    {[0.9276]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p: {[0.8182]}    {[0.0167]}    {[0.3381]}

[assrD_MaleFemale, assrD_FxsHcs] = cohensD_mainEffects2x2(assrFXS_M_roiMean', assrFXS_F_roiMean', assrHCS_M_roiMean', assrHCS_F_roiMean');
assrD_F2                          = cohensF2(tStats_assr{3}, df_assr{3}(1), df_assr{3}(2));
disp(sprintf('%.3f %.3f %.3f',assrD_MaleFemale, assrD_FxsHcs, assrD_F2))
    % Cohen 0.040 0.490 0.011

nexttile(nestedLayout4, 3)
barHandle = bar([mean(offsetFXS_M_roiMean) mean(offsetFXS_F_roiMean) mean(offsetHCS_M_roiMean) mean(offsetHCS_F_roiMean)])
hold on
errorbar(1:4, [mean(offsetFXS_M_roiMean) mean(offsetFXS_F_roiMean) mean(offsetHCS_M_roiMean) mean(offsetHCS_F_roiMean)], ...
              [0 0 0 0], ...
              [ std(offsetFXS_M_roiMean)  std(offsetFXS_F_roiMean)  std(offsetHCS_M_roiMean)  std(offsetHCS_F_roiMean)], ...
              'linestyle', 'none', 'color', [0 0 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10(\muV^2/Hz)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})
title('OffsetER at Cz')
ylim([0 1.5])
[tStats_offset, df_offset, pvals_offset] = statcond({offsetFXS_M_roiMean' offsetFXS_F_roiMean'; offsetHCS_M_roiMean' offsetHCS_F_roiMean'});
    % t: {[0.2557]}    {[8.0984]}    {[1.1652]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p:  {[0.6144]}    {[0.0055]}    {[0.2833]}

[offsetD_MaleFemale, offsetD_FxsHcs] = cohensD_mainEffects2x2(offsetFXS_M_roiMean', offsetFXS_F_roiMean', offsetHCS_M_roiMean', offsetHCS_F_roiMean');
offsetD_F2                          = cohensF2(tStats_offset{3}, df_offset{3}(1), df_offset{3}(2));
disp(sprintf('%.3f %.3f %.3f',offsetD_MaleFemale, offsetD_FxsHcs, offsetD_F2))
    % Cohen 0.091 0.567 0.013

exportgraphics(gcf, '/srv/Makoto/ASSR/p1300_erspMaleFemale/ERSPs.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')