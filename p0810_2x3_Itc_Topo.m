% 12/02/2024 Makoto. Used again.
% 10/04/2024 Makoto. Modified.
% 08/07/2024 Makoto. Modified.
% 07/10/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools/frank-pk-DataViz-3.2.3.0/daboxplot')


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


%%

% Load dummy EEG.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0100_upToDipfit/', 'loadmode', 'info');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/p0000_raw/Final_SSCT_Adults_Adolescents_editedLAD_11152024_modified.xlsx');
groupIdx = demographicData.Var6;

% Sort the subjIDs.
[subjID_sorted, sortingIdx] = sort(demographicData.Var1);
groupIdx = groupIdx(sortingIdx);

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
freqIdx = [freqIdx1 freqIdx2 freqIdx3 freqIdx4 freqIdx5 freqIdx6 freqIdx7 freqIdx8];
freqLabels = [1 2 4 8 13 20 40 80];
timeBins = -1000:3999;


%% Median ITC.

allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');

stackData_ITC = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData_ITC(:,:,:,matIdx) = elecItcMedian;
end
itcHSC = stackData_ITC(:,:,:,hcsIdx);
itcFXS = stackData_ITC(:,:,:,fxsIdx);
itcFXS_M      = stackData_ITC(:,:,:,fxsM_Idx);
itcFXS_F      = stackData_ITC(:,:,:,fxsF_Idx);
itcHCS_M      = stackData_ITC(:,:,:,hcsM_Idx);
itcHCS_F      = stackData_ITC(:,:,:,hcsF_Idx);

% Prepare the plotting data.
thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
gammaIdx = freqIdx7;

onsetIdx  = find(timeBins>=0 & timeBins<=500);
assrIdx   = find(timeBins>=0 & timeBins<=3000);
offsetIdx = find(timeBins>=3000 & timeBins<=3500);

fxsThetaOnset = squeeze(mean(mean(itcFXS(:, thetaIdx, onsetIdx, :),2),3));
hscThetaOnset = squeeze(mean(mean(itcHSC(:, thetaIdx, onsetIdx, :),2),3));

fxs40Hz = squeeze(mean(mean(itcFXS(:, gammaIdx, assrIdx, :),2),3));
hsc40Hz = squeeze(mean(mean(itcHSC(:, gammaIdx, assrIdx, :),2),3));

fxsThetaOffset = squeeze(mean(mean(itcFXS(:, thetaIdx, offsetIdx, :),2),3));
hscThetaOffset = squeeze(mean(mean(itcHSC(:, thetaIdx, offsetIdx, :),2),3));

%%
FCzIdx = 6; % FCz.
CzIdx  = 55;

figure('position', [800   100   700   800])
subplot(3,3,1)
topoplot(mean(fxsThetaOnset,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('FXS Onset Theta')

subplot(3,3,2)
topoplot(mean(fxs40Hz,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('FXS 40Hz ASSR')

subplot(3,3,3)
topoplot(mean(fxsThetaOffset,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('FXS Offset Theta')
currentPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Coh.');
set(gca,'position', currentPosition);

subplot(3,3,4)
topoplot(mean(hscThetaOnset,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('HSC Onset Theta')

subplot(3,3,5)
topoplot(mean(hsc40Hz,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('HSC 40Hz ASSR')

subplot(3,3,6)
topoplot(mean(hscThetaOffset,2), EEG.chanlocs, 'maplimits', [0 0.6])
title('HSC Offset Theta')

subplot(3,3,7)
fxsmOnset = squeeze(mean(mean(itcFXS_M(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
fxsfOnset = squeeze(mean(mean(itcFXS_F(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
hcsmOnset = squeeze(mean(mean(itcHCS_M(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
hcsfOnset = squeeze(mean(mean(itcHCS_F(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
[itcStats_onset, itcDf_onset, itcPvals_onset] = statcond({fxsmOnset' fxsfOnset'; hcsmOnset' hcsfOnset'});
    % T  {[0.1004]}    {[0.8405]}    {[0.1106]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.7521]}    {[0.3618]}    {[0.7402]}
h = daboxplot({fxsmOnset, fxsfOnset, hcsmOnset, hcsfOnset},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at Cz')
xlim([0.5 4.5])
ylim([0 1])


subplot(3,3,8)
fxsmAssr = squeeze(mean(mean(itcFXS_M(FCzIdx, gammaIdx, assrIdx, :),2),3));
fxsfAssr = squeeze(mean(mean(itcFXS_F(FCzIdx, gammaIdx, assrIdx, :),2),3));
hcsmAssr = squeeze(mean(mean(itcHCS_M(FCzIdx, gammaIdx, assrIdx, :),2),3));
hcsfAssr = squeeze(mean(mean(itcHCS_F(FCzIdx, gammaIdx, assrIdx, :),2),3));
[itcStats_assr, itcDf_assr, itcPvals_assr] = statcond({fxsmAssr' fxsfAssr'; hcsmAssr' hcsfAssr'});
    % T  {[6.2372e-05]}    {[2.8007]}    {[2.8308]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.9937]}    {[0.0978]}    {[0.0960]}
h = daboxplot({fxsmAssr, fxsfAssr, hcsmAssr, hcsfAssr},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at FCz')
xlim([0.5 4.5])
ylim([0 1])

subplot(3,3,9)
fxsmOffset = squeeze(mean(mean(itcFXS_M(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
fxsfOffset = squeeze(mean(mean(itcFXS_F(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
hcsmOffset = squeeze(mean(mean(itcHCS_M(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
hcsfOffset = squeeze(mean(mean(itcHCS_F(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
[itcStats_offset, itcDf_offset, itcPvals_offset] = statcond({fxsmOffset' fxsfOffset'; hcsmOffset' hcsfOffset'});
    % T   {[1.4061]}    {[6.0404]}    {[6.8789]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.2389]}    {[0.0159]}    {[0.0103]}
h = daboxplot({fxsmOffset, fxsfOffset, hcsmOffset, hcsfOffset},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at Cz')
xlim([0.5 4.5])
ylim([0 1])

sgtitle('ITC')

exportgraphics(gcf, '/srv/Makoto/ASSR/p0810_2x3_Itc_Topo/ITC.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')
% print('/srv/Makoto/ASSR/p0810_2x3_Itc_Topo/ERSP', '-r200', '-djpeg95')
% print('/srv/Makoto/ASSR/p0810_2x3_Itc_Topo/ERSP', '-r200', '-dsvg')