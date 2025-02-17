% 12/10/2024 Makoto. Added male-female comparison.
% 12/02/2024 Makoto. Used again.
% 10/30/2024 Makoto. Delta band power analysis added.
% 08/07/2024 Makoto. Modified.
% 07/10/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all

addpath('/srv/Makoto/Tools/frank-pk-DataViz-3.2.3.0/daboxplot')


%% Preselect datasets recorded from 13+ yrs old subjects.
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
timeBins = -1000:4000;


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
erspHCS        = stackData_base(:,:,:,hcsIdx);
erspFXS        = stackData_base(:,:,:,fxsIdx);
erspFXS_M      = stackData_base(:,:,:,fxsM_Idx);
erspFXS_F      = stackData_base(:,:,:,fxsF_Idx);
erspHCS_M      = stackData_base(:,:,:,hcsM_Idx);
erspHCS_F      = stackData_base(:,:,:,hcsF_Idx);

% Prepare the plotting data.
deltaIdx = find(wtFreqBins>=1 & wtFreqBins<=2);
thetaIdx = find(wtFreqBins>=3 & wtFreqBins<=8);
deltaThetaAlphaIdx = find(wtFreqBins>=2 & wtFreqBins<=13);
gammaIdx = freqIdx7;

onsetIdx  = find(timeBins>=0 & timeBins<=500);
assrIdx   = find(timeBins>=0 & timeBins<=3000);
offsetIdx = find(timeBins>=3000 & timeBins<=3500);

fxsThetaOnset = squeeze(mean(mean(erspFXS(:, thetaIdx, onsetIdx, :),2),3));
hcsThetaOnset = squeeze(mean(mean(erspHCS(:, thetaIdx, onsetIdx, :),2),3));

fxs40Hz = squeeze(mean(mean(erspFXS(:, gammaIdx, assrIdx, :),2),3));
hcs40Hz = squeeze(mean(mean(erspHCS(:, gammaIdx, assrIdx, :),2),3));

% fxsDelta = squeeze(mean(mean(erspFXS(:, deltaIdx, assrIdx, :),2),3));
% hcsDelta = squeeze(mean(mean(erspHCS(:, deltaIdx, assrIdx, :),2),3));

fxsThetaOffset = squeeze(mean(mean(erspFXS(:, thetaIdx, offsetIdx, :),2),3));
hcsThetaOffset = squeeze(mean(mean(erspHCS(:, thetaIdx, offsetIdx, :),2),3));

%%
FCzIdx = 6; % FCz.
CzIdx  = 55;

figure('position', [800   100   700   800])
subplot(3,3,1)
topoplot(mean(fxsThetaOnset,2), EEG.chanlocs, 'maplimits', [-1 1])
title('FXS Onset Theta')

subplot(3,3,2)
topoplot(mean(fxs40Hz,2), EEG.chanlocs, 'maplimits', [-1 1])
title('FXS 40Hz ASSR')

% subplot(3,3,3)
% topoplot(mean(fxsDelta,2), EEG.chanlocs, 'maplimits', [-1 1])
% title('FXS Delta')

subplot(3,3,3)
topoplot(mean(fxsThetaOffset,2), EEG.chanlocs, 'maplimits', [-1 1])
title('FXS Offset Theta')
currentPosition = get(gca,'position');
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB');
set(gca,'position', currentPosition);


subplot(3,3,4)
topoplot(mean(hcsThetaOnset,2), EEG.chanlocs, 'maplimits', [-1 1])
title('HCS Onset Theta')

subplot(3,3,5)
topoplot(mean(hcs40Hz,2), EEG.chanlocs, 'maplimits', [-1 1])
title('HCS 40Hz ASSR')

% subplot(3,3,7)
% topoplot(mean(hcsDelta,2), EEG.chanlocs, 'maplimits', [-1 1])
% title('HCS Delta')

subplot(3,3,6)
topoplot(mean(hcsThetaOffset,2), EEG.chanlocs, 'maplimits', [-1 1])
title('HCS Offset Theta')


subplot(3,3,7)
fxsmOnset = squeeze(mean(mean(erspFXS_M(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
fxsfOnset = squeeze(mean(mean(erspFXS_F(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
hcsmOnset = squeeze(mean(mean(erspHCS_M(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
hcsfOnset = squeeze(mean(mean(erspHCS_F(CzIdx, deltaThetaAlphaIdx, onsetIdx, :),2),3));
[erspStats_onset, erspDf_onset, erspPvals_onset] = statcond({fxsmOnset' fxsfOnset'; hcsmOnset' hcsfOnset'});
    % T  {[10.4877]}    {[10.2977]}    {[2.6117]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.0017]}    {[0.0019]}    {[0.1097]}
h = daboxplot({fxsmOnset, fxsfOnset, hcsmOnset, hcsfOnset},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at Cz')
xlim([0.5 4.5])
ylim([-0.5 2])


subplot(3,3,8)
fxsmAssr = squeeze(mean(mean(erspFXS_M(FCzIdx, gammaIdx, assrIdx, :),2),3));
fxsfAssr = squeeze(mean(mean(erspFXS_F(FCzIdx, gammaIdx, assrIdx, :),2),3));
hcsmAssr = squeeze(mean(mean(erspHCS_M(FCzIdx, gammaIdx, assrIdx, :),2),3));
hcsfAssr = squeeze(mean(mean(erspHCS_F(FCzIdx, gammaIdx, assrIdx, :),2),3));
[erspStats_assr, erspDf_assr, erspPvals_assr] = statcond({fxsmAssr' fxsfAssr'; hcsmAssr' hcsfAssr'});
    % T  {[0.0530]}    {[5.9534]}    {[0.9262]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.8184]}    {[0.0167]}    {[0.3385]}
h = daboxplot({fxsmAssr, fxsfAssr, hcsmAssr, hcsfAssr},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at FCz')
xlim([0.5 4.5])
ylim([0 3.5])


subplot(3,3,9)
fxsmOffset = squeeze(mean(mean(erspFXS_M(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
fxsfOffset = squeeze(mean(mean(erspFXS_F(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
hcsmOffset = squeeze(mean(mean(erspHCS_M(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
hcsfOffset = squeeze(mean(mean(erspHCS_F(CzIdx, deltaThetaAlphaIdx, offsetIdx, :),2),3));
[erspStats_offset, erspDf_offset, erspPvals_offset] = statcond({fxsmOffset' fxsfOffset'; hcsmOffset' hcsfOffset'});
    % T  {[0.2520]}    {[8.0885]}    {[1.1687]}
    % DF {[1 88]}    {[1 88]}    {[1 88]}
    % P  {[0.6169]}    {[0.0055]}    {[0.2826]}
h = daboxplot({fxsmOffset, fxsfOffset, hcsmOffset, hcsfOffset},'linkline', 1, 'xtlabels', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'}, 'whiskers',0,'outliers',0,'outsymbol','r*','scatter',2,'boxalpha',0.6);
h.bx(1).FaceColor = [0.8500 0.3250 0.0980];
h.bx(2).FaceColor = [0.8500 0.3250 0.0980];
h.bx(3).FaceColor = [0 0.4470 0.7410];
h.bx(4).FaceColor = [0 0.4470 0.7410];
ylabel('dB Power at Cz')
xlim([0.5 4.5])
ylim([-0.5 2.5])

sgtitle('ERSP')

exportgraphics(gcf, '/srv/Makoto/ASSR/p0710_2x3_Ersp_Topo/ERSP.pdf', 'ContentType', 'vector')
% print('/srv/Makoto/ASSR/p0710_2x3_ERSP_Topo/ERSP', '-r200', '-djpeg95')
% print('/srv/Makoto/ASSR/p0710_2x3_ERSP_Topo/ERSP', '-r200', '-dsvg')