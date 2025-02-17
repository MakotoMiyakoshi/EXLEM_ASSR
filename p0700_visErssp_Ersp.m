% 12/02/2024 Makoto. Used again.
% 11/29/2024 Makoto. Used again.
% 10/25/2024 Makoto. Used.
% 10/04/2024 Makoto. Modified.
% 08/07/2024 Makoto. Modified.
% 07/09/2024 Makoto. Created.
clear
close all


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

hcsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'TDC'));
hcsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'TDC'));
fxsM_Idx = find(strcmp(sexGroupList(:,1), 'M') & strcmp(sexGroupList(:,2), 'FXS'));
fxsF_Idx = find(strcmp(sexGroupList(:,1), 'F') & strcmp(sexGroupList(:,2), 'FXS'));
hcsIdx = sort([hcsM_Idx; hcsF_Idx]);
fxsIdx = sort([fxsM_Idx; fxsF_Idx]);

% Obtain EEG IDs.
hcsM_Names = subjNames(hcsM_Idx);
hcsF_Names = subjNames(hcsF_Idx);
fxsM_Names = subjNames(fxsM_Idx);
fxsF_Names = subjNames(fxsF_Idx);
hcsNames = subjNames(hcsIdx);
fxsNames = subjNames(fxsIdx);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Plot ERSSP on power. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Median ERSP.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecErspMedian.mat');

% % Apply preselection for 13+ yrs olds.
% allMats = allMats(allAbove13Idx);

stackData = zeros(128, 100, 5001, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecErspMedian;
end
stackData      = 10*log10(stackData);
stackData_base = bsxfun(@minus, stackData, mean(stackData(:,:,1:1000,:),3));
erspHCS        = mean(stackData_base(:,:,:,hcsIdx),4);
erspFXS        = mean(stackData_base(:,:,:,fxsIdx),4);


% FXS.
pTensor = zeros(size(erspFXS));
for chIdx = 1:size(erspFXS, 1)
    for freqIdxIdx = 1:size(erspFXS,2)
        currentFreqP = squeeze(normcdf(erspFXS(chIdx,freqIdxIdx,:), squeeze(mean(erspFXS(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspFXS(chIdx,freqIdxIdx,baselineIdx),0,3))));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_FXS = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0700_visErssp_Ersp/ersp_pTensor_FXS', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')

% HCS.
pTensor = zeros(size(erspHCS));
for chIdx = 1:size(erspHCS, 1)
    for freqIdxIdx = 1:size(erspFXS,2)
        currentFreqP = squeeze(normcdf(erspHCS(chIdx,freqIdxIdx,:), squeeze(mean(erspHCS(chIdx,freqIdxIdx,baselineIdx),3)), squeeze(std(erspHCS(chIdx,freqIdxIdx,baselineIdx),0,3))));
        pTensor(chIdx,freqIdxIdx,:) = currentFreqP;
    end
end

% Determine functional ROI.
pTensor_rTale = 1-pTensor;
pTensor_sigMask = pTensor_rTale<0.01;
pTensorSum_HCS = squeeze(sum(pTensor_sigMask));

% Dump pTensor_rTale.
pTensor_rTale = single(pTensor_rTale);
save('/srv/Makoto/ASSR/p0700_visErssp_Ersp/ersp_pTensor_HCS', 'pTensor_rTale', 'timeBins', 'wtFreqBins', '-nocompression')

% Apply multiple comparison.
pTensorSum_FXS(pTensorSum_FXS<=6) = 0;
pTensorSum_HCS(pTensorSum_HCS<=6) = 0;

% Plot figures.
figure('position', [150 150 2000 500])
subplot(1,3,1)
imagesc(timeBins, 1:100, pTensorSum_FXS, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FXS grand mean (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh')

subplot(1,3,2)
imagesc(timeBins, 1:100, pTensorSum_HCS, [0 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('HCS grand mean (n=%d)', length(hcsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');

subplot(1,3,3)
diffTensor = pTensorSum_FXS-pTensorSum_HCS;
maxColorScale = max(abs(diffTensor(:)));
%imagesc(timeBins, 1:100, pTensorSum_FXS-pTensorSum_HCS, [-maxColorScale maxColorScale])
imagesc(timeBins, 1:100, pTensorSum_FXS-pTensorSum_HCS, [-128 128])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title('FXS-HCS')
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'NumCh');

print('/srv/Makoto/ASSR/p0700_visErssp_Ersp/erspSigMaskSum_groupComparison', '-dsvg')
print('/srv/Makoto/ASSR/p0700_visErssp_Ersp/erspSigMaskSum_groupComparison', '-djpeg95', '-r200')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Plot an example ERSP map at FCz to explain how ERSSP works. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FCzERSP_HCS         = squeeze(erspHCS(6,:,:));
FCzERSP_sigMask_HCS = squeeze(pTensor_sigMask(6,:,:));

figure('position', [150 150 1200 500])

subplot(1,2,1)
imagesc(timeBins, 1:100, FCzERSP_HCS, [-2 2])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('FCz, HCS grand mean (n=%d)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'dB')

subplot(1,2,2)
imagesc(timeBins, 1:100, FCzERSP_sigMask_HCS, [-1 1])
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
colormap jet
axis xy
axis square
title(sprintf('Significance mask (p<0.01)', length(fxsIdx)))
set(gca, 'ytick', freqIdx, 'yticklabel', freqLabels)
line([0 0],       ylim, 'color', [0 0 0], 'linewidth', 2)
line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
xlim([-800 3800])
xlabel('Latency (ms)')
ylabel('Frequency (Hz)')
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', 'Sig.')

print('/srv/Makoto/ASSR/p0700_visErssp_Ersp/HcsFCzERSP', '-dsvg')
print('/srv/Makoto/ASSR/p0700_visErssp_Ersp/HcsFCzERSP', '-djpeg95', '-r200')


%% Ordinary two-sample t-test on 40-Hz ASSR at FCz.
fdxGamma = squeeze(mean(stackData_base(6,80,assrIdx,fxsIdx),3));
tdcGamma = squeeze(mean(stackData_base(6,80,assrIdx,hcsIdx),3));

[H,P,CI,STATS] = ttest2(fdxGamma, tdcGamma);
disp(sprintf('t(%d)=%.1f, p=%.4f', STATS.df, STATS.tstat, P)) % t(66)=-2.1, p=0.0422


%% Very well, then it's LME. First, with Gamma.
        % fdxGamma = squeeze(mean(stackData_base(6,80,assrIdx,fxsIdx),3));
        % tdcGamma = squeeze(mean(stackData_base(6,80,assrIdx,hcsIdx),3));
        % 
        % % allAbove13Idx     = find(ageVec>12);
        % % subjNames_13plus  = subjNames( allAbove13Idx)';
        % % groupNames_13plus = groupNames(allAbove13Idx)';
        % % ageVec_13plus     = ageVec(    allAbove13Idx)';
        % 
        % fxsIdx_preselected = find(contains(groupNames_13plus, 'FXS'));
        % hcsIdx_preselected = find(contains(groupNames_13plus, 'TDC'));
        % 
        % fxsAge_preselected = ageVec_13plus(fxsIdx_preselected);
        % hcsAge_preselected = ageVec_13plus(hcsIdx_preselected);
        % 
        % % Combine data to a table.
        % SubjectIdx    = [1:length(fdxGamma)+length(tdcGamma)]';
        % Group         = [repmat({'FXS'}, [length(fdxGamma), 1]); repmat({'HCS'}, [length(tdcGamma), 1])];
        % Age           = [fxsAge_preselected'; hcsAge_preselected'];
        % meanGamma     = [fdxGamma; tdcGamma];
        % gammaTable = table(SubjectIdx, Group, Age, meanGamma, 'variableNames', {'SubjectIdx' 'Group', 'Age', 'meanGamma'});
        % gammaTable.SubjectIdx = categorical(gammaTable.SubjectIdx);
        % gammaTable.Group      = categorical(gammaTable.Group);
        % 
        % % Start LME.
        % lmeModel1 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + (1|SubjectIdx)'); % t(25) = 2.1743, p=0.039
        % lmeModel2 = fitlme(gammaTable, 'meanGamma ~ 1 + Group * Age + (1|SubjectIdx)');
        % % lmeModel3 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + Group * Age + (1|SubjectIdx)'); % meaningless.
        % 
        % compare(lmeModel1, lmeModel2)
        % 
        % lmeModel1.Rsquared % R^2 0.7525 Yey!
        % sqrt(lmeModel1.Rsquared.Adjusted) % r=0.8675.
        % 
        % % The stats structure contains p-values and other statistics
        % [fixedEffectsEstimates, fixedEffectsNames, stats] = fixedEffects(lmeModel1);
        % groupDiffIdx = find(contains(fixedEffectsNames.Name, 'Group_HCS'));
        % pValues = stats.pValue(groupDiffIdx);
        % disp(pValues) % 0.0394.


%% Second, with onset-evoked theta.
fxsTheta1 = squeeze(mean(mean(stackData_base(6,thetaIdx,onsetIdx, fxsIdx),2),3));
hcsTheta1 = squeeze(mean(mean(stackData_base(6,thetaIdx,onsetIdx, hcsIdx),2),3));

[H,P,CI,STATS] = ttest2(fxsTheta1, hcsTheta1);
disp(sprintf('t(%d)=%.1f, p=%.4f', STATS.df, STATS.tstat, P)) % t(66)=1.9, p=0.0647


        % % Combine data to a table.
        % meanTheta1 = [fxsTheta1; hcsTheta1];
        % 
        % Theta1Table = table(SubjectIdx, Group, Age, meanTheta1, 'variableNames', {'SubjectIdx' 'Group', 'Age', 'meanTheta1'});
        % Theta1Table.SubjectIdx = categorical(Theta1Table.SubjectIdx);
        % Theta1Table.Group      = categorical(Theta1Table.Group);
        % 
        % % Start LME.
        % lmeModel1 = fitlme(Theta1Table, 'meanTheta1 ~ 1 + Group + Age + (1|SubjectIdx)');
        % lmeModel2 = fitlme(Theta1Table, 'meanTheta1 ~ 1 + Group * Age + (1|SubjectIdx)');
        % % lmeModel3 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + Group * Age + (1|SubjectIdx)'); % meaningless.
        % 
        % compare(lmeModel1, lmeModel2)
        % 
        % lmeModel1.Rsquared % R^2 0.6163
        % sqrt(lmeModel1.Rsquared.Adjusted) % r=0.7850
        % 
        % lmeModel2.Rsquared % R^2 0.7362
        % sqrt(lmeModel2.Rsquared.Adjusted) % r=0.8580
        % 
        % % The stats structure contains p-values and other statistics
        % [fixedEffectsEstimates, fixedEffectsNames, stats] = fixedEffects(lmeModel1);
        % groupDiffIdx = find(contains(fixedEffectsNames.Name, 'Group_HCS'));
        % pValues = stats.pValue(groupDiffIdx);
        % disp(pValues) % 0.0271.


%% Finally, with offset-evoked theta.
fxsTheta2 = squeeze(mean(mean(stackData_base(6,thetaIdx,offsetIdx, fxsIdx),2),3));
hcsTheta2 = squeeze(mean(mean(stackData_base(6,thetaIdx,offsetIdx, hcsIdx),2),3));

[H,P,CI,STATS] = ttest2(fxsTheta2, hcsTheta2);
disp(sprintf('t(%d)=%.1f, p=%.4f', STATS.df, STATS.tstat, P)) % t(66)=2.0, p=0.0484


        % % Combine data to a table.
        % meanTheta2 = [fxsTheta2; hcsTheta2];
        % 
        % Theta2Table = table(SubjectIdx, Group, Age, meanTheta2, 'variableNames', {'SubjectIdx' 'Group', 'Age', 'meanTheta2'});
        % Theta2Table.SubjectIdx = categorical(Theta2Table.SubjectIdx);
        % Theta2Table.Group      = categorical(Theta2Table.Group);
        % 
        % % Start LME.
        % lmeModel1 = fitlme(Theta2Table, 'meanTheta2 ~ 1 + Group + Age + (1|SubjectIdx)');
        % lmeModel2 = fitlme(Theta2Table, 'meanTheta2 ~ 1 + Group * Age + (1|SubjectIdx)');
        % % lmeModel3 = fitlme(gammaTable, 'meanGamma ~ 1 + Group + Age + Group * Age + (1|SubjectIdx)'); % meaningless.
        % 
        % compare(lmeModel1, lmeModel2)
        % 
        % lmeModel1.Rsquared % R^2 0.6076
        % sqrt(lmeModel1.Rsquared.Adjusted) % r=0.7795
        % 
        % lmeModel2.Rsquared % R^2 0.8028
        % sqrt(lmeModel2.Rsquared.Adjusted) % r=0.8960
        % 
        % % The stats structure contains p-values and other statistics
        % [fixedEffectsEstimates, fixedEffectsNames, stats] = fixedEffects(lmeModel1);
        % groupDiffIdx = find(contains(fixedEffectsNames.Name, 'Group_HCS'));
        % pValues = stats.pValue(groupDiffIdx);
        % disp(pValues) % 0.0575.