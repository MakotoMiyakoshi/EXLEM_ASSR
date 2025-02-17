% 01/17/2025 Makoto. Changing from bar graphs to box plots.
% 01/10/2025 Makoto. Eureka proven.
% 01/09/2025 Makoto. Eureka.
% 01/08/2025 Makoto. Updated.
% 12/11/2024 Makoto. Created.
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
    sexGroupList{setIdx,1} = EEG.etc.Sex{1,1};
    sexGroupList{setIdx,2} = EEG.etc.Dx{1,1};
end

subjNames  = cellfun(@(x) x(1:4), {allSets.name},    'UniformOutput', false);
groupNames = cellfun(@(x) x(),    sexGroupList(:,2), 'UniformOutput', false);
    %{
    figure
    bar(sort(ageVec))
    %}

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
fxsIdx   = find(strcmp(sexGroupList(:,2), 'FXS'));
hcsIdx   = find(strcmp(sexGroupList(:,2), 'TDC'));

fxs_demographics = [subjNames(fxsIdx)' sexGroupList(fxsIdx,:)];
hcs_demographics = [subjNames(hcsIdx)' sexGroupList(hcsIdx,:)];

% Obtain EEG IDs.
fxsM_Names = subjNames(fxsM_Idx);
fxsF_Names = subjNames(fxsF_Idx);
hcsM_Names = subjNames(hcsM_Idx);
hcsF_Names = subjNames(hcsF_Idx);
fxsNames = subjNames(fxsIdx);
hcsNames = subjNames(hcsIdx);

%% Load EEG data.
load('/srv/Makoto/ASSR/p0210_epochForERP/meanMedianERPs.mat')
meanERP   = meanERP(:,:,preselectingIdx2);
medianERP = medianERP(:,:,preselectingIdx2);

czIdx = 55; % Cz.
erpTimes = -1000:3999;

% Identifying the artifact within the first birst.
close all

figure('position', [50 50 800 1200])
mainLayout = tiledlayout(6,3,"TileSpacing", "tight", "Padding", "tight");

% ERP.
czIdx  = 55; % Cz.
fczIdx = 6;
fxsM = squeeze(medianERP(:, :, fxsM_Idx));
fxsF = squeeze(medianERP(:, :, fxsF_Idx));
hcsM = squeeze(medianERP(:, :, hcsM_Idx));
hcsF = squeeze(medianERP(:, :, hcsF_Idx));
fxsERP = medianERP(:,:,fxsIdx);
hcsERP = medianERP(:,:,hcsIdx);

% fxsM = squeeze(meanERP(:, :, fxsM_Idx));
% fxsF = squeeze(meanERP(:, :, fxsF_Idx));
% hcsM = squeeze(meanERP(:, :, hcsM_Idx));
% hcsF = squeeze(meanERP(:, :, hcsF_Idx));
% fxsERP = meanERP(:,:,fxsIdx);
% hcsERP = meanERP(:,:,hcsIdx);

nexttile(mainLayout, [1 3])
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsERP(czIdx,:,:)),2) mean(squeeze(hcsERP(czIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([-200 3500])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
title('Cz Grand-mean ERP')
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS' 'HCS'}, 'Location', 'best')

nexttile(mainLayout, 4)
plotData1 = mean(squeeze(fxsERP(czIdx,:,:)),2);
plotData2 = mean(squeeze(hcsERP(czIdx,:,:)),2);

baselineIdx = find(erpTimes >= -5 & erpTimes<=0);
plotData1 = plotData1-mean(plotData1(baselineIdx));
plotData2 = plotData2-mean(plotData2(baselineIdx));

linePlotHandles = plot(erpTimes, [plotData1 plotData2]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([-7 13])
ylim([-3.2 3.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([2 2], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS' 'HCS'}, 'Location', 'best')

nestedLayout1 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 5;

plotData1 = mean(fxsERP,3);
plotData2 = mean(hcsERP,3);
baselineIdx = find(erpTimes >= -5 & erpTimes<=0);
plotData1 = plotData1-mean(plotData1(:,baselineIdx),2);
plotData2 = plotData2-mean(plotData2(:,baselineIdx),2);

nexttile(nestedLayout1, 1)
topoplot(plotData1(:,1003), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout1, 2)
topoplot(plotData2(:,1003), EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 6)
meanForBar = [mean(squeeze(fxsERP(55,1003,:))) mean(squeeze(hcsERP(55,1003,:)))];
seForBar  = [std(squeeze(fxsERP(55,1003,:)))/sqrt(size(fxsERP,3))  std(squeeze(hcsERP(55,1003,:)))/sqrt(size(hcsERP,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, seForBar, [0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-3 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS' 'HCS'})

statsData1 = fxsERP;
statsData2 = hcsERP;
statsData1 = bsxfun(@minus, statsData1, mean(statsData1(:,baselineIdx,:),2));
statsData2 = bsxfun(@minus, statsData2, mean(statsData2(:,baselineIdx,:),2));
[H1,P1,CI1,STATS1] = ttest2(squeeze(statsData1(55,1003,:))', squeeze(statsData2(55,1003,:))');
d1 = cohensD(squeeze(fxsERP(55,1003,:))', squeeze(hcsERP(55,1003,:))');
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS1.df, STATS1.tstat, P1, d1))
    %t(60)=1.712, p=0.092, Cohen's d = 0.501

nexttile(mainLayout, 7)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsERP(czIdx,:,:)),2) mean(squeeze(hcsERP(czIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([-100 400])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([98 98], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS' 'HCS'}, 'Location', 'best')

nestedLayout2 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 8;
nexttile(nestedLayout2, 1)
topoplot(mean(squeeze(fxsERP(:,1100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout2, 2)
topoplot(mean(squeeze(hcsERP(:,1100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 9)
meanForBar = [mean(squeeze(fxsERP(55,1100,:))) mean(squeeze(hcsERP(55,1100,:)))];
seForBar  = [std(squeeze(fxsERP(55,1100,:)))/sqrt(size(fxsERP,3))  std(squeeze(hcsERP(55,1100,:)))/sqrt(size(hcsERP,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, seForBar, [0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-2.2 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS' 'HCS'})

[H2,P2,CI2,STATS2] = ttest2(squeeze(fxsERP(55,1100,:))', squeeze(hcsERP(55,1100,:))');
d2 = cohensD(squeeze(fxsERP(55,1100,:))', squeeze(hcsERP(55,1100,:))');
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS2.df, STATS2.tstat, P2, d2))
    % t(60)=-1.750, p=0.085, Cohen's d = -0.445.

% 01/09/2025 Added. 
nexttile(mainLayout, 10)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsERP(czIdx,:,:)),2) mean(squeeze(hcsERP(czIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([500 3000])
ylim([-1 1])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([98 98], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS' 'HCS'}, 'Location', 'best')

fxsAssrStdTopo = squeeze(mean(std(fxsERP(:,1501:4000,:),0,2),3));
hcsAssrStdTopo = squeeze(mean(std(hcsERP(:,1501:4000,:),0,2),3));

nestedLayout2 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 11;
nexttile(nestedLayout2, 1)
topoplot(fxsAssrStdTopo, EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout2, 2)
topoplot(hcsAssrStdTopo, EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 12)
meanForBar = [mean(squeeze(std(fxsERP(55,1501:4000,:),0,2))) mean(squeeze(std(hcsERP(55,1501:4000,:),0,2)))];
seForBar   = [std(squeeze(std(fxsERP(55,1501:4000,:),0,2)))/sqrt(length(fxsIdx)) std(squeeze(std(hcsERP(55,1501:4000,:),0,2)))/sqrt(length(hcsIdx))];
barHandle  = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, [0 0], seForBar, 'linestyle', 'none', 'color', [0 0 0])
ylim([0 1.2])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS' 'HCS'})

fxsAssrStd = squeeze(std(fxsERP(55,1501:4000,:),0,2));
hcsAssrStd = squeeze(std(hcsERP(55,1501:4000,:),0,2));

[H3,P3,CI3,STATS3] = ttest2(fxsAssrStd, hcsAssrStd);
d3 = cohensD(fxsAssrStd, hcsAssrStd);
disp(sprintf('t(%d)=%.3f, p=%.6f, Cohen''s d = %.3f.', STATS3.df, STATS3.tstat, P3, d3))
    % t(60)=4.937, p=0.000007, Cohen's d = 1.257.


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute FCz ERP. %%%
%%%%%%%%%%%%%%%%%%%%%%%%
fxsAssrErp = squeeze(fxsERP(6,1501:4000,:));
hcsAssrErp = squeeze(hcsERP(6,1501:4000,:));
fxsAssrErp = detrend(fxsAssrErp, 'linear');
hcsAssrErp = detrend(hcsAssrErp, 'linear');

% Calculate PSDs.
[spectra1,freqs] = spectopo(fxsAssrErp', 0, EEG.srate, 'winsize', 1000, 'freqfac', 1, 'plot', 'off');
[spectra2,freqs] = spectopo(hcsAssrErp', 0, EEG.srate, 'winsize', 1000, 'freqfac', 1, 'plot', 'off');

% Plot PSD.
nestedLayout5 = tiledlayout(mainLayout, 1, 4, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout5.Layout.Tile     = 13;
nestedLayout5.Layout.TileSpan = [1 3];
nexttile(nestedLayout5, 1)

linePlotHandles = plot(freqs, [mean(spectra1)' mean(spectra2)']);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([1 55])
xlabel('Frequency (Hz)')
ylabel('10*log10\muV^2/Hz')
legend({'FXS' 'HCS'}, 'location','southwest')

nexttile(nestedLayout5, 2)
assrFreqBinIdx = find(freqs==40);
fxs40HzPower = spectra1(:,assrFreqBinIdx);
hcs40HzPower = spectra2(:,assrFreqBinIdx);

meanForBar = [mean(fxs40HzPower) mean(hcs40HzPower)];
seForBar   = [std(fxs40HzPower)/sqrt(length(fxs40HzPower)) std(hcs40HzPower)/sqrt(length(hcs40HzPower))];
barHandle  = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, seForBar, [0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-12.5 -8.5])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10\muV^2/Hz')
xticklabels({'FXS' 'HCS'})
title('40Hz Peak')

[H4,P4,CI4,STATS4] = ttest2(fxs40HzPower, hcs40HzPower);
d4 = cohensD(fxs40HzPower, hcs40HzPower);
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS4.df, STATS4.tstat, P4, d4))
    % t(60)=2.546, p=0.013, Cohen's d = 0.648.

nexttile(nestedLayout5, 3)
belowLineFreqIdx = find(freqs>=1 & freqs<=55);
fxsSlopePower = mean(spectra1(:,belowLineFreqIdx),2);
hcsSlopePower = mean(spectra2(:,belowLineFreqIdx),2);
meanForBar = [mean(fxsSlopePower) mean(hcsSlopePower)];
seForBar   = [std(fxsSlopePower)/sqrt(length(fxsSlopePower)) std(hcsSlopePower)/sqrt(length(hcsSlopePower))];
barHandle  = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, seForBar, [0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-27 -23])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10\muV^2/Hz')
xticklabels({'FXS' 'HCS'})
title('1-55Hz Slope')

[H5,P5,CI5,STATS5] = ttest2(fxsSlopePower, hcsSlopePower);
d5 = cohensD(fxsSlopePower, hcsSlopePower);
disp(sprintf('t(%d)=%.3f, p=%.7f, Cohen''s d = %.3f.', STATS5.df, STATS5.tstat, P5, d5))
    % t(60)=5.590, p=0.0000006, Cohen's d = 1.423.

nexttile(nestedLayout5, 4)
belowLineFreqIdx = find(freqs>=1 & freqs<=55);
fxsSlopePower = mean(spectra1(:,belowLineFreqIdx),2);
hcsSlopePower = mean(spectra2(:,belowLineFreqIdx),2);
meanForBar = [mean(fxs40HzPower-fxsSlopePower) mean(hcs40HzPower-hcsSlopePower)];
seForBar   = [std(fxs40HzPower-fxsSlopePower)/sqrt(length(fxs40HzPower-fxsSlopePower)) std(hcs40HzPower-hcsSlopePower)/sqrt(length(hcs40HzPower-hcsSlopePower))];
barHandle  = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, [0 0], seForBar, 'linestyle', 'none', 'color', [0 0 0])
ylim([12 16])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('10*log10\muV^2/Hz')
xticklabels({'FXS' 'HCS'})
title('Peak-Slope')

[H6,P6,CI6,STATS6] = ttest2(fxs40HzPower-fxsSlopePower, hcs40HzPower-hcsSlopePower);
d6 = cohensD(fxs40HzPower-fxsSlopePower, hcs40HzPower-hcsSlopePower);
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS6.df, STATS6.tstat, P6, d6))
    % t(60)=-0.882, p=0.381, Cohen's d = -0.225.

nexttile(mainLayout, 16)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsERP(czIdx,:,:)),2) mean(squeeze(hcsERP(czIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [0 0 1 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = '-';
xlim([2900 3400])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([3100 3100], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS' 'HCS'}, 'Location', 'best')

nestedLayout3 = tiledlayout(mainLayout, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout3.Layout.Tile = 17;
nexttile(nestedLayout3, 1)
topoplot(mean(squeeze(fxsERP(:,4096,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout3, 2)
topoplot(mean(squeeze(hcsERP(:,4096,:)),2), EEG.chanlocs, 'maplimits', [-3 3])


nexttile(mainLayout, 18)
meanForBar = [mean(squeeze(fxsERP(55,4096,:))) mean(squeeze(hcsERP(55,4096,:)))];
seForBar  = [std(squeeze(fxsERP(55,4096,:)))/sqrt(size(fxsM,3))  std(squeeze(hcsERP(55,4096,:)))/sqrt(size(fxsF,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:2, meanForBar, seForBar, [0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-3 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS' 'HCS'})

[H7,P7,CI7,STATS7] = ttest2(squeeze(fxsERP(55,4096,:))', squeeze(hcsERP(55,4096,:))');
d7 = cohensD(squeeze(fxsERP(55,4096,:))', squeeze(hcsERP(55,4096,:))');
disp(sprintf('t(%d)=%.3f, p=%.3f, Cohen''s d = %.3f.', STATS7.df, STATS7.tstat, P7, d7))
    % t(60)=-3.139, p=0.003, Cohen's d = -0.799.

%% Export the results.
combinedDemographics = [fxs_demographics; hcs_demographics];
fmrp     = [fmrpVec(fxsIdx); fmrpVec(hcsIdx)];
subjID   = combinedDemographics(:,1);
sex      = combinedDemographics(:,2);
group    = combinedDemographics(:,3);
onsetER  = [squeeze(fxsERP(55,1100,:)); squeeze(hcsERP(55,1100,:))];
ASSR     = [fxsAssrStd; hcsAssrStd];
offsetER = [squeeze(fxsERP(55,4096,:)); squeeze(hcsERP(55,4096,:))];
psdPeak  = [fxs40HzPower; hcs40HzPower];
psdSlope = [fxsSlopePower; hcsSlopePower];
psdPeakSlopeDiff = [fxs40HzPower-fxsSlopePower; hcs40HzPower-hcsSlopePower];
resultTableErp = table(subjID, sex, group, onsetER, ASSR, offsetER, psdPeak, psdSlope, psdPeakSlopeDiff);

exportgraphics(gcf, '/srv/Makoto/ASSR/p1412_ERP/ERPs.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')
save('/srv/Makoto/ASSR/p1412_ERP/resultTableErp.mat', 'resultTableErp')