% 12/11/2024 Makoto. Created.
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

%% Load EEG data.
load('/srv/Makoto/ASSR/p0210_epochForERP/meanMedianERPs.mat')

chIdx = 55; % Cz.
erpTimes = -1000:3999;

% Identifying the artifact within the first birst.
close all

figure('position', [50 50 800 1200])
mainLayout = tiledlayout(4,3,"TileSpacing", "tight", "Padding", "tight");

% ERP.
chIdx = 55; % Cz.
fxsM = squeeze(meanERP(:, :, fxsM_Idx));
fxsF = squeeze(meanERP(:, :, fxsF_Idx));
hcsM = squeeze(meanERP(:, :, hcsM_Idx));
hcsF = squeeze(meanERP(:, :, hcsF_Idx));

nexttile(mainLayout, [1 3])
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsM(chIdx,:,:)),2) mean(squeeze(fxsF(chIdx,:,:)),2) mean(squeeze(hcsM(chIdx,:,:)),2) mean(squeeze(hcsF(chIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [1 0 0 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = ':';
linePlotHandles(3).Color = [0 0 1 0.5];
linePlotHandles(3).LineWidth = 1;
linePlotHandles(3).LineStyle = '-';
linePlotHandles(4).Color = [0 0 1 0.5];
linePlotHandles(4).LineWidth = 1;
linePlotHandles(4).LineStyle = ':';
xlim([-200 3500])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
title('Cz Grand-mean ERP')
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS M' 'FXS F' 'HCS M' 'HCS F'}, 'Location', 'best')

nexttile(mainLayout, 4)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsM(chIdx,:,:)),2) mean(squeeze(fxsF(chIdx,:,:)),2) mean(squeeze(hcsM(chIdx,:,:)),2) mean(squeeze(hcsF(chIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [1 0 0 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = ':';
linePlotHandles(3).Color = [0 0 1 0.5];
linePlotHandles(3).LineWidth = 1;
linePlotHandles(3).LineStyle = '-';
linePlotHandles(4).Color = [0 0 1 0.5];
linePlotHandles(4).LineWidth = 1;
linePlotHandles(4).LineStyle = ':';
xlim([-7 13])
ylim([-3.2 3.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([2 2], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS M' 'FXS F' 'HCS M' 'HCS F'}, 'Location', 'best')

nestedLayout1 = tiledlayout(mainLayout, 2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout1.Layout.Tile = 5;
nexttile(nestedLayout1, 1)
topoplot(mean(squeeze(fxsM(:,1003,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout1, 2)
topoplot(mean(squeeze(fxsF(:,1003,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout1, 3)
topoplot(mean(squeeze(hcsM(:,1003,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout1, 4)
topoplot(mean(squeeze(hcsF(:,1003,:)),2), EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 6)
meanForBar = [mean(squeeze(fxsM(55,1003,:))) mean(squeeze(fxsF(55,1003,:))) mean(squeeze(hcsM(55,1003,:))) mean(squeeze(hcsF(55,1003,:)))];
seForBar  = [std(squeeze(fxsM(55,1003,:)))/sqrt(size(fxsM,3))  std(squeeze(fxsF(55,1003,:)))/sqrt(size(fxsF,3))  std(squeeze(hcsM(55,1003,:)))/sqrt(size(hcsM,3))  std(squeeze(hcsF(55,1003,:)))/sqrt(size(hcsF,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:4, meanForBar, seForBar, [0 0 0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-5 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})

[stats1, df1, pvals1] = statcond({squeeze(fxsM(55,1003,:))' squeeze(fxsF(55,1003,:))'; squeeze(hcsM(55,1003,:))' squeeze(hcsF(55,1003,:))'});
    % tStats: {[3.9160]}    {[11.1690]}    {[0.6506]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p:  {[0.0510]}    {[0.0012]}    {[0.4221]}

[stats1_MaleFemale, stats1_FxsHcs] = cohensD_mainEffects2x2(squeeze(fxsM(55,1003,:))', squeeze(fxsF(55,1003,:))', squeeze(hcsM(55,1003,:))', squeeze(hcsF(55,1003,:))');
stats1_F2                          = cohensF2(stats1{3}, df1{3}(1), df1{3}(2));
disp(sprintf('%.3f %.3f %.3f', stats1_MaleFemale, stats1_FxsHcs, stats1_F2))
    % Cohen 0.341 0.544 0.007

nexttile(mainLayout, 7)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsM(chIdx,:,:)),2) mean(squeeze(fxsF(chIdx,:,:)),2) mean(squeeze(hcsM(chIdx,:,:)),2) mean(squeeze(hcsF(chIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [1 0 0 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = ':';
linePlotHandles(3).Color = [0 0 1 0.5];
linePlotHandles(3).LineWidth = 1;
linePlotHandles(3).LineStyle = '-';
linePlotHandles(4).Color = [0 0 1 0.5];
linePlotHandles(4).LineWidth = 1;
linePlotHandles(4).LineStyle = ':';
xlim([-100 400])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([98 98], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS M' 'FXS F' 'HCS M' 'HCS F'}, 'Location', 'best')

nestedLayout2 = tiledlayout(mainLayout, 2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout2.Layout.Tile = 8;
nexttile(nestedLayout2, 1)
topoplot(mean(squeeze(fxsM(:,1097,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout2, 2)
topoplot(mean(squeeze(fxsF(:,1097,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout2, 3)
topoplot(mean(squeeze(hcsM(:,1097,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout2, 4)
topoplot(mean(squeeze(hcsF(:,1097,:)),2), EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 9)
meanForBar = [mean(squeeze(fxsM(55,1097,:))) mean(squeeze(fxsF(55,1097,:))) mean(squeeze(hcsM(55,1097,:))) mean(squeeze(hcsF(55,1097,:)))];
seForBar  = [std(squeeze(fxsM(55,1097,:)))/sqrt(size(fxsM,3))  std(squeeze(fxsF(55,1097,:)))/sqrt(size(fxsF,3))  std(squeeze(hcsM(55,1097,:)))/sqrt(size(hcsM,3))  std(squeeze(hcsF(55,1097,:)))/sqrt(size(hcsF,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:4, meanForBar, seForBar, [0 0 0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-3 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})

[stats2, df2, pvals2] = statcond({squeeze(fxsM(55,1097,:))' squeeze(fxsF(55,1097,:))'; squeeze(hcsM(55,1097,:))' squeeze(hcsF(55,1097,:))'});
    % tStats: {[0.0990]}    {[6.2311]}    {[0.6387]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p:  {[0.7538]}    {[0.0144]}    {[0.4263]}

[stats2_MaleFemale, stats2_FxsHcs] = cohensD_mainEffects2x2(squeeze(fxsM(55,1097,:))', squeeze(fxsF(55,1097,:))', squeeze(hcsM(55,1097,:))', squeeze(hcsF(55,1097,:))');
stats2_F2                          = cohensF2(stats2{3}, df2{3}(1), df2{3}(2));
disp(sprintf('%.3f %.3f %.3f', stats2_MaleFemale, stats2_FxsHcs, stats2_F2))
    % Cohen 0.065 0.398 0.007

nexttile(mainLayout, 10)
linePlotHandles = plot(erpTimes, [mean(squeeze(fxsM(chIdx,:,:)),2) mean(squeeze(fxsF(chIdx,:,:)),2) mean(squeeze(hcsM(chIdx,:,:)),2) mean(squeeze(hcsF(chIdx,:,:)),2)]);
linePlotHandles(1).Color = [1 0 0 0.5];
linePlotHandles(1).LineWidth = 1;
linePlotHandles(1).LineStyle = '-';
linePlotHandles(2).Color = [1 0 0 0.5];
linePlotHandles(2).LineWidth = 1;
linePlotHandles(2).LineStyle = ':';
linePlotHandles(3).Color = [0 0 1 0.5];
linePlotHandles(3).LineWidth = 1;
linePlotHandles(3).LineStyle = '-';
linePlotHandles(4).Color = [0 0 1 0.5];
linePlotHandles(4).LineWidth = 1;
linePlotHandles(4).LineStyle = ':';
xlim([2900 3400])
ylim([-2.2 2.2])
line([0 0], ylim, 'color', [0 0 0])
line([3000 3000], ylim, 'color', [0 0 0], 'linestyle', ':')
line([3100 3100], ylim, 'color', [0.8 0.8 0.8])
xlabel('Latency (ms)')
ylabel('Amplitude (\muV)')
legend({'FXS M' 'FXS F' 'HCS M' 'HCS F'}, 'Location', 'best')

nestedLayout3 = tiledlayout(mainLayout, 2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
nestedLayout3.Layout.Tile = 11;
nexttile(nestedLayout3, 1)
topoplot(mean(squeeze(fxsM(:,4100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout3, 2)
topoplot(mean(squeeze(fxsF(:,4100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout3, 3)
topoplot(mean(squeeze(hcsM(:,4100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])
nexttile(nestedLayout3, 4)
topoplot(mean(squeeze(hcsF(:,4100,:)),2), EEG.chanlocs, 'maplimits', [-3 3])

nexttile(mainLayout, 12)
meanForBar = [mean(squeeze(fxsM(55,4100,:))) mean(squeeze(fxsF(55,4100,:))) mean(squeeze(hcsM(55,4100,:))) mean(squeeze(hcsF(55,4100,:)))];
seForBar  = [std(squeeze(fxsM(55,4100,:)))/sqrt(size(fxsM,3))  std(squeeze(fxsF(55,4100,:)))/sqrt(size(fxsF,3))  std(squeeze(hcsM(55,4100,:)))/sqrt(size(hcsM,3))  std(squeeze(hcsF(55,4100,:)))/sqrt(size(hcsF,3))];
barHandle = bar(meanForBar);
hold on
errorbar(1:4, meanForBar, seForBar, [0 0 0 0], 'linestyle', 'none', 'color', [0 0 0])
ylim([-3 0])
barHandle.FaceColor = "flat";
barHandle.CData     = [1 0 0; 1 0 0; 0 0 1; 0 0 1];
barHandle.FaceAlpha = 0.2;
ylabel('Amplitude (\muV)')
xticklabels({'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})

[stats3, df3, pvals3] = statcond({squeeze(fxsM(55,4100,:))' squeeze(fxsF(55,4100,:))'; squeeze(hcsM(55,4100,:))' squeeze(hcsF(55,4100,:))'});
    % tStats: {[0.1496]}    {[17.7403]}    {[0.4489]}
    % df: {[1 88]}    {[1 88]}    {[1 88]}
    % p:  {[0.6999]}    {[6.1035e-05]}    {[0.5046]}

[stats3_MaleFemale, stats3_FxsHcs] = cohensD_mainEffects2x2(squeeze(fxsM(55,4100,:))', squeeze(fxsF(55,4100,:))', squeeze(hcsM(55,4100,:))', squeeze(hcsF(55,4100,:))');
stats3_F2                          = cohensF2(stats3{3}, df3{3}(1), df3{3}(2));
disp(sprintf('%.3f %.3f %.3f', stats3_MaleFemale, stats3_FxsHcs, stats3_F2))
    % Cohen 0.081 0.711 0.005

exportgraphics(gcf, '/srv/Makoto/ASSR/p1400_ERP_maleFemale/ERPs.pdf', 'ContentType', 'vector', 'Resolution', 300, 'Colorspace', 'rgb')