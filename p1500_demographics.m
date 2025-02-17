% 01/07/2024 Makoto. Created.

% Load the latest demographic data.
demoTable = readtable('/srv/Makoto/ASSR/p0000_raw/Final_List_SSCT_Adults_Adolescents_01062025.xlsx');
eegIdVec = demoTable.eegid;

% Exclusion list by Craig.
exclusionList = [44 532 554 879 2537];
[C, IA] = setdiff(eegIdVec, exclusionList);
demoTable = demoTable(IA, :);

% Obtain FMRP.
fmrpVec = demoTable.fmrp_replicate_trim_mean;
fmrpLogVec = log10(fmrpVec+1);

% FXS, HCS, M, F.
fxsmIdx = find(strcmp(demoTable.Sex, 'M') & demoTable.Dx_binary==1);
fxsfIdx = find(strcmp(demoTable.Sex, 'F') & demoTable.Dx_binary==1);
hcsmIdx = find(strcmp(demoTable.Sex, 'M') & demoTable.Dx_binary==0);
hcsfIdx = find(strcmp(demoTable.Sex, 'F') & demoTable.Dx_binary==0);

% Plot.
figure('position', [100 100 800 400])
subplot(1,2,1)
bar([nanmean(fmrpLogVec(fxsmIdx)) nanmean(fmrpLogVec(fxsfIdx)) nanmean(fmrpLogVec(hcsmIdx)) nanmean(fmrpLogVec(hcsfIdx))])
hold on
errorbar(1:4, [nanmean(fmrpLogVec(fxsmIdx)) nanmean(fmrpLogVec(fxsfIdx)) nanmean(fmrpLogVec(hcsmIdx)) nanmean(fmrpLogVec(hcsfIdx))], ...
         [0 0 0 0], [nanstd(fmrpLogVec(fxsmIdx)) nanstd(fmrpLogVec(fxsfIdx)) nanstd(fmrpLogVec(hcsmIdx)) nanstd(fmrpLogVec(hcsfIdx))], ...
         'LineStyle', 'none', 'Color', [0 0 0])
title('FMRP, log-transformed')
ylabel('log10(picomol)')
set(gca,'XTickLabel', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})

subplot(1,2,2)
bar([nanmean(fmrpVec(fxsmIdx)) nanmean(fmrpVec(fxsfIdx)) nanmean(fmrpVec(hcsmIdx)) nanmean(fmrpVec(hcsfIdx))])
hold on
errorbar(1:4, [nanmean(fmrpVec(fxsmIdx)) nanmean(fmrpVec(fxsfIdx)) nanmean(fmrpVec(hcsmIdx)) nanmean(fmrpVec(hcsfIdx))], ...
         [0 0 0 0], [nanstd(fmrpVec(fxsmIdx)) nanstd(fmrpVec(fxsfIdx)) nanstd(fmrpVec(hcsmIdx)) nanstd(fmrpVec(hcsfIdx))], ...
         'LineStyle', 'none', 'Color', [0 0 0])
title('FMRP')
ylabel('picomol')
set(gca,'XTickLabel', {'FXS-M' 'FXS-F' 'HCS-M' 'HCS-F'})

print('/srv/Makoto/ASSR/p1500_demographics', '-djpeg95', '-r200')