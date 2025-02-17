% 01/10/2025 Makoto. Created.
clear
close all

resultTable = readtable('/srv/Makoto/ASSR/code/statsResults002.xlsx');

column1 = resultTable.Var1;
pIdx = find(strcmp(column1, 'p'));
dIdx = find(strcmp(column1, 'd'));
pVals = resultTable.FXS_HCS(pIdx);
measureLabels = column1(pIdx-2);
fdrP = fdr(pVals(2:end));

ERP_ERSP_ITC_PowerERSSP_ITCERSSP = measureLabels(2:end);
uncorrP = pVals(2:end);
d       = resultTable.FXS_HCS(dIdx(2:end));

table(ERP_ERSP_ITC_PowerERSSP_ITCERSSP, uncorrP, fdrP, d)