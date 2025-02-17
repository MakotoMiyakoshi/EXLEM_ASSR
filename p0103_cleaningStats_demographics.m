% 11/29/2024 Makoto. Used again.
% 09/06/2024 Makoto. Created.
clear
clc

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

figure
bar(sort(ageVec))
title('Age (sorted)')

tooYoungIdx = find(ageVec<6);
subjNames( tooYoungIdx)'
groupNames(tooYoungIdx)'

earlyAdolescenceIdx = find(ageVec>=6 & ageVec<=12);
subjNames( earlyAdolescenceIdx)'
groupNames(earlyAdolescenceIdx)'

allAbove13Idx = find(ageVec>12);
subjNames_13plus  = subjNames( allAbove13Idx)';
groupNames_13plus = groupNames(allAbove13Idx)';

sum(contains(groupNames_13plus, 'FXS')) % 35 FXS, 33 HCS.
sum(contains(sexGroupList(:,2), 'TDC'))


mean(cleaningMatrix(:,[1 2])) % Ch rej, 0.0441 (0.2069) [0-1]; dB reduce, -2.57 (1.86), [-8.05 -0.13]. % Old record: Ch rej, 2.4 (1.84) [0-7]; dB reduce, -2.8 (2.1), [-7.8 - -0.1].
std(cleaningMatrix(:,[1 2]))
min(cleaningMatrix(:,[1 2]))
max(cleaningMatrix(:,[1 2]))

