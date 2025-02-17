% 07/12/2024 Makoto. Created.

close all
clear
clc

allSets = dir('/srv/Makoto/ASSR/p0410_selectBrainClassICs/brainClass/*.set');
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    % Epoch the data.
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1  4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Save the data.
    pop_saveset(EEG, 'filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0440_epochBrainForERP');
end


%% Compute grand-mean and -median ERP.

allSets = dir('/srv/Makoto/ASSR/p0440_epochBrainForERP/*.set');
meanERP   = single(zeros(128, 5000, length(allSets)));
medianERP = single(zeros(128, 5000, length(allSets)));
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    meanERP(  :,:,setIdx) = mean(  EEG.data,3);
    medianERP(:,:,setIdx) = median(EEG.data,3);
end

save('/srv/Makoto/ASSR/p0440_epochBrainForERP/meanMedianERPs_brainClass', 'meanERP', 'medianERP', '-nocompression')

figure
plot(mean(meanERP(55,:,:),3));