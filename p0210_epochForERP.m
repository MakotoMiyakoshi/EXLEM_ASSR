% 01/28/2025 Makoto. Used again.
% 12/11/2024 Makoto. Used again.
% 12/02/2024 Makoto. Used again.
% 11/29/2024 Makoto. Used again.
% 07/12/2024 Makoto. Created.

close all
clear
clc

% allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
% for setIdx = 1:length(allSets)
% 
%     loadName = allSets(setIdx).name;
%     dataName = loadName(1:4);
%     loadPath = allSets(setIdx).folder;
%     EEG = pop_loadset('filename', loadName, 'filepath', loadPath);
% 
%     % Epoch the data.
%     EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1  4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');
% 
%     % Save the data.
%     pop_saveset(EEG, 'filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0210_epochForERP');
% end


%% Compute grand-mean and -median ERP.

allSets = dir('/srv/Makoto/ASSR/p0210_epochForERP/*.set');
meanERP   = single(zeros(128, 5000, length(allSets)));
medianERP = single(zeros(128, 5000, length(allSets)));
numTrialsVec = zeros(length(allSets),1);
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    numTrialsVec(setIdx) = EEG.trials;
    meanERP(  :,:,setIdx) = mean(  EEG.data,3);
    medianERP(:,:,setIdx) = median(EEG.data,3);
end

save('/srv/Makoto/ASSR/p0210_epochForERP/meanMedianERPs', 'meanERP', 'medianERP')
save('/srv/Makoto/ASSR/p0210_epochForERP/numTrialsVec', 'numTrialsVec')

figure
plot(EEG.times, mean(meanERP(  55,:,:),3)); hold on
plot(EEG.times, mean(medianERP(55,:,:),3))