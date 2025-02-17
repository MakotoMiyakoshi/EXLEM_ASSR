% 08/05/2024 Makoto. Modified.
% 07/12/2024 Makoto. Created.

close all
clear
clc

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    % 40-Hz band-pass filter.
    EEG = pop_firws(EEG, 'fcutoff', [37.5 42.5], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', 660, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Epoch the data.
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [0.5  3], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Nakayoshi special.
    data2D = permute(EEG.data, [3 1 2]);
    data2D = double(data2D(:,:));
    [V,D] = eigs(data2D*data2D', size(data2D, 1));
    firstEigenVector = V(:,1);
    numPositiveEpochs = sum(firstEigenVector>0);
    if numPositiveEpochs>=50
        [~, sortingIdx] = sort(firstEigenVector, 'ascend');
    else
        [~, sortingIdx] = sort(firstEigenVector, 'descend');
    end
    exclusionIdx = sortingIdx(1:round(length(sortingIdx)*0.2));

    % Re-load the data.
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1  4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Reject epochs.
    EEG = pop_select(EEG,'notrial', exclusionIdx);

        %{
        figure
        czData = mean(EEG.data(55,:,:),3);
        plot(czData);
        hold on
        czData3D_clean = squeeze(EEG.data(55,:,:));
        czData3D_clean(:,exclusionIdx) = [];
        plot(mean(czData3D_clean,2))
        %}

        % % Trial selection.
        % medianError = bsxfun(@minus, EEG.data, median(EEG.data,3));
        % figure
        % imagesc(squeeze(mean(medianError,2)))
        % 
        % figure
        % hist(medianError(:), 100)
        % meanError = mean(medianError(:));
        % stdError  = std(medianError(:));
        % 
        % outlierMask = medianError<meanError-stdError*2 | medianError>meanError+stdError*2;
        % 
        % sum(outlierMask(:))/(128*5000*100)
        % 
        % figure
        % czData = mean(EEG.data(55,:,:),3);
        % plot(czData);
        % hold on
        % czData3D = squeeze(EEG.data(55,:,:));
        % czMask   = squeeze(outlierMask(55,:,:));
        % czMasked = czData3D;
        % czMasked(czMask) = NaN;
        % plot(nanmean(czMasked,2))
        % 
        % figure
        % imagesc(czMasked)

    % Save the data.
    pop_saveset(EEG, 'filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0600_epochForERP_trialRej');
end


%% Compute grand-mean and -median ERP.

allSets = dir('/srv/Makoto/ASSR/p0600_epochForERP_trialRej/*.set');
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

save('/srv/Makoto/ASSR/p0600_epochForERP_trialRej', 'meanERP', 'medianERP', '-nocompression')

figure
plot(mean(meanERP(55,:,:),3));