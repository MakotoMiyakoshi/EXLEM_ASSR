% 09/06/2024 Makoto. Concluded that rejecting 20% of trials made only trivial differences.
% 07/12/2024 Makoto. itcMedian was mean. Fixed.
% 07/08/2024 Makoto. Created.
close all
clear
clc

% Set the rejection rate.
rejRate = 0.2;

% Set wavelet parameters.
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

% Compute Siyi Deng's fast Morlet wavelet told by Ramesh.
scal = hztoscale(wtFreqBins, 'morlet', 1000);

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
uniqueEventMatrix = cell(length(allSets), 2);
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);
    EEGorig = EEG;

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Reject bad epochs. %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 40-Hz band-pass filter.
    EEG = pop_firws(EEG, 'fcutoff', [37.5 42.5], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', 660, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Epoch the data.
    EEG = pop_epoch(EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1 4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Obtain the epoch quality index for each ch (larger, better).
    sortIdxChanMatrix = zeros(EEG.trials, EEG.nbchan);
    assrIdx = find(EEG.times>=0 & EEG.times<=3000);
    for chIdx = 1:EEG.nbchan
        data2D = double(squeeze(EEG.data(chIdx,assrIdx,:)))';
        [V,D] = eigs(data2D*data2D', size(data2D, 1));
        firstEigenVector = V(:,1);
            %{
            figure
            bar(sort(firstEigenVector))
            %}
        [~, sortingIdx] = sort(firstEigenVector, 'ascend');
        sortIdxChanMatrix(:,chIdx) = sortingIdx;
    end


    %%

    % Epoch the data.
    EEG = EEGorig;
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1  4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Electrode data. %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    elecData2D = double(EEG.data(:,:));
    elecErspMean   = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecErspMedian = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecItcMean    = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecItcMedian  = single(zeros(EEG.nbchan, 100, EEG.pnts));
    for chIdx = 1:EEG.nbchan

        disp(sprintf('Electrode %d/%d...', chIdx, EEG.nbchan))

        % Wavelet transformation using the conventional Morlet.
        coeffs = [fcwt(elecData2D(chIdx,:)', scal, 'morlet')]';
        coeffs3D = reshape(coeffs, [size(coeffs,1) size(EEG.data, 2) size(EEG.data, 3)]);

        % Reject the worst epochs.
        epochRejIdx = sortIdxChanMatrix(1:round(EEG.trials*rejRate), chIdx);
        coeffs3D(:,:,epochRejIdx) = [];

        % Calculate ERSP.
        erspMean   = mean(  abs(coeffs3D).^2, 3);
        erspMedian = median(abs(coeffs3D).^2, 3);

        % Calculate ITPC.
        normalizedData = coeffs3D./sqrt(abs(coeffs3D).^2);
        itcMean     = abs(mean(normalizedData,3));
        median_real = median(real(normalizedData),3);
        median_imag = median(imag(normalizedData),3);
        itcMedian   = abs(median_real+1i*median_imag);
            %{
            figure
            subplot(1,2,1)
            imagesc(itcMean); axis xy; colormap jet
            subplot(1,2,2)
            imagesc(itcMedian); axis xy; colormap jet
            %}

        % Store the results.
        elecErspMean(  chIdx,:,:) = erspMean;
        elecErspMedian(chIdx,:,:) = erspMedian;
        elecItcMean(   chIdx,:,:) = itcMean;
        elecItcMedian( chIdx,:,:) = itcMedian;
    end


    % Save data.
    save(sprintf('/srv/Makoto/ASSR/p0230_epoch_EigRej/%s_elecErspMean',   dataName), 'elecErspMean',   '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0230_epoch_EigRej/%s_elecErspMedian', dataName), 'elecErspMedian', '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0230_epoch_EigRej/%s_elecItcMean',    dataName), 'elecItcMean',    '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0230_epoch_EigRej/%s_elecItcMedian',  dataName), 'elecItcMedian',  '-nocompression')
end