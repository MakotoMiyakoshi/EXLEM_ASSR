% 11/30/2024 Makoto. Epoch after continuous WT.
% 11/29/2024 Makoto. Used again.
% 07/12/2024 Makoto. itcMedian was mean. Fixed.
% 07/08/2024 Makoto. Created.
close all
clear
clc

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
scal = hztoscale(wtFreqBins, 'morlet', 1000); % For 1000 Hz sampling.

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
uniqueEventMatrix = cell(length(allSets), 2);
for setIdx = length(allSets):-2:64

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    % Obtain latencies to epoch data.
    allEventTypes = {EEG.event.type}';
    onsetIdx     = find(contains(allEventTypes, 'DI66') | contains(allEventTypes, 'DIN3') | contains(allEventTypes, 'DIN7') | contains(allEventTypes, 'DIN8'));
    onsetLatency = [EEG.event(onsetIdx).latency]';
    edgeLatency  = [onsetLatency-1000 onsetLatency+4000]; % -1000 to 4000 ms.
    if any(edgeLatency(:,1)<1 | edgeLatency(:,2)>EEG.pnts)
        exclusionIdx = find(edgeLatency(:,1)<1 | edgeLatency(:,2)>EEG.pnts);
        edgeLatency(exclusionIdx,:) = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Electrode data. %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    elecData2D = EEG.data(:,:);
    elecErspMean   = single(zeros(EEG.nbchan, 100, length(-1000:4000)));
    elecErspMedian = single(zeros(EEG.nbchan, 100, length(-1000:4000)));
    elecItcMean    = single(zeros(EEG.nbchan, 100, length(-1000:4000)));
    elecItcMedian  = single(zeros(EEG.nbchan, 100, length(-1000:4000)));
    for elecIdx = 1:EEG.nbchan

        disp(sprintf('Electrode %d/%d...', elecIdx, EEG.nbchan))

        % Wavelet transformation using the conventional Morlet.
        coeffs = [fcwt(elecData2D(elecIdx,:)', scal, 'morlet')]';

        % Epoch the time-freq coefficients.
        coeffs3D = single(zeros(numFreqBins, length(-1000:4000), size(edgeLatency, 1)));
        for epochIdx = 1:size(coeffs3D,3)
            currentEpoch = coeffs(:, edgeLatency(epochIdx,1):edgeLatency(epochIdx,2));
            coeffs3D(:,:,epochIdx) = currentEpoch;
        end

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
        elecErspMean(  elecIdx,:,:) = erspMean;
        elecErspMedian(elecIdx,:,:) = erspMedian;
        elecItcMean(   elecIdx,:,:) = itcMean;
        elecItcMedian( elecIdx,:,:) = itcMedian;
    end


    %%%%%%%%%%%%%%%%
    %%% IC data. %%%
    %%%%%%%%%%%%%%%%
    icaact2D     = EEG.icaweights*EEG.icasphere*EEG.data;
    icErspMean   = single(zeros(size(EEG.icaweights,1), 100, length(-1000:4000)));
    icErspMedian = single(zeros(size(EEG.icaweights,1), 100, length(-1000:4000)));
    icItcMean    = single(zeros(size(EEG.icaweights,1), 100, length(-1000:4000)));
    icItcMedian  = single(zeros(size(EEG.icaweights,1), 100, length(-1000:4000)));
    for icIdx = 1:size(EEG.icaweights,1)

        disp(sprintf('IC %d/%d...', icIdx, size(EEG.icaweights,1)))

        % Wavelet transformation using the conventional Morlet.
        coeffs = [fcwt(icaact2D(icIdx,:)', scal, 'morlet')]';

        % Epoch the time-freq coefficients.
        coeffs3D = single(zeros(numFreqBins, length(-1000:4000), size(edgeLatency, 1)));
        for epochIdx = 1:size(coeffs3D,3)
            currentEpoch = coeffs(:, edgeLatency(epochIdx,1):edgeLatency(epochIdx,2));
            coeffs3D(:,:,epochIdx) = currentEpoch;
        end

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
            imagesc(erspMedian); axis xy; colormap jet
            subplot(1,2,2)
            imagesc(itcMedian); axis xy; colormap jet
            %}

        % Store the results.
        icErspMean(  icIdx,:,:) = erspMean;
        icErspMedian(icIdx,:,:) = erspMedian;
        icItcMean(   icIdx,:,:) = itcMean;
        icItcMedian( icIdx,:,:) = itcMedian;
    end


    %%%%%%%%%%%%%%%%%%
    %%% Save data. %%%
    %%%%%%%%%%%%%%%%%%
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_elecErspMean',   dataName), 'elecErspMean',   '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_elecErspMedian', dataName), 'elecErspMedian', '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_elecItcMean',    dataName), 'elecItcMean',    '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_elecItcMedian',  dataName), 'elecItcMedian',  '-nocompression')

    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_icErspMean',   dataName), 'icErspMean',   '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_icErspMedian', dataName), 'icErspMedian', '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_icItcMean',    dataName), 'icItcMean',    '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0200_epoch/%s_icItcMedian',  dataName), 'icItcMedian',  '-nocompression')
end