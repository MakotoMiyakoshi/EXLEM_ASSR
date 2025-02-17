% 08/06/2024 Makoto. Lauren Ethridge wrote me this for homework.
%
%    "Cz referenced values from the fronto-central area (where we normally pull ITC),
%    which shouldn"t include any influence from temporal electrodes, show expected
%    ITC values consistent with what we%ve found previously in other
%    (cleaned, average referenced) datasets." (In her email on 7/31/2024)
%
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
scal = hztoscale(wtFreqBins, 'morlet', 1000);

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
uniqueEventMatrix = cell(length(allSets), 2);
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Re-reference to Cz. %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    E55Idx = find(strcmp({EEG.chanlocs.labels}, 'E55'));
    EEG = pop_reref(EEG, E55Idx);



    % Epoch the data.
    EEG = pop_epoch( EEG, {'DI66' 'DIN3' 'DIN7' 'DIN8'}, [-1  4], 'newname', 'EGI file epochs', 'epochinfo', 'yes');

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Electrode data. %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    elecData2D = EEG.data(:,:);
    elecErspMean   = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecErspMedian = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecItcMean    = single(zeros(EEG.nbchan, 100, EEG.pnts));
    elecItcMedian  = single(zeros(EEG.nbchan, 100, EEG.pnts));
    for elecIdx = 1:EEG.nbchan

        disp(sprintf('Electrode %d/%d...', elecIdx, EEG.nbchan))

        % Wavelet transformation using the conventional Morlet.
        coeffs = [fcwt(elecData2D(elecIdx,:)', scal, 'morlet')]';
        coeffs3D = reshape(coeffs, [size(coeffs,1) size(EEG.data, 2) size(EEG.data, 3)]);

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

    %%%%%%%%%%%%%%%%%%
    %%% Save data. %%%
    %%%%%%%%%%%%%%%%%%
    save(sprintf('/srv/Makoto/ASSR/p0220_epoch_CzRef/%s_elecErspMean',   dataName), 'elecErspMean',   '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0220_epoch_CzRef/%s_elecErspMedian', dataName), 'elecErspMedian', '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0220_epoch_CzRef/%s_elecItcMean',    dataName), 'elecItcMean',    '-nocompression')
    save(sprintf('/srv/Makoto/ASSR/p0220_epoch_CzRef/%s_elecItcMedian',  dataName), 'elecItcMedian',  '-nocompression')
end