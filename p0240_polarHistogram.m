% 08/08/2024 Makoto. Created.
% 07/12/2024 Makoto. itcMedian was mean. Fixed.
% 07/08/2024 Makoto. Created.
close all
clear
clc

addpath('/srv/Makoto/Tools/CircStat2012a')

% Prepare video output parameters.
outputVideo = VideoWriter('/srv/Makoto/ASSR/p0240_polarHistgram/singleSubj40HzAssrItc_polarHist', 'Motion JPEG AVI');
outputVideo.FrameRate = 30;
outputVideo.Quality   = 95;
open(outputVideo)

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
for setIdx = 30 %1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

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
    for elecIdx = 6 %1:EEG.nbchan

        disp(sprintf('Electrode %d/%d...', elecIdx, EEG.nbchan))

        % Wavelet transformation using the conventional Morlet.
        coeffs = [fcwt(elecData2D(elecIdx,:)', scal, 'morlet')]';
        coeffs3D = reshape(coeffs, [size(coeffs,1) size(EEG.data, 2) size(EEG.data, 3)]);

        inputData = squeeze(coeffs3D(freqIdx7,1001:4000,:));
        instPhase = angle(inputData);

        normalizedData = coeffs3D./sqrt(abs(coeffs3D).^2);
        itc            = squeeze(abs(mean(normalizedData(freqIdx7,1001:4000,:),3)));

        figHandle = figure;
        for timeIdx = 1:1000
            polarhistogram(instPhase(timeIdx,:), 18, 'Normalization', 'count')
            title(sprintf('Subj2287 FCz at 40Hz, %d ms, ITC=%.2f, Angle %.0f (SD %.0f)', timeIdx, itc(timeIdx), ...
                circ_rad2ang(circ_mean([instPhase(timeIdx,:)]')), ...
                circ_rad2ang(circ_std( [instPhase(timeIdx,:)]'))))

            % Write the current figure to a frame of the video.
            currentMovieFrame = getframe(figHandle);
            writeVideo(outputVideo, currentMovieFrame);

            cla
        end
    end
end

close(outputVideo)