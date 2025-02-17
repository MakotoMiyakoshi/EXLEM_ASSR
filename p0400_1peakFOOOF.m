% 07/11/2024 Makoto. Modified.
% 06/24/2024 makoto. Modified.
% 05/30/2024 Makoto. Showing three results: 0, 1, and 3 peak models.
% 05/29/2024 Makoto. Fitting 3 Gaussians allowing broad speaks is not optimal for evaluating 1/f.
% 01/10/2024 Makoto. Error detected.
% 01/09/2024 Makoto. Used.
% 01/05/2024 Makoto. Created.
clear
close all
clc

addpath(genpath('/srv/Makoto/Ribbon/code/ribbonanalysistoolbox'))
allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
for setIdx = 1:2:length(allSets)

    % Load data.
    currentSet  = allSets(setIdx).name;
    currentPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', currentSet, 'filepath', currentPath);

    %%%%%%%%%%%%%
    %%% FOOOF %%%
    %%%%%%%%%%%%%
    % Calculate and store PSD.
    EEG.icaact = EEG.icaweights*EEG.icasphere*double(EEG.data);
    [spectra,psdFreqs] = spectopo(EEG.icaact, EEG.pnts, EEG.srate, 'freqfac', 10, 'overlap', EEG.srate/2, 'plot', 'off');

    % Get the valid freq bin indices
    goodFreqIdx = find(psdFreqs>=1 & psdFreqs <= 55);
    goodFreqsHz = double(psdFreqs(goodFreqIdx));

    fooofOutputs = cell(size(spectra,1),1);
    for icIdx = 1:size(spectra,1)

        currentPSD = spectra(icIdx,:);

        % Trim and Mask the psdMatrix.
        psdTrimmed = double(currentPSD(goodFreqIdx));

        %       Iteratively fit peaks to flattened spectrum.
        %
        %       Parameters
        %       ----------
        %       freqs : 1xn array
        %           Frequency values for the power spectrum, in linear scale.
        %       flat_iter : 1xn array
        %           Flattened (aperiodic removed) power spectrum.
        %       max_n_peaks : double
        %           Maximum number of gaussians to fit within the spectrum.
        %       peak_threshold : double
        %           Threshold (in standard deviations of noise floor) to detect a peak.
        %       min_peak_height : double
        %           Minimum height of a peak (in log10).
        %       gauss_std_limits : 1x2 double
        %           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
        %       proxThresh : double
        %           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
        %       peakType : {'gaussian', 'cauchy', 'both'}
        %           Which types of peaks are being fitted
        %       guess_weight : {'none', 'weak', 'strong'}
        %           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
        %       hOT : 0 or 1
        %           Defines whether to use constrained optimization, fmincon, or
        %           basic simplex, fminsearch.
        %
        %       Returns
        %       -------
        %       gaussian_params : mx3 array, where m = No. of peaks.
        %           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].

        opt.freq_range          = [goodFreqsHz(1) goodFreqsHz(end)];
        opt.aperiodic_mode      = 'fixed';   % aperiodic_mode : {'fixed', 'knee'} Defines absence or presence of knee in aperiodic component.
        opt.max_peaks           = 1;        % Maximum number of gaussians to fit within the spectrum.
        opt.peak_threshold      = 1;        % 2 std dev: parameter for interface simplification
        opt.min_peak_height     = 0;        % Minimum height of a peak (in log10).
        opt.peak_width_limits   = [0.5 3];    % Values were taken from https://neuroimage.usc.edu/brainstorm/Tutorials/Fooof
        opt.proximity_threshold = 2;        % Minimum distance between two peaks, in st. dev. (gamma) of peaks.
        %opt.peak_type           = 'cauchy'; % {'gaussian', 'cauchy', 'both'}
        opt.peak_type           = 'best'; % {'gaussian', 'cauchy', 'both'}
        opt.guess_weight        = 'weak';   % Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
        hOT                     = 1;        % Defines whether to use constrained optimization, fmincon, or basic simplex, fminsearch.
        opt.thresh_after        = true;     % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)

        % Convert dB to microVolt^2 (the FOOOF requirement)
        psdMatrixTrimmed_microVoltSquare = 10.^(psdTrimmed/10);

        %%%%%%%%%%%%%%%%%%
        %%% Run FOOOF. %%%
        %%%%%%%%%%%%%%%%%%
        [fs, fg] = FOOOF_matlab(psdMatrixTrimmed_microVoltSquare, goodFreqsHz', opt, hOT); % currentIcPSD MUST be a row vector.
        fg.fs    = fs;

        % Store FOOOF results.
        fooofOutputs{icIdx} = fg;
    end
    EEG.etc.FOOOF_concatenated       = fooofOutputs;
    EEG.etc.PSD_concatenated.spectra = spectra;
    EEG.etc.PSD_concatenated.freqs   = psdFreqs;

    % Save the data.
    pop_saveset(EEG, 'filename', currentSet, 'filepath', '/srv/Makoto/ASSR/p0400_1peakFOOOF')
end