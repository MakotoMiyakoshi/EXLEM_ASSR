% 12/02/2024 Makoto. Used again.
% 10/31/2024 Makoto. Calculate delta-gamma power-power coupling.
% 07/12/2024 Makoto. itcMedian was mean. Fixed.
% 07/08/2024 Makoto. Created.
close all
clear
clc


%% Preselect datasets recorded from 13+ yrs old subjects.
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
    %{
    figure
    bar(sort(ageVec))
    %}

% Load dummy data.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0210_epochForERP/');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/p0000_raw/Final_SSCT_Adults_Adolescents_editedLAD_11152024_modified.xlsx');
groupIdx = demographicData.Var6;

% Sort the subjIDs.
[subjID_sorted, sortingIdx] = sort(demographicData.Var1);
groupIdx = groupIdx(sortingIdx);

hcsIdx = find(strcmp(groupIdx, 'TDC'));
fxsIdx = find(strcmp(groupIdx, 'FXS'));

% Obtain EEG IDs.
hcsNames = subjNames(hcsIdx);
fxsNames = subjNames(fxsIdx);


%% Set wavelet parameters.
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

% Set parameters for ASSR time bins, delta freq bins, and gamma freq bins. 
assrIdx  = 1001:4000;
deltaIdx = find(wtFreqBins>=1 & wtFreqBins<=2);
gammaIdx = freqIdx7;

% Compute Siyi Deng's fast Morlet wavelet told by Ramesh.
scal = hztoscale(wtFreqBins, 'morlet', 1000);

% Apply preselection for 13+ yrs olds.
allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
rMatrix = zeros(EEG.nbchan, length(allSets));
for setIdx = 1:length(allSets)

    disp(sprintf('%d/%d...', setIdx, length(allSets)))

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

        % Calculate delta-gamma power-power coupling.
        coeffs3D_power_log = 10*log10(coeffs3D.*conj(coeffs3D));
        coeffs3D_ERSP      = bsxfun(@minus, coeffs3D_power_log, mean(coeffs3D_power_log(:,1:1000,:), 2)) ;
        %{
        figure
        imagesc(mean(coeffs3D_ERSP,3))
        %}
        deltaTimeSeries = squeeze(mean(coeffs3D_ERSP(deltaIdx, assrIdx, :),1));
        gammaTimeSeries = squeeze(     coeffs3D_ERSP(gammaIdx, assrIdx, :));

        rVec = zeros(size(deltaTimeSeries,2),1);
        for epochIdx = 1:size(deltaTimeSeries,2)
            r = corrcoef([deltaTimeSeries(:,epochIdx) gammaTimeSeries(:,epochIdx)]);
            rVec(epochIdx) = r(2,1);
        end

        % Store the results.
        rMatrix(elecIdx, setIdx) = median(rVec);
    end
end

%%%%%%%%%%%%%%%%%%
%%% Save data. %%%
%%%%%%%%%%%%%%%%%%
save('/srv/Makoto/ASSR/p0900_deltaGammaCoupling/rMatrix', 'rMatrix', 'groupIdx', 'hcsIdx', 'fxsIdx', 'hcsNames', 'fxsNames',   '-nocompression')