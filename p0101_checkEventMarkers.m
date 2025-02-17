% 11/29/2024 Makoto. Used again.
% 07/08/2024 Makoto. Created.
close all
clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wavelet on electrode data to extract just 40 Hz data. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/srv/Makoto/Tools/siyisCodeFromRamesh')
scal = hztoscale(40,'morlet', 1000);

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
uniqueEventMatrix = cell(length(allSets), 2);
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath);

    uniqueEvents = unique({EEG.event.type});
    
    % Count how many events are there.
    numEvents        = zeros(length(uniqueEvents),1);
    eventLatencyCell = cell(length(uniqueEvents),1);
    for uniqueEventIdx = 1:length(uniqueEvents)
        currentEventIdx = find(strcmp(uniqueEvents{uniqueEventIdx}, {EEG.event.type}));
        eventLatencyCell{uniqueEventIdx,1} = [EEG.event(currentEventIdx).latency];
        numEvents(uniqueEventIdx) = sum(strcmp(uniqueEvents{uniqueEventIdx}, {EEG.event.type}));
    end

    % Remove <50 events.
    excludingEventIdx = find(numEvents<50);
    uniqueEvents(    excludingEventIdx) = [];
    eventLatencyCell(excludingEventIdx) = [];
    numEvents(       excludingEventIdx) = [];

    % Determine which event is trigger onset.
    czItcComparison  = zeros(1,2);
    czErspComparison = zeros(1,2);
    for eventTypeIdx = 1:2

        % Find time of the events.
        currentEventTypeOnsetLatencyVec = [eventLatencyCell{eventTypeIdx}]';
        onsetOffsetMatrix = [currentEventTypeOnsetLatencyVec-EEG.srate currentEventTypeOnsetLatencyVec+2*EEG.srate]; 

        if any(find(onsetOffsetMatrix<0))
            excludingIdx = find(onsetOffsetMatrix(:,1)<0);
            onsetOffsetMatrix(excludingIdx,:) = [];
        end
        if any(find(onsetOffsetMatrix>EEG.pnts))
            excludingIdx = find(onsetOffsetMatrix(:,2)>EEG.pnts);
            onsetOffsetMatrix(excludingIdx,:) = [];
        end

        % Calcuate Cz time-frequency data.
        coeffs = [fcwt(EEG.data(6,:)', scal, 'morlet')]';

        erspStack = zeros(1, 3001, size(onsetOffsetMatrix,1));
        for epochIdx = 1:size(onsetOffsetMatrix,1)
            erspStack(:,:,epochIdx) = coeffs(:, onsetOffsetMatrix(epochIdx,1):onsetOffsetMatrix(epochIdx,2));
        end

        % ITC.
        normalizedData = erspStack./sqrt(abs(erspStack).^2);
        median_real = median(real(normalizedData),3);
        median_imag = median(imag(normalizedData),3);
        itc         = abs(median_real+1i*median_imag);
        %itc = abs(mean(erspStack./sqrt(abs(erspStack).^2), 3));
            %{
            figure
            plot(itc)
            %}

        % ERSP.
        medianErsp            = median(abs(erspStack).^2, 3);
        baselineCorrectedErsp = bsxfun(@minus, medianErsp, mean(medianErsp(:, 1:1000), 2));
            %{
            figure
            imagesc(baselineCorrectedErsp)
            %}

        czItcComparison( eventTypeIdx) = mean(itc(1, 1001:end))-mean(itc(1, 1:1000));
        czErspComparison(eventTypeIdx) = mean(baselineCorrectedErsp(1, 1001:end))-mean(baselineCorrectedErsp(1, 1:1000));
    end
    EEG.etc.uniqueEvents     = uniqueEvents;
    EEG.etc.czItcComparison  = czItcComparison;
    EEG.etc.czErspComparison = czErspComparison;

    % Save the data.
    pop_saveset(EEG, 'filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0100_upToDipfit');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if all the stimulus onset types are correct. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear

allSets = dir('/srv/Makoto/ASSR/p0100_upToDipfit/*.set');
stimOnsetType = cell(length(allSets), 1);
difference    = nan(length(allSets), 2);
for setIdx = 1:length(allSets)

    % if setIdx == 5 | setIdx == 8 | setIdx == 20 | setIdx == 21 | setIdx == 40
    %     continue
    % end

    loadName = allSets(setIdx).name;
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath, 'loadmode', 'info');

    [~, itcMaxIdx]  = max(EEG.etc.czItcComparison);
    [~, erspMaxIdx] = max(EEG.etc.czErspComparison);
    if itcMaxIdx~=erspMaxIdx
        error('ERSP and ITC results do not match.')
    end
    stimOnsetType{setIdx} = EEG.etc.uniqueEvents{itcMaxIdx};
    difference(setIdx,:)  = [max(EEG.etc.czItcComparison) max(EEG.etc.czErspComparison)];
end

stimOnsetTypeNoSpace = stimOnsetType;
emptyCellIdx = find(cellfun(@isempty, stimOnsetTypeNoSpace));
stimOnsetTypeNoSpace(emptyCellIdx) = [];
unique(stimOnsetTypeNoSpace)

    % {'DI66'}
    % {'DIN3'}
    % {'DIN7'}
    % {'DIN8'}