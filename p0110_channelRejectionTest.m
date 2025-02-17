% 11/29/2024 Makoto. Used again.
% 07/08/2024 Makoto. 0558 and 1774 error.
% 07/06/2024 Makoto. Used. Originally for the Habituation data.
% 06/18/2024 Makoto. Modified.
% 06/07/2024 Makoto. Created.

clear global
close all
clc

datasetTable = readtable('/srv/Makoto/ASSR/p0000_raw/Final_SSCT_Adults_Adolescents_editedLAD_11152024_modified.xlsx');

% Find subjects with >= 13 years old. Looks like everyone is above 13!? (11/27/2024)
allListedNamesDouble = datasetTable.Var1;
ageListCell = datasetTable.Var3;
yearMonthMatrix = cell2mat(cellfun(@(x) sscanf(x, '%d years , %d months')', ageListCell, 'UniformOutput', false));
yearVec = yearMonthMatrix(:,1)+yearMonthMatrix(:,2)/12;

% Move all datasets that are not included in this list.
allRaws = dir('/srv/Makoto/ASSR/p0000_raw/*.raw'); % EGI .raw format?
allRawNames = cellfun(@(x) x(1:4), {allRaws.name}, 'uniformoutput', false)';
allRawNamesDouble = str2num(cell2mat(allRawNames));
% [exclusionDataID,IA] = setdiff(allRawNamesDouble, allListedNamesDouble);
% for excludeIdx = 1:length(excludedDataID)
%     currentTarget = exclusionDataID(excludeIdx);
%     allRaws = dir('/srv/Makoto/ASSR/p0000_raw/*.raw'); % EGI .raw format?
%     allRawNames = cellfun(@(x) x(1:4), {allRaws.name}, 'uniformoutput', false)';
%     allRawNamesDouble = str2num(cell2mat(allRawNames));
%     exclusionIdx = find(allRawNamesDouble==currentTarget);
%     targetName = allRaws(exclusionIdx).name;
%     SOURCE      = [allRaws(exclusionIdx).folder filesep allRaws(exclusionIdx).name];
%     DESTINATION = [allRaws(exclusionIdx).folder filesep 'excluded' filesep allRaws(exclusionIdx).name];
%     movefile(SOURCE, DESTINATION,'f');
% end


%%

numExcludedChVec = zeros(length(allRaws),1);
for rawIdx = [62] %1:length(allRaws)

    % Import EGI data.
    loadName = allRaws(rawIdx).name;
    dataName = loadName(1:4);
    EEG = pop_readegi(['/srv/Makoto/ASSR/p0000_raw' filesep loadName], [],[],'auto');

    % Populate demographic data.
    excelTableIdx = find(datasetTable.Var1 == str2num(dataName));
    EEG.etc.AgeAtVisit = datasetTable.Var3(excelTableIdx);
    EEG.etc.Sex  = datasetTable.Var4(excelTableIdx);
    EEG.etc.Dx   = datasetTable.Var6(excelTableIdx);
    EEG.etc.MosaicStatus = datasetTable.Var8(excelTableIdx);

    % Copy some of the info to EEG structure.
    EEG.subject = dataName;
    EEG.group   = EEG.etc.Dx;

    % Import channel data.
    EEG = pop_chanedit(EEG, 'load',{'/srv/Makoto/Habituation/code/GSN-HydroCel-128.sfp','filetype','autodetect'});

    % % Resample the data to 250 Hz.
    EEG = pop_resample(EEG, 100);

    % Apply 0.5-Hz HPF (TBW 1.0)
    EEG = pop_firws(EEG, 'fcutoff', 0.5, 'ftype', 'highpass', 'wtype', 'blackman', 'forder', 550, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

    % Remove line noise up to Nyquist.
    %EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist', 1:EEG.nbchan,'computepower',1,'linefreqs',60:60:480,'newversion',1,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

    % Discard 3 second of data from both edges.
    EEG = pop_select( EEG, 'time', [3 EEG.xmax-3]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Detect outlier channels (experimental, simple) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Re-reference to median.
    medianPotential = median(EEG.data,1);
    medianRerefData = bsxfun(@minus, EEG.data, medianPotential);



    %%%%%%%%%%%%
    %%% ASR. %%%
    %%%%%%%%%%%%
    EEG_beforeASR = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off', ...
        'LineNoiseCriterion','off','Highpass','off', ...
        'BurstCriterion',25, 'WindowCriterion', 0.75, ...
        'BurstRejection','off','Distance','Euclidian', ...
        'BurstCriterionRefMaxBadChns', 0, ...
        'BurstCriterionRefTolerances', [-inf 8], ...
        'WindowCriterionTolerances',[-Inf 8]);
    %{
        vis_artifacts(EEG, EEG_beforeASR);
    %}








    % Segment the data into 1000 epochs.
    permillWindowLength = floor(EEG.pnts/1000);
    trimmingLength      = mod(EEG.pnts, permillWindowLength);
    trimmedData         = EEG.data(:,1:end-trimmingLength);
    numEpochs           = size(trimmedData,2)/permillWindowLength;
    trimmedData3D       = reshape(trimmedData, [size(trimmedData,1) permillWindowLength numEpochs]);

    % Calculate std for each.
    std3D = squeeze(std(trimmedData3D, 0, 2));
    %{
        figure; imagesc(std3D, [0 100])
    
        figure
        tiledlayout('flow')
        nexttile
        hist(std3D(1,:), 100)
        nexttile
        hist(std3D(126,:), 100)
    
        figure
        bar(mean(std3D,2))
    %}

    % Calculate robust 30 SD channel rejection threshold for large amplitude.
    channelRejectionThreshold = median(median(std3D,2))+30*mad(median(std3D,2),1)*1.4826;

    %channelRejectionThreshold = median(mean(std3D,2))+20*mad(mean(std3D,2),1)*1.4826;

    largeAmplitudeChannelIdx  = find(mean(std3D,2)>channelRejectionThreshold);

    % Find flat channels.
    flatChannelIdx = find(mean(std3D,2)<median(mean(std3D,2))/10);

    % Combine channel indices to exclude.
    excludingChannelIdx = [flatChannelIdx; largeAmplitudeChannelIdx];

    numExcludedChVec(rawIdx) = length(excludingChannelIdx);
end

figure
bar(numExcludedChVec)