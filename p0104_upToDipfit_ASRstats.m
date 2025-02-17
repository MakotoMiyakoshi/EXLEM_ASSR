% 01/28/2025 Makoto. Modified.
% 12/02/2024 Makoto. medianReref was not used. Fixed.
% 11/27/2024 Makoto. Used again.
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


%% Re-do 0292, 2384.

subjNames(6);  % 0292
subjNames(42); % 2384 

for rawIdx = 1:length(allRaws)

    try

        % Import EGI data.
        loadName = allRaws(rawIdx).name;
        dataName = loadName(1:4);

        % if strcmp(dataName, '0292') | strcmp(dataName, '2384')
        if strcmp(dataName, '0292') | strcmp(dataName, '2384')
            EEG = pop_readegi(['/srv/Makoto/ASSR/p0000_raw' filesep loadName], [],[],'auto');
        else
            continue
        end

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
        % EEG = pop_resample(EEG, 250);

        % Apply 0.5-Hz HPF (TBW 1.0)
        EEG = pop_firws(EEG, 'fcutoff', 0.5, 'ftype', 'highpass', 'wtype', 'blackman', 'forder', 5500, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

        % % Apply 52.5-Hz LPF (TBW 5.0)
        EEG = pop_firws(EEG, 'fcutoff', 52.5, 'ftype', 'lowpass', 'wtype', 'blackman', 'forder', 276, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);

        % Remove line noise up to Nyquist.
        %EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist', 1:EEG.nbchan,'computepower',1,'linefreqs',60:60:480,'newversion',1,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

        % Discard 3 second of data from both edges.
        EEG = pop_select( EEG, 'time', [3 EEG.xmax-3]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Detect outlier channels (experimental, simple) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the initial flatness.
        chStd = std(EEG.data,0,2);
        flatChThreshold  = median(chStd)-mad(chStd,1)*4; % Robust 4SD.
        if flatChThreshold<=0
            flatChThreshold = 0.1; % heuristic choice.
        end
        flatChannelIdx = find(chStd<flatChThreshold);

        % Re-reference to median.
        medianPotential = median(EEG.data,1);
        medianRerefData = bsxfun(@minus, EEG.data, medianPotential);

        % Segment the data into 1000 epochs.
        permillWindowLength = floor(EEG.pnts/1000);
        trimmingLength      = mod(EEG.pnts, permillWindowLength);
        %trimmedData         = EEG.data(:,1:end-trimmingLength);
        trimmedData         = medianRerefData(:,1:end-trimmingLength);
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

        % Calculate robust 3SD channel rejection threshold for large amplitude.
        channelRejectionThreshold = median(median(std3D,2))+20*mad(median(std3D,2),1)*1.4826; % 25 or 30; For subject 1, 25 is about the largest eye channel value.         
        largeAmplitudeChannelIdx  = find(median(std3D,2)>channelRejectionThreshold);

        % % Find flat channels.
        % flatChannelIdx = find(mean(std3D,2)<median(mean(std3D,2))/10);

        % Combine channel indices to exclude.
        %excludingChannelIdx = union(flatChannelIdx, largeAmplitudeChannelIdx); % 11/27/2024 Makoto. Suggested by Hyeonseok.
        excludingChannelIdx = flatChannelIdx; % Empirically confirmed.

        % Exclude outlier channels.
        EEG.etc.rejectedChannelIdx = excludingChannelIdx;
        EEGorig = EEG;
        EEG = pop_select(EEG, 'nochannel', excludingChannelIdx);


        %%%%%%%%%%%%
        %%% ASR. %%%
        %%%%%%%%%%%%
        EEG_beforeASR = EEG;
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off', ...
            'LineNoiseCriterion','off','Highpass','off', ...
            'BurstCriterion',25, 'WindowCriterion', 'off', ...
            'BurstRejection','off','Distance','Euclidian', ...
            'BurstCriterionRefMaxBadChns', 0, ...
            'BurstCriterionRefTolerances', [-inf 8], ...
            'WindowCriterionTolerances',[-Inf 8], 'MaxMem', 65536);
        if isfield(EEG.etc, 'clean_sample_mask')
            clean_sampleIdx = find(EEG.etc.clean_sample_mask);
        else
            clean_sampleIdx = 1:EEG.pnts;
            EEG.etc.clean_sample_mask = ones(1,EEG.pnts);
        end
        survivedDataIdx           = find(EEG.etc.clean_sample_mask);
        rejectedDataIdx           = find(~EEG.etc.clean_sample_mask);
        asrBeforeAfterDiff        = sum(EEG_beforeASR.data(:,survivedDataIdx)-EEG.data,1);
        unchangedDataIdx          = find(asrBeforeAfterDiff==0);
        changedDataIdx            = find(asrBeforeAfterDiff~=0);
        windowRejRate             = 1-EEG.pnts/EEG_beforeASR.pnts;
        windowInterpolationRate   = 1-(EEG.pnts-length(unchangedDataIdx))/EEG.pnts;
        asrPowerReductionDb       = 10*log10(var(EEG.data(:,changedDataIdx ),0,2)./var(EEG_beforeASR.data(:,changedDataIdx),0,2));
        windowRejPowReducDb       = 10*log10(var(EEG.data(:,changedDataIdx ),0,2)./var(EEG_beforeASR.data(:,rejectedDataIdx),0,2));
        EEG.etc.ASR.windowRejectionRate     = windowRejRate;
        EEG.etc.ASR.windowInterpolationRate = windowInterpolationRate;
        EEG.etc.ASR.varianceReductionInDbByWinRej = windowRejPowReducDb;
        EEG.etc.ASR.varianceReductionInDbByAsr    = asrPowerReductionDb;

        % Interpolate the rejected channels.
        EEG = pop_interp(EEG, EEGorig.chanlocs, 'spherical');

        % Perform average reference.
        averefData = bsxfun(@minus, EEG.data, sum(EEG.data,1)/(EEG.nbchan+1));
        averefData = bsxfun(@minus, averefData, mean(averefData,2));
        EEG.data   = averefData;
        EEG.ref    = 'average';

        % % AMICA.
        % EEG_downsampled = pop_resample(EEG, 100);
        % numICs = sum(eig(cov(double(EEG_downsampled.data')))>1E-6);
        % runamica15(double(EEG_downsampled.data), 'outdir', sprintf('/srv/Makoto/ASSR/p0100_upToDipfit/%s', dataName), ...
        %             'max_iter', 2000, 'max_threads', 16, ...
        %             'pcakeep', numICs, 'writestep', 2000, ...
        %             'do_reject', 1, 'numrej', 10, 'rejsig', 3, 'rejstart', 1, 'rejint', 1);
        % EEG.etc.amica  = loadmodout15(sprintf('/srv/Makoto/ASSR/p0100_upToDipfit/%s', dataName));
        % EEG.icaweights = EEG.etc.amica.W;
        % EEG.icasphere  = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs,:);
        % EEG.icawinv    = EEG.etc.amica.A;
        % EEG.icachansind = 1:EEG.nbchan;
        % 
        % % Run ICLabel.
        % EEG = pop_iclabel(EEG, 'default');
        % 
        % % Run dipfit.
        % EEG = pop_dipfit_settings( EEG, 'hdmfile', '/srv/TOOLKITS/eeglab-2022.0/plugins/dipfit4.3/standard_BEM/standard_vol.mat', ...
        %                                 'mrifile', '/srv/TOOLKITS/eeglab-2022.0/plugins/dipfit4.3/standard_BEM/standard_mri.mat', ...
        %                                 'chanfile','/srv/TOOLKITS/eeglab-2022.0/plugins/dipfit4.3/standard_BEM/elec/standard_1005.elc', ...
        %                                 'coord_transform',[0.05476 -17.3653 -8.1318 0.075502 0.0031836 -1.5696 11.7138 12.7933 12.213] , ...
        %                                 'coordformat','MNI', 'chansel', 1:EEG.nbchan);
        % EEG = pop_multifit(EEG, 1:numICs ,'threshold',100,'plotopt',{'normlen','on'});
        % 
        % % Run bilateral dipole fitting.
        % EEG = fitTwoDipoles(EEG, 'LRR', 35);

        % Save.
        EEG = pop_saveset( EEG, 'filename', dataName, 'filepath', ['/srv/Makoto/ASSR/p0104_upToDipfit_ASRstats/']);

    catch

        dummyData = 0;
        save(sprintf('/srv/Makoto/ASSR/p0104_upToDipfit_ASRstats/ERROR_%s', dataName), 'dummyData')
        continue
    end
end


%% Summarize the results.
allSets = dir('/srv/Makoto/ASSR/p0104_upToDipfit_ASRstats/*.set');
winInterpRateVec = zeros(length(allSets),1);
varReducAsr    = nan(128, length(allSets));
varReducWinRej = nan(128, length(allSets));
for setIdx = 1:length(allSets)

    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    loadPath = allSets(setIdx).folder;
    EEG = pop_loadset('filename', loadName, 'filepath', loadPath, 'loadmode', 'info');

    winInterpRateVec(setIdx) = EEG.etc.ASR.windowInterpolationRate;

    if isempty(EEG.etc.rejectedChannelIdx)
        varReducAsr(:,setIdx) = EEG.etc.ASR.varianceReductionInDbByAsr;
    else
        presentChIdx = setdiff(1:EEG.nbchan, EEG.etc.rejectedChannelIdx);
        varReducAsr(   presentChIdx,setIdx) = EEG.etc.ASR.varianceReductionInDbByAsr;
        varReducWinRej(presentChIdx,setIdx) = EEG.etc.ASR.varianceReductionInDbByWinRej;
    end
end

% Replace Inf with a group-median value.
[I,J] = ind2sub(size(varReducAsr), find(varReducAsr==inf))
varReducAsr(I,J) = median(varReducAsr(I,:));

save('/srv/Makoto/ASSR/p0104_upToDipfit_ASRstats/asrStats', 'winInterpRateVec', 'varReducAsr', 'varReducWinRej')

figure
bar(nanmedian(varReducAsr,2))

figure
topoplot(median(varReducAsr,2), EEG.chanlocs)
topoplot(median(varReducAsr,2), EEG.chanlocs, 'maplimits', [-12 0])

sort(varReducAsr,1)