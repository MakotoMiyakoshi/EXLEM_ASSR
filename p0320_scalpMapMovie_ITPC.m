% 07/09/2024 Makoto. Created.
clear
close all

% Load dummy EEG.
EEG = pop_loadset('filename','0009.set','filepath','/srv/Makoto/ASSR/p0100_upToDipfit/', 'loadmode', 'info');

% Load the demographic data.
demographicData = readtable('/srv/Makoto/ASSR/code/Subset_SSCT_for_Ernie.xlsx');
groupIdx = demographicData.Dx;
tdcIdx = find(strcmp(groupIdx, 'TDC'));
fxsIdx = find(strcmp(groupIdx, 'FXS'));

% Generate frequency bins.
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
timeBins = -1000:3999;


%% Median Gamma-band power topo.
allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat])
    stackData(:,:,:,matIdx) = elecItcMedian;
end
itcFXS        = squeeze(mean(stackData(:,freqIdx7,:,fxsIdx),4));
itcTDC        = squeeze(mean(stackData(:,freqIdx7,:,tdcIdx),4));

%% Prepare video output parameters.
outputVideo = VideoWriter('/srv/Makoto/ASSR/p0320_scalpMapMovie_ITPC/gammaItcTopo', 'Motion JPEG AVI');
outputVideo.FrameRate = 20;
outputVideo.Quality   = 95;
open(outputVideo)

figHandle = figure('position', [200 200 800 800], 'visible', 'on');
plottingTime = -200:3500;
for timeIdx = 1:10:length(plottingTime)

    disp(sprintf('%d/%d...', timeIdx, length(plottingTime)))

    ax1 = subplot(2,2,1);
    topoplot(itcFXS(:,800+timeIdx), EEG.chanlocs, 'maplimits', [-0.4 0.4]); title('FXS')

    ax2 = subplot(2,2,2);
    topoplot(itcTDC(:,800+timeIdx), EEG.chanlocs, 'maplimits', [-0.4 0.4]); title('TDC')
    originalPosition = get(gca,'position');
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
    set(gca,'position', originalPosition)

    if timeIdx == 1
        ax3 = subplot(2,2,3);
        plot(timeBins, itcFXS', 'color', [0.85 0.85 0.85])
        xlim([plottingTime(1) plottingTime(end)])
        ylim([0.05 0.45])
        line([0    0],    ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', '-')
        line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')

        ax4 = subplot(2,2,4);
        plot(timeBins, itcTDC', 'color', [0.85 0.85 0.85])
        xlim([plottingTime(1) plottingTime(end)])
        ylim([0.05 0.45])
        line([0    0],    ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', '-')
        line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    end

    axes(ax3)
    movingLineHandle1 = line([plottingTime(timeIdx) plottingTime(timeIdx)], ylim, 'color', [1 0 0], 'linewidth', 1, 'linestyle', '-');
    axes(ax4)
    movingLineHandle2 = line([plottingTime(timeIdx) plottingTime(timeIdx)], ylim, 'color', [1 0 0], 'linewidth', 1, 'linestyle', '-');

    sgtitle(sprintf('40Hz ITPC at %d ms', plottingTime(timeIdx)))

    % Write the current figure to a frame of the video.
    currentMovieFrame = getframe(figHandle);
    writeVideo(outputVideo, currentMovieFrame);

    % Clear the current axes.
    cla(ax1)
    cla(ax2)
    delete(movingLineHandle1)
    delete(movingLineHandle2)
end
close(outputVideo)



%% Median Theta-band itc topo.
thetaIdx = find(wtFreqBins>3 & wtFreqBins<8);

allMats = dir('/srv/Makoto/ASSR/p0200_epoch/*_elecItcMedian.mat');
stackData = zeros(128, 100, 5000, length(allMats));
for matIdx = 1:length(allMats)
    disp(sprintf('%d/%d', matIdx, length(allMats)))
    currentMat = allMats(matIdx).name;
    load(['/srv/Makoto/ASSR/p0200_epoch/' currentMat]);
    stackData(:,:,:,matIdx) = elecItcMedian;
end
itcFXS = squeeze(mean(mean(stackData(:,thetaIdx,:,fxsIdx),4),2));
itcTDC = squeeze(mean(mean(stackData(:,thetaIdx,:,tdcIdx),4),2));

%% Prepare video output parameters.
outputVideo = VideoWriter('/srv/Makoto/ASSR/p0320_scalpMapMovie_ITPC/thetaItcTopo', 'Motion JPEG AVI');
outputVideo.FrameRate = 20;
outputVideo.Quality   = 95;
open(outputVideo)

figHandle = figure('position', [200 200 800 800], 'visible', 'off');
plottingTime = -200:3500;
for timeIdx = 1:10:length(plottingTime)

    disp(sprintf('%d/%d...', timeIdx, length(plottingTime)))

    ax1 = subplot(2,2,1);
    topoplot(itcFXS(:,800+timeIdx), EEG.chanlocs, 'maplimits', [-0.3 0.3]); title('FXS')

    ax2 = subplot(2,2,2);
    topoplot(itcTDC(:,800+timeIdx), EEG.chanlocs, 'maplimits', [-0.3 0.3]); title('TDC')
    originalPosition = get(gca,'position');
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', 'dB')
    set(gca,'position', originalPosition)

    if timeIdx == 1
        ax3 = subplot(2,2,3);
        plot(timeBins, itcFXS', 'color', [0.85 0.85 0.85])
        xlim([plottingTime(1) plottingTime(end)])
        ylim([0.05 0.35])
        line([0    0],    ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', '-')
        line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')

        ax4 = subplot(2,2,4);
        plot(timeBins, itcTDC', 'color', [0.85 0.85 0.85])
        xlim([plottingTime(1) plottingTime(end)])
        ylim([0.05 0.35])
        line([0    0],    ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', '-')
        line([3000 3000], ylim, 'color', [0 0 0], 'linewidth', 2, 'linestyle', ':')
    end

    axes(ax3)
    movingLineHandle1 = line([plottingTime(timeIdx) plottingTime(timeIdx)], ylim, 'color', [1 0 0], 'linewidth', 1, 'linestyle', '-');
    axes(ax4)
    movingLineHandle2 = line([plottingTime(timeIdx) plottingTime(timeIdx)], ylim, 'color', [1 0 0], 'linewidth', 1, 'linestyle', '-');

    sgtitle(sprintf('3-8Hz ITC at %d ms', plottingTime(timeIdx)))

    % Write the current figure to a frame of the video.
    currentMovieFrame = getframe(figHandle);
    writeVideo(outputVideo, currentMovieFrame);

    % Clear the current axes.
    cla(ax1)
    cla(ax2)
    delete(movingLineHandle1)
    delete(movingLineHandle2)
end
close(outputVideo)