% 06/26/2024 Makoto. Created.

clear global
clear all
close all

allSets = dir('/srv/Makoto/ASSR/p0400_1peakFOOOF/*.set');
for setIdx = 1:length(allSets)

    % Load the current data.
    loadName = allSets(setIdx).name;
    dataName = loadName(1:4);
    EEG = pop_loadset('filename', loadName, 'filepath', '/srv/Makoto/ASSR/p0400_1peakFOOOF');
    EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data(:,:);
    EEG.icaact = reshape(EEG.icaact, [size(EEG.icaweights,1), EEG.pnts, size(EEG.data,3)]);

    % Obtain Brain ICs.
    [classProb,classIdx] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
    brainClassIdx = find(classIdx==1 & classProb>0.5);

    % Obtain FOOOF-certified ICs.
    expSlopeCoeffVec = zeros(size(EEG.icaweights,1),1);
    for icIdx = 1:size(EEG.icaweights,1)
        expSlopeCoeffVec(icIdx) = EEG.etc.FOOOF_concatenated{icIdx}.aperiodic_params(2);
    end
    oneOverFIdx = find(expSlopeCoeffVec>0.5);
        %{
        figure
        scatter(classIdx, expSlopeCoeffVec)
        %}

    % Find good brain-class ICs.
    brainIcIdx    = intersect(brainClassIdx, oneOverFIdx);
    nonBrainIcIdx = setdiff(1:size(EEG.icaweights,1), brainIcIdx);

    % Extract brain- and non-brain-class ICs.
    EEGorig = EEG;
    EEG = pop_subcomp(EEG, brainIcIdx, 0, 1);
    pop_saveset(EEG, 'filename', dataName, 'filepath', '/srv/Makoto/ASSR/p0410_selectBrainClassICs/brainClass')

    EEG = pop_subcomp(EEGorig, nonBrainIcIdx, 0, 1);
    pop_saveset(EEG, 'filename', dataName, 'filepath', '/srv/Makoto/ASSR/p0410_selectBrainClassICs/nonBrainClass')
end

