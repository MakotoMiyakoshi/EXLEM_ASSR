onsetIdx = find(strcmp({EEG.event.type}, 'DI66'));

latencyVec = [EEG.event(onsetIdx).latency];

diffLatency = diff(latencyVec);

diffInMs = diffLatency*(1000/EEG.srate);

figure; hist(diffInMs)