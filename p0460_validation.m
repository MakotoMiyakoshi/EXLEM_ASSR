% 07/12/2024 Makoto. Created.

originalData      = load('/srv/Makoto/ASSR/p0210_epochForERP/meanMedianERPs.mat');
brainClassData    = load('/srv/Makoto/ASSR/p0440_epochBrainForERP/meanMedianERPs_brainClass.mat');
nonBrainClassData = load('/srv/Makoto/ASSR/p0450_epochNonBrainForERP/meanMedianERPs_nonBrainClass.mat');

brainPlusNonBrain = brainClassData.meanERP + nonBrainClassData.meanERP;

residual = originalData.meanERP - brainPlusNonBrain;


%%
figure('position', [100         600        2350         435])
subplot(1,4,1)
plot(-1000:3999, [mean(originalData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([-1000 4000])
ylim([-3.5 2.5])
title('Original ERP at Cz (n=44)')
ylabel('Amplitude (\muV)')
xlabel('Latency (ms)')
title('Original ERP at Cz (n=44)')

subplotPosition = get(gca,'position');
insetRelativePosition = [0.6 0.05 0.35 0.3]; % SouthEast. 
deltax = subplotPosition(3)*insetRelativePosition(1);
deltay = subplotPosition(4)*insetRelativePosition(2);
insetPosition = [subplotPosition(1:2) 0 0] + [deltax deltay subplotPosition(3:4).*insetRelativePosition(3:4)];
axes('Position',insetPosition)
plot(-1000:3999, [mean(originalData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([400 600])
ylim([-1.5 1.5])


subplot(1,4,2)
plot(-1000:3999, [mean(brainClassData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([-1000 4000])
ylim([-3.5 2.5])
title('Brain-class IC ERPs at Cz')
ylabel('Amplitude (\muV)')
xlabel('Latency (ms)')

subplotPosition = get(gca,'position');
insetRelativePosition = [0.6 0.05 0.35 0.3]; % SouthEast. 
deltax = subplotPosition(3)*insetRelativePosition(1);
deltay = subplotPosition(4)*insetRelativePosition(2);
insetPosition = [subplotPosition(1:2) 0 0] + [deltax deltay subplotPosition(3:4).*insetRelativePosition(3:4)];
axes('Position',insetPosition)
plot(-1000:3999, [mean(brainClassData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([400 600])
ylim([-1.5 1.5])


subplot(1,4,3)
plot(-1000:3999, [mean(nonBrainClassData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([-1000 4000])
ylim([-3.5 2.5])
title('Non-Brain-class IC ERPs at Cz')
ylabel('Amplitude (\muV)')
xlabel('Latency (ms)')

subplotPosition = get(gca,'position');
insetRelativePosition = [0.6 0.05 0.35 0.3]; % SouthEast. 
deltax = subplotPosition(3)*insetRelativePosition(1);
deltay = subplotPosition(4)*insetRelativePosition(2);
insetPosition = [subplotPosition(1:2) 0 0] + [deltax deltay subplotPosition(3:4).*insetRelativePosition(3:4)];
axes('Position',insetPosition)
plot(-1000:3999, [mean(nonBrainClassData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([400 600])
ylim([-1.5 1.5])


subplot(1,4,4)
plot(-1000:3999, [mean(originalData.meanERP(:,:,:),3)-mean(brainClassData.meanERP(:,:,:),3)-mean(nonBrainClassData.meanERP(:,:,:),3)]', 'color', [0.7 0.7 0.7])
xlim([-1000 4000])
ylim([-3.5 2.5])
title('Original - (Brain + Non-Brain)')

print('/srv/Makoto/ASSR/p0460_validation/validation_ERP', '-r200', '-djpeg95')