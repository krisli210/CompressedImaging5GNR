close all

PdStructure = load("data\OMP_GridSearch\PdGridCut.mat");
PfStructure = load("data\OMP_GridSearch\PfGridCut.mat");

PdGrid = PdStructure.PdGrid;
PfGrid = PfStructure.PfGrid;

KRange = PdStructure.KRangeCut;
PdGrid = PdGrid(1:length(KRange), :, :, :);
PfGrid = PfGrid(1:length(KRange), :, :, :);

packetRange = PdStructure.packetRange;
userRange = PdStructure.userRange;
NTargetRange = PdStructure.NTargetRange;

KStr = string([repmat('K = ', length(KRange), 1), num2str(KRange.')]);
packetString = string([repmat('T = ', length(packetRange), 1), num2str(packetRange.')]);
userString = string([repmat('U = ', length(userRange), 1), num2str(userRange.')]);
NTargetString = string([repmat('N = ', length(NTargetRange), 1), num2str(NTargetRange.')]);

figure; hold on; grid on;

KFixed = 272; KIndex = find(KRange == KFixed);
packetsFixed = 9; packetIndex = find(packetRange == packetsFixed);
usersFixed= 3; userIndex = find(userRange == usersFixed);
NTargetsFixed = 9; NTargetsIndex = find(NTargetRange == NTargetsFixed);

yyaxis left;
plot(NTargetRange, squeeze(PdGrid(KIndex, packetIndex, userIndex, :)), '-x');
ylim([0 1]); ylabel('P_D')

yyaxis right;
plot(NTargetRange, squeeze(PfGrid(KIndex, packetIndex, userIndex, :)), '-x');
ylim([0 1]); ylabel('P_F');

