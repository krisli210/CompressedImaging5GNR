function [waveform, info] = genMU5GNRWaveform
    % taken from nrwaveformgenerator doc

carriers = {
    nrSCSCarrierConfig('SubcarrierSpacing',60,'NStartGrid',10,'NSizeGrid',100), ...
    nrSCSCarrierConfig('SubcarrierSpacing',120,'NStartGrid',0,'NSizeGrid',70)};

bwp = {
    nrWavegenBWPConfig('BandwidthPartID',1,'SubcarrierSpacing',60,'NStartBWP',30,'NSizeBWP',80), ...
    nrWavegenBWPConfig('BandwidthPartID',2,'SubcarrierSpacing',120,'NStartBWP',0,'NSizeBWP',60)};

numTransmitted = 0;
ssb = nrWavegenSSBurstConfig('BlockPattern','Case D', 'TransmittedBlocks', [ones(1, numTransmitted), zeros(1, 64-numTransmitted)], 'SubcarrierSpacingCommon', 120);

pdcch = {
    nrWavegenPDCCHConfig('SearchSpaceID',1,'BandwidthPartID',1,'RNTI',1,'DMRSScramblingID',1), ...
    nrWavegenPDCCHConfig('SearchSpaceID',2,'BandwidthPartID',2,'RNTI',2,'DMRSScramblingID',2, ...
        'AggregationLevel',4)};

coreset = {
    nrCORESETConfig('CORESETID',1,'FrequencyResources',[1 1 1 1 1 0 0 0 0 0 1],'Duration',3), ...
    nrCORESETConfig('CORESETID',2,'FrequencyResources',[0 0 0 0 0 0 0 0 1 1])};

ss = {
    nrSearchSpaceConfig('SearchSpaceID',1,'CORESETID',1,'StartSymbolWithinSlot',4), ...
    nrSearchSpaceConfig('SearchSpaceID',2,'CORESETID',2,'NumCandidates',[8 8 4 0 0])};

pdsch = {
    nrWavegenPDSCHConfig('BandwidthPartID',1,'Modulation','16QAM','RNTI',1,'NID',1), ...
    nrWavegenPDSCHConfig('BandwidthPartID',2,'Modulation','QPSK','RNTI',2,'NID',2, ...
            'PRBSet', 50:59)};

csirs = {
     nrWavegenCSIRSConfig('BandwidthPartID',1,'RowNumber',2,'RBOffset',10), ... 
     nrWavegenCSIRSConfig('BandwidthPartID',2,'Density','three','RowNumber',4)};

cfgDL = nrDLCarrierConfig( ...
    'FrequencyRange','FR2', ...
    'ChannelBandwidth',200, ...
    'NumSubframes',20, ...
    'SCSCarriers',carriers, ...
    'BandwidthParts',bwp, ...
    'CORESET',coreset, ...
    'SSBurst', ssb, ...
    'SearchSpaces',ss, ...
    'PDCCH',pdcch, ...
    'PDSCH',pdsch, ...
    'CSIRS',csirs);

[waveform, info] = nrWaveformGenerator(cfgDL);
end