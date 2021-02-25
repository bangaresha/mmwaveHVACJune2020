function [errStats,pathLoss,ebnoEst,bernumber,berratio,Capacity] = BPSKwithRFImpairmentsSim(...
    dist,freq,atmosCond,eff,diaTx,diaRx,hpaBackoff,rxTemp, ...
    showConstellation,showSpectrum,showHPA,resetBER,channel,channelTime)

% dist                  - link distance (km)
% freq                  - carrier frequency (MHz)
% eff                   - antenna efficiency (0<eff<1)
% diaTx                 - transmit antenna diameter (m)
% diaRx                 - receive antenna diameter (m)
% hpaBackoff            - HPA backoff from saturation (dB)
% rxTemp                - receiver noise temperature (K)
% fD                    - Doppler shift (Hz)
% ampImb                - I/Q amplitude imbalance (dB)
% phImb                 - I/Q phase imbalance (deg)
% dcOffset              - complex DC offset (V)
% phaseCorrect          - {0, -12, -16} (deg)
% correctDCOffset       - true to enable, false to disable
% correctIQImbalance    - true to enable, false to disable
% correctDoppler        - true to enable, false to disable
% showConstellation     - true to enable, false to disable
% showSpectrum          - true to enable, false to disable
% showHPA               - true to enable, false to disable
% resetBER              - true if error counter should be reset


% Establish persistent System objects
persistent TX_FILT RX_FILT HPA AGC PF_OFFSET ...
    THERMALNOISE PHASENOISE ERRCNT IQCOMP DCBLOCK CARRIERSYNC

% Exclude scopes from code generation
coder.extrinsic('rfPlotConstellation')
coder.extrinsic('rfPlotSpectrum')
coder.extrinsic('rfPlotHPA')
coder.extrinsic('rng')
modOrder = 2;  % BPSK modulation
freq = (str2double(freq))*1E9;
fs = 2*freq;
% Create System objects the first time the function is run
if isempty(TX_FILT)
    
    % Raised root cosine transmit and receive filter objects having 8
    % samples per symbol (sps)
    sps = 2;
    TX_FILT = comm.RaisedCosineTransmitFilter(...
        'RolloffFactor',0.3, ...
        'FilterSpanInSymbols',1, ...
        'OutputSamplesPerSymbol',sps, ...
        'Gain',sqrt(sps));
    RX_FILT = comm.RaisedCosineReceiveFilter(...
        'RolloffFactor',0.3, ...
        'FilterSpanInSymbols',1, ...
        'InputSamplesPerSymbol',sps,'DecimationFactor',sps, ...
        'Gain',1/sqrt(sps));
    
    % Create a memoryless nonlinearity object to model a high power
    % amplifier (HPA). Call the supporting function rfHPAGains
    % to determine the input and output scaling for the HPA.
    [hpaGainIn,hpaGainOut] = rfHPAGains(hpaBackoff);
    HPA = comm.MemorylessNonlinearity(...
        'Method','Saleh model', ...
        'InputScaling',hpaGainIn, ...
        'OutputScaling',hpaGainOut);
    
    % Create an AGC object.
    AGC = comm.AGC('AveragingLength',16, 'MaxPowerGain',10);
    
    % Create phase frequency offset objects to model the frequency offset
    % that arises from a Doppler shift and to apply a frequency and/or
    % phase correction to the received signal.
    PF_OFFSET = comm.PhaseFrequencyOffset(...
        'FrequencyOffsetSource','Input port','SampleRate',fs);
    
    % Create a thermal noise object to add receiver noise.
    THERMALNOISE = comm.ThermalNoise('SampleRate',fs, ...
       'NoiseTemperature',rxTemp);
    
    % Create an error rate object to calculate the error statistics.
    % Calculate the receive delay based on the filter span and modulation
    % order. Set a computation delay of 5000 to allow the AGC sufficient
    % time to converge.
    delay = TX_FILT.FilterSpanInSymbols*log2(modOrder);
    ERRCNT = comm.ErrorRate('ReceiveDelay',delay,'ComputationDelay',5000);

    % Create a carrier synchronizer object to compensate for Doppler
    CARRIERSYNC = comm.CarrierSynchronizer( ...
        'DampingFactor', 0.707, ...
        'NormalizedLoopBandwidth', 0.002, ...
        'SamplesPerSymbol', 2, ...
        'Modulation', 'BPSK', ...
        'ModulationPhaseOffset', 'Auto');
      
    % Set the random number generator to default to be able to repeat
    % simulation results.
    rng default
    % Set flag to indicate first time the function is called
    firstCall = true;
else
    % Set the tunable System object properties to new values (if
    % applicable). The properties can be changed during the simulation by
    % adjusting the GUI.
    [hpaGainIn,hpaGainOut] = rfHPAGains(hpaBackoff);
    HPA.InputScaling = hpaGainIn;
    HPA.OutputScaling = hpaGainOut;
    
    THERMALNOISE.NoiseTemperature = rxTemp;
    
    firstCall = false;
    
    % Reset the error rate counter if flag is true. This is done if an
    % input value is changed, since BER estimates are meaningless without
    % consistent input conditions. 
    if resetBER
        reset(ERRCNT)
    end
end

% Calculate the transmit and receive antenna gains as well as the path
% loss.
%lightSpeed = physconst('light')
lightSpeed = 2.98E8;
%freq = 60E9;

waveLength = lightSpeed/freq;
txAntGain = sqrt(eff)*pi*str2double(diaTx)/waveLength;
rxAntGain = sqrt(eff)*pi*str2double(diaRx)/waveLength;
%screen 100 -> 0.1km
dist = dist*1e-2;
freeSpacePL = fspl(dist, waveLength);
T = 15; % Temperature in degree C
pathLoss = freeSpacePL;  
% switch atmosCond % Get path loss in dB
%     case 'fg'   % Fog
%         den = .05; % Liquid water density in g/m^3
%         % Approximate maximum 18km for fog/cloud
%         pathLoss = freeSpacePL + ...
%             fogpl( min(dist, 18)*1000, freq*1e6, T, den); 
%     case 'fs'   % Free space
%         pathLoss = freeSpacePL;  
%     case 'gs'   % Gas
%         P = 101.325e3; % Dry air pressure in Pa
%         den = 7.5;     % Water vapor density in g/m^3
%         % Approximate maximum 100km for atmospheric gases
%         pathLoss = freeSpacePL + ...
%             gaspl( min(dist, 100)*1000, freq*1e6, T, P, den);
%     otherwise   % Rain
%         RR = 3; % Rain rate in mm/h
%         % Approximate maximum 2km for rain
%         pathLoss = freeSpacePL + ...
%             rainpl(min(dist, 2)*1000, freq*1e6, RR);
% end

% Main Processing Loop
% Random data is modulated, filtered, amplified, and transmitted through a
% satellite communications channel. Impairments such as a Doppler shift,
% phase noise, DC offset, and an I/Q phase imbalance are introduced.
% Compensation techniques are applied and the bit error rate (BER) is
% calculated.

% Transmit operations
dataIn = randi([0 1],16000,1); % Generate binary data
modData = pskmod(dataIn, modOrder);
txFiltOut = TX_FILT(modData); % Filter signal
hpaOut = HPA(txFiltOut); % Amplify signal with HPA
txSig = txAntGain*hpaOut; % Apply transmit antenna gain

% rxSig = conv(txSig,channelTime);
% rxSig = rxSig(1:32000);

% Channel operations
rxSig = txSig/10^(pathLoss/20); % Apply free-space loss
rxSig = rxSig(1:32000);


rxSigGain = rxAntGain*rxSig; % Apply receive antenna gain
rxSigAWGN = THERMALNOISE(rxSigGain); % Add AWGN % rxSigGain;
% Estimate Eb/No
snrdB = 10*log10(var(rxSigGain)./var(rxSigAWGN-rxSigGain));
ebnoEst = snrdB + ...
    10*log10(TX_FILT.OutputSamplesPerSymbol) - ...
    10*log10(log2(modOrder));

AGCOut = AGC(rxSigAWGN); % Apply AGC
rxFiltOut = RX_FILT(AGCOut); % Receive filter
dataOut = pskdemod(rxFiltOut, modOrder);  % Demodulate
Pse = 1E-5;
Kn = 2;
deltaF = 0.3125E6;
bandwidth = 80E6;
noisePower = 80e-6;       % AWGN noise power
% deltaF = 40E6;
% bandwidth = 2.16E9;
% noisePower = 2.16e-3;       % AWGN noise power

snrTruedB = 30;
snr = (10.^(0.1*snrTruedB));
txPower = noisePower*snr/(dist^2);  
S = svd(channel);
channel = S.^2/noisePower;

% channel = abs(channel).^2;
% cnr = channel./noisePower;

K0 = (3*bandwidth*txPower)/((erfinv(Pse/Kn))^2);
% K0 = (3*bandwidth*(10.^(0.1*snrTruedB)))/((erfinv(Pse/Kn))^2);
allocatedPower = waterFilling(channel,deltaF,K0);
for k = 1:length(channel)
    CapacityVec(k) = deltaF*log2((allocatedPower*K0*channel(k))/log(2));
end
Capacity = sum(CapacityVec);
% Display visuals
rfPlotSpectrum(txFiltOut,rxSigAWGN,showSpectrum,firstCall)
%rfPlotHPA(txFiltOut,hpaOut,showHPA,firstCall)
% release(coder)
% rfPlotConstellation(rxFiltOut,showConstellation,firstCall)

% Calculate the error statistics.
[bernumber,berratio] = biterr(dataIn,dataOut(1:size(dataIn,1)));
errStats = ERRCNT(dataIn,dataOut(1:size(dataIn,1)));
errStats(1) = errStats(1)/10000;
errStats(2) = errStats(2)/10000;


% EOF