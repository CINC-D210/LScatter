function rxWaveformCorrected = waveformCorrect(rxWaveform,enb)
%waveformCorrect 修正频偏和时延
    %% Derived parameters
    samplesPerFrame = 10e-3*enb.SamplingRate; 
    %% Channel estimation configuration structure
    cec.PilotAverage = 'UserDefined';  % Type of pilot symbol averaging
    cec.FreqWindow = 9;                % Frequency window size in REs
    cec.TimeWindow = 9;                % Time window size in REs
    cec.InterpType = 'Cubic';          % 2D interpolation type
    cec.InterpWindow = 'Centered';     % Interpolation window type
    cec.InterpWinSize = 3;             % Interpolation window size
    %% Receiver processing    
    % Perform frequency offset correction for known cell ID
    frequencyOffset = lteFrequencyOffset(enb,rxWaveform);%%%
    rxWaveform = lteFrequencyCorrect(enb,rxWaveform,frequencyOffset);
%     fprintf('\nCorrected a frequency offset of %i Hz.\n',frequencyOffset)
    
    % Perform the blind cell search to obtain cell identity and timing offset 
    % Use 'PostFFT' SSS detection method to improve speed
    cellSearch.SSSDetection = 'PostFFT'; cellSearch.MaxCellCount = 1;
    NCellID = lteCellSearch(enb,rxWaveform,cellSearch);
%     fprintf('Detected a cell identity of %i.\n', NCellID);
    %     enb.NCellID = NCellID; % From lteCellSearch
    
    [frameOffset,corr] = lteDLFrameOffset(enb,rxWaveform);%%%
    
%     subplot(2,2,1);
%     plot(corr);
%     hold on;
%     % subplot(1,3,1),plot(1:length(rxWaveform),[zeros(frameOffset,1);0.18;zeros(length(rxWaveform)-frameOffset-1,1)]);
%     hold off;
%     axis square;
%     title('Packet Detection');
%     drawnow;
    
    % Sync the captured samples to the start of an LTE frame, and trim off any samples that are part of an incomplete frame.
    rxWaveform2 = rxWaveform(frameOffset+1:frameOffset+19200); % [307200x1] -> [153600x1]
    tailSamples = mod(length(rxWaveform2),samplesPerFrame);
    rxWaveformCorrected = rxWaveform2(1:end-tailSamples,:);
    enb.NSubframe = 0;
%     fprintf('Corrected a timing offset of %i samples.\n',frameOffset)
end

