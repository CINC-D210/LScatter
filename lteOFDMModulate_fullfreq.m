%lteOFDMModulate OFDM modulation
%   [WAVEFORM,INFO] = lteOFDMModulate(...) performs DC subcarrier
%   insertion, IFFT calculation, cyclic prefix insertions and optional
%   raised-cosine windowing and overlapping of adjacent OFDM symbols.
%
%   [WAVEFORM,INFO] = lteOFDMModulate(ENB,GRID) performs DC subcarrier
%   insertion, IFFT calculation, cyclic prefix insertion and optional
%   raised-cosine windowing and overlapping of adjacent OFDM symbols of the
%   complex symbols in the resource array GRID. GRID is a 3-dimensional
%   array of the Resource Elements for a number of subframes across all
%   configured antenna ports, as described in the <a href="matlab: web([docroot '/lte/gs/data-structures.html'])">Data Structures</a> 
%   documentation, or multiple such matrices concatenated to give multiple
%   subframes (concatenation across the columns or 2nd dimension). The
%   antenna planes in GRID are each OFDM modulated to give the columns of
%   WAVEFORM.
%   ENB must be a structure including the following fields:
%   CyclicPrefix - Optional. Cyclic prefix length 
%                  ('Normal'(default),'Extended')
%   Windowing    - Optional. The number of time-domain samples over which 
%                  windowing and overlapping of OFDM symbols is applied. 
%                  (default depends on parameters above, see <a href="matlab: doc('lteOFDMModulate')">lteOFDMModulate</a>
%                  for details)
%   GRID is an MxNxP array where M is the number of subcarriers, N is the
%   number of OFDM symbols and P is the number of antennas. Dimension M
%   must be an integer multiple of 12 RE per RB (NRB=M/12), up to a maximum
%   of 2048. Dimension N must be a multiple of number of symbols in a
%   subframe L, where L=14 for normal cyclic prefix and L=12 for extended
%   cyclic prefix. 
%
%   WAVEFORM is a T-by-P matrix where P is the number of antennas and T is
%   the number of the time domain samples. T=K*30720/2048*Nfft where Nfft
%   is the IFFT size and K is the number of subframes in the input GRID.
%   Nfft is a function of the Number of Resource Blocks (NRB):
%      NRB     Nfft
%        6      128
%       15      256
%       25      512
%       50     1024
%       75     2048
%      100     2048
%   In general, Nfft is the smallest power of 2 greater than or equal to
%   12*NRB/0.85, i.e. the smallest FFT that spans all subcarriers and
%   results in a bandwidth occupancy (12*NRB/Nfft) of no more than 85%.
%   INFO is a structure containing information about the OFDM modulated 
%   waveform. This structure contains the fields:
%   SamplingRate        - The sampling rate of the time domain WAVEFORM, 
%                         given by SamplingRate = 30.72MHz / 2048 * Nfft.
%   Nfft                - The number of FFT points as defined above. 
%   Windowing           - The number of time-domain samples over which 
%                         windowing and overlapping of OFDM symbols is 
%                         applied. 
%   CyclicPrefixLengths - Cyclic prefix length (in samples) of each OFDM 
%                         symbol in a subframe.
%   Note that GRID can span multiple subframes and windowing and
%   overlapping is applied between all adjacent OFDM symbols, including the
%   last of one subframe and the first of the next. Therefore a different
%   result is obtained than if lteOFDMModulate(...) is called on individual
%   subframes and then those time-domain waveforms concatenated. The
%   resulting waveform in that case would have discontinuities at the
%   start/end of each subframe. Therefore it is recommended that all
%   subframes for OFDM modulation first be concatenated prior to calling
%   lteOFDMModulate(...) on the resulting multi-subframe array. However,
%   individual subframes can be OFDM modulated and the resulting
%   multi-subframe time-domain waveform created by manually overlapping. 
%   If ENB.Windowing is absent, a default value for the number of windowed
%   and overlapped samples is used, chosen as a function of NRB to
%   compromise between the effective duration of cyclic prefix (and
%   therefore the channel delay spread tolerance) and the spectral
%   characteristics of the transmitted signal (not considering any
%   additional FIR filtering). See <a href="matlab:
%   doc('lteOFDMModulate')">lteOFDMModulate</a> for details. The value used
%   will be returned in INFO.Windowing. Note that with a value of zero, the
%   issues above concerning concatenation of subframes before OFDM
%   modulation do not apply.
%   When INFO.Nfft=2048, INFO.CyclicPrefixLengths is
%   [160 144 144 144 144 144 144 160 144 144 144 144 144 144] for normal
%   cyclic prefix and [512 512 512 512 512 512 512 512 512 512 512 512] for
%   extended cyclic prefix. For other values of INFO.Nfft, these lengths
%   are scaled by INFO.Nfft/2048.
%
%   [WAVEFORM,INFO] = lteOFDMModulate(ENB,GRID,WINDOWING) allows control of
%   the number of windowed and overlapped samples used in the time-domain
%   windowing, specified by the WINDOWING parameter; the value in
%   ENB.Windowing, if present, is ignored and the output INFO.Windowing
%   will equal WINDOWING.
%
%   Example:
%   % OFDM modulation of one subframe of random uniformly-distributed 
%   % noise, using a 10MHz 2-antenna configuration.
%   
%   enb = struct('NDLRB',50,'CyclicPrefix','Normal','CellRefP',2);
%   dims = lteDLResourceGridSize(enb);
%   reGrid = reshape(lteSymbolModulate(randi([0,1],prod(dims)*2,1), ...
%            'QPSK'),dims);
%   waveform = lteOFDMModulate(enb,reGrid);
%
%   See also lteOFDMDemodulate, lteOFDMInfo, lteDLResourceGrid,
%   lteFadingChannel, lteHSTChannel, lteMovingChannel.

%   Copyright 2009-2020 The MathWorks, Inc.

function [waveform,info,freqGrid] = lteOFDMModulate_fullfreq(enb,grid,varargin)  %grid:时频栅格 72*140

    if (~isnumeric(grid))
        error('lte:error','The input resource grid must be a numeric array.');
    end

    % establish some dimensionality information
    nSC = size(grid,1);
    if (nSC > 2048)
        error('lte:error','The number of subcarriers (%d) cannot be greater than 2048.',nSC);
    end
    % calculate the effective number of RB so that we can get the CP lengths
    % and FFT sizes out of lteOFDMInfo
    log2nsc = log2(nSC);
    if (log2nsc == fix(log2nsc) && log2nsc > 6)
        nrb = fix(0.85*nSC/12);
    else   
        nrb = fix(nSC/12);
        nSC = 12*nrb;
    end
    if (nrb < 6)
        nrb = 6;
    end
    if nrb > 110
        nrb = 110;
    end
    %nrb: RB的个数，nSC：子载波数  （12倍的关系）
    enb.NDLRB = nrb;
    % Handle optional windowing argument.
    if (nargin==3)
        enb.Windowing = varargin{1};
    end
    info = lteOFDMInfo(enb);
    nFFT = double(info.Nfft);
    info.Windowing = double(info.Windowing);
    % N is the number of windowed samples.
    N = info.Windowing;
    nAnts = size(grid,3); %天线数
    cpLengths = double(info.CyclicPrefixLengths);
    symbolsPerSlot = numel(cpLengths)/2;
    nSubframes = floor(size(grid,2)/(symbolsPerSlot*2));
    
    if (N>(nFFT-info.CyclicPrefixLengths(1)))
        error('lte:error','For the Windowing parameter the value (%d) must be less than or equal to %d (the IFFT size (%d) minus the longest cyclic prefix length (%d))',N,nFFT-info.CyclicPrefixLengths(1),nFFT,info.CyclicPrefixLengths(1));
    end
    
    % index of first subcarrier in IFFT input
    firstSC = (nFFT/2) - (nSC/2) + 1; %128/2-72/2+1=29
        
    % number of active subcarriers in IFFT
    if (~any(size(grid,1)==[0 nSC info.Nfft]))
        error('lte:error','The input resource grid must contain a whole number of resource blocks i.e. number of rows must be an integer multiple of 12.');
    end
    
    % validate number of OFDM symbols in the input
    if (mod(size(grid,2),symbolsPerSlot*2)~=0)
        error('lte:error','The input resource grid must contain a whole number of subframes i.e. number of columns must be an integer multiple of 14 for normal cyclic prefix or 12 for extended cyclic prefix.');
    end
    
    % handle various empty input cases
    if(size(grid,1)==0)
        nSubframes = 0;
    end
    if (nSubframes==0)
       head = 0;
       N = 0;
    end
    
    % pre-allocate output  预定义输出时域waveform的维度（全0）
    samplesPerSubframe = sum(cpLengths) + (nFFT*symbolsPerSlot*2);  %每个帧的采样点数
    waveform = zeros(nSubframes*samplesPerSubframe,nAnts);   %*帧数*天线
    
    % pre-calculate windowing; there are two different windows required:
    % one for the first OFDM symbol of a slot and one for other OFDM
    % symbols, because the first OFDM symbol of a slot has a different
    % cyclic prefix length.    
    window0 = lte.internal.raisedCosineWindow(nFFT+cpLengths(1),N);  %对于加长CP的窗
    window1 = lte.internal.raisedCosineWindow(nFFT+cpLengths(2),N);  %对于正常CP的窗
    
    pos = 0;
    % for each OFDM symbol:
    %初始化IFFT前频域栅格：
    freqGrid=zeros(nFFT,nSubframes*symbolsPerSlot*2);
    for i = 1:nSubframes*symbolsPerSlot*2

        % create IFFT input (map subcarriers)
        if (size(grid,1)==nFFT)  %若无空闲子载波（频域值=IFFT数）
            ifftin = squeeze(grid(:,i,:));
        else 
            ifftin = zeros(nFFT,nAnts); % (IFFT的两边填0 )
            ifftin(firstSC+(0:nSC/2-1),:) = grid(1:nSC/2,i,:);  %原频域的前一半值 赋值给IFFT的中间前一半
            ifftin(firstSC+nSC/2+1+(0:nSC/2-1),: ) = grid(nSC/2+1:end,i,:);  %原频域的后一半值 赋值给IFFT的中间后一半
        end

        % perform IFFT
        freq=fftshift(ifftin,1); %freq是ifft前频域的值
        iffout = ifft(freq); %iffout就是OFDM的时域

        freqGrid(:,i)=freq;%构建全频率栅格


        % add cyclic prefix  +CP
        cpLength = cpLengths(mod(i-1,length(cpLengths))+1);
        extended = [iffout(end-(cpLength+N)+1:end,:); iffout];  % CP+OFDM

        % perform windowing, using the appropriate window (first OFDM
        % symbol of a slot or otherwise)
        if (mod(i-1,symbolsPerSlot)==0)
            windowed = extended .* window0;
        else
            windowed = extended .* window1;
        end

        % perform overlapping and creation of output signal. Note that with
        % windowing the signal "head" gets chopped off the start of the
        % waveform and finally superposed on the end. This means the
        % overall signal can be seamlessly looped when output from an
        % arbitrary waveform generator.
        if (i==1)
            % first OFDM symbol: chop the "head" and then output the rest
            head = windowed(1:N,:);
            L = cpLength + nFFT;
            waveform(1:L,:) = windowed(N+1:end,:);
        else
            % subsequent OFDM symbols: add the windowed part to the end of
            % the previous OFDM symbol (overlapping them) and then output
            % the rest; 'pos' points to the end of the previous OFDM symbol
            L = cpLength + nFFT + N;
            waveform(pos-N+(1:L),:) = waveform(pos-N+(1:L),:) + windowed;
        end

        % Update 'pos' to point to the end of the current OFDM symbol
        pos = pos + cpLength + nFFT;
       
    end

    % finally, overlap the "head" with the very end of the signal
    waveform(end-N+1:end,:) = waveform(end-N+1:end,:) + head;
    
end
