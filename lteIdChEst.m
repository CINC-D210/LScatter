function [hD,H] = lteIdChEst(prmLTEPDSCH, prmMdl, chPathG, nS)
%LTEIDCHEST Ideal channel estimation for LTE subframes.
%
%   Given the system parameters and the MIMO channel path Gains, provide
%   the ideal channel estimates for the RE corresponding to the data.

%   Copyright 2012 The MathWorks, Inc.

% Pseudo-logic
% Option 1: Implemented below
%   Limitation - will work for path delays to be multiple of channel sample
%   time and largest pathDelay < nfft
%       add interpolation for unrestricted use
%
%   Given the 1920 x numPaths x Tx x Rx pathaGains, convert these
%   into single tap gains by way of OFDM CP
%   i.e. 1920 x Tx x Rx - in size
%   Per link processing
%       1920 x numPaths -> take mean over each OFDM symbol-> get 14 x numPaths,
%           acct for CP, take FFT, reorder and scale
%   Then similar to OFDM Rx and ExtData, get the Data RE based tap Gains
%   And use these going forward
%       size(chPathG)   1920 7 2 2
%       rxFade     1920x2           
%       rxSig      1920x2           
%       rxGrid     72x14x2         
%       chEst      912x2x2           
%       dataRx     648x2            
%       hD         648x2x2 
%
% Option 2: future consideration
%   Get channel gains, fft & average
%   For a 2x2 scheme first - try all channel EPA, EVA, ETU options
%   Work off pilots and interpolate between them - provide options

persistent hFFT; 
if isempty(hFFT) 
   hFFT = dsp.FFT; 
end 

% get parameters
numDataTones = prmLTEPDSCH.Nrb*12; % Nrb_sc = 12
N            = prmLTEPDSCH.N;
cpLen0       = prmLTEPDSCH.cpLen0;
cpLenR       = prmLTEPDSCH.cpLenR;

slotLen = (N*7 + cpLen0 + cpLenR*6);

if strncmp(prmMdl.chanMdl, 'Fre', 3)
    pathDelays = 0;
elseif strncmp(prmMdl.chanMdl, 'EPA', 3)
    pathDelays = [0 30 70 90 110 190 410]*1e-9;
elseif strncmp(prmMdl.chanMdl, 'EVA', 3)
    pathDelays = [0 30 150 310 370 710 1090 1730 2510]*1e-9;
elseif strncmp(prmMdl.chanMdl, 'ETU', 3)
    pathDelays = [0 50 120 200 230 500 1600 2300 5000]*1e-9;
elseif strncmp(prmMdl.chanMdl, 'f', 1)
    pathDelays = [0 10 20 30 100]*(1/prmLTEPDSCH.chanSRate);
elseif strncmp(prmMdl.chanMdl, 'User', 4)
    pathDelays = prmMdl.pathDelays *(1/prmLTEPDSCH.chanSRate);
end
% Delays, in terms of number of channel samples, +1 for indexing
sampIdx = round(pathDelays/(1/prmLTEPDSCH.chanSRate)) + 1;

[~, numPaths, numTx, numRx] = size(chPathG);

H = complex(zeros(numDataTones, 14, numTx, numRx));
for i= 1:numTx
    for j = 1:numRx
        link_PathG = chPathG(:, :, i, j);
        % Split this per OFDM symbol
        g = complex(zeros(2*7, numPaths));
        for jj = 1:2 % over two slots
            % First OFDM symbol
            g((jj-1)*7+1, :) = mean(link_PathG((jj-1)*slotLen + (1:(N+cpLen0)), :), 1);
            
            % Next 6 OFDM symbols
            for k = 1:6
                g((jj-1)*7+k+1, :) = mean(link_PathG((jj-1)*slotLen+cpLen0+k*N+(k-1)*cpLenR + (1:(N+cpLenR)), :), 1);
            end
        end
        hImp = complex(zeros(2*7, N));
        hImp(:, sampIdx) = g; % assign pathGains at sample locations
        % FFT processing
        h = step(hFFT, hImp.'); 

        % Reorder, remove DC, Unpack channel gains
        h = [h(N/2+1:N, :); h(1:N/2, :)];
        H(:, :, i, j) = [h(N/2-numDataTones/2+1:N/2, :); h(N/2+2:N/2+1+numDataTones/2, :)];
    end
end
% H - 72x14x2x2

% Now, align these with the data RE per antenna link and reuse them
[H_rx1, csrRx1] = lteExtData( squeeze(H(:,:,:,1)), nS, prmLTEPDSCH, 'chan');
hD =  complex(zeros(size(H_rx1,1), numTx, numRx));
hD(:,:,1) = H_rx1;
csrRx =  complex(zeros(size(csrRx1,1), numTx, numRx));
csrRx(:,:,1) = csrRx1;
for i = 2:numRx
    [hD(:,:,i), csrRx(:,:,i)] = lteExtData( squeeze(H(:,:,:,i)), nS, prmLTEPDSCH, 'chan');
end

% [EOF]
