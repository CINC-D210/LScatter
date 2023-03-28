function cbIdx = lteCbSelect(h, enPMIfback, numTx, numLayers, snrdB)
%#codegen

% Copyright 2012 The MathWorks, Inc.

% LTE Codebook selection using:
%   maximum capacity, or minimum MSE, or minimum singular value criteria.
%
% Get all Wn, For each Wn, calculate metric value and sum over all k
% Pick Wn for max/min sum based on criteria, send back cbIdx.
%
% For now, MMSE is the active criteria, 
%   uncomment lines below to switch to one of other two criteria.

if (enPMIfback)
    
    snr = 10^(0.1*snrdB);    
    numk = size(h, 1);
    
    if (numTx == 2)
        cbLen = 2; % Only indices 1 and 2 are used for 2-layer closed-loop Spatial MUX
%         Ikcb = zeros(cbLen, 1);
        MSEcb = zeros(cbLen, 1);
%         SVcb = zeros(cbLen, 1);
        
        for cbIdx = 1:cbLen
            Wn = getWn(cbIdx, numLayers, numTx);
            
            for idx = 1:numk
                hk = h(idx, :, :);                            % numTx x numRx
                hk = reshape(hk(:), numLayers, numLayers).';  % numRx x numTx
                
%                 currI = log2( real( det(eye(numLayers) + (snr/numLayers)*(Wn'*hk'*hk*Wn)) ) );
%                 Ikcb(cbIdx) = Ikcb(cbIdx) + currI;
                currMSE = real(trace((snr/numLayers)*pinv(eye(numLayers) + (snr/numLayers)*(Wn'*(hk'*hk)*Wn)))); 
                MSEcb(cbIdx) = MSEcb(cbIdx) + currMSE;
%                 currSV = svd(hk*Wn);
%                 SVcb(cbIdx) = SVcb(cbIdx) + min(currSV);
            end
        end
%         if (max(Ikcb)-min(Ikcb) > sqrt(eps))
%             [~, cbIdx] = max(Ikcb); % 0-based, note 0 and 3 are not used
%         else
%             cbIdx = 1;
%         end
        if (max(MSEcb)-min(MSEcb) > sqrt(eps))
            % Only take minimum if significantly different
            [~, cbIdx] = min(MSEcb); % 0-based, note 0 and 3 are not used
        else % no real difference
            % Keep the same cbIdx as before, i.e. no change
            cbIdx = 1;
        end
%         if (max(SVcb)-min(SVcb) > sqrt(eps))
%             [~, cbIdx] = max(SVcb); % 0-based, note 0 and 3 are not used
%         else
%             cbIdx = 1;
%         end

    else % for numTx=4
        cbLen = 2^numLayers;
        
%         Ikcb = zeros(cbLen, 1);
        MSEcb = zeros(cbLen, 1);
%         SVcb = zeros(cbLen, 1);
        for cbIdx = 1:cbLen
            Wn = getWn(cbIdx-1, numLayers, numTx);
            
            for idx = 1:numk
                hk = h(idx, :, :);                            % numTx x numRx
                hk = reshape(hk(:), numLayers, numLayers).';  % numRx x numTx
                
%                 currI = log2( real( det(eye(numLayers) + (snr/numLayers)*(Wn'*hk'*hk*Wn)) ) );
%                 Ikcb(cbIdx) = Ikcb(cbIdx) + currI;
                currMSE = real(trace((snr/numLayers)*pinv(eye(numLayers) + (snr/numLayers)*(Wn'*(hk'*hk)*Wn)))); 
                MSEcb(cbIdx) = MSEcb(cbIdx) + currMSE;
%                 currSV = svd(hk*Wn);
%                 SVcb(cbIdx) = SVcb(cbIdx) + min(currSV);
            end
        end
        
%         if (max(Ikcb)-min(Ikcb) > sqrt(eps))
%             [~, cbIdx] = max(Ikcb); % 1-based
%             cbIdx = cbIdx-1;        % 0-based
%         else
%             cbIdx = 1;
%         end
        if (max(MSEcb)-min(MSEcb) > sqrt(eps))
            % Only take minimum if significantly different
            [~, cbIdx] = min(MSEcb); % 1-based
            cbIdx = cbIdx-1;        % 0-based
        else % no real difference
            % Keep the same cbIdx as before, i.e. no change
            cbIdx = 1;
        end
%         if (max(SVcb)-min(SVcb) > sqrt(eps))
%             [~, cbIdx] = max(SVcb); % 1-based
%             cbIdx = cbIdx-1;        % 0-based
%         else
%             cbIdx = 1;
%         end
    end
    
else
    cbIdx = 1;
end

%--------------------------------------------------------------------------
function Wn = getWn(cbIdx, numLayers, numTx)
% Get Wn for a specific codebook index.
%   Similar code as the precoder.

j = complex(0,1);

switch numTx
    case 2  % for numLayers = 2 only
        Wn = complex(ones(numTx, numLayers));
        switch cbIdx
            case 1
                Wn = (1/2)*[1 1; 1 -1];
            case 2
                Wn = (1/2)*[1 1; j -j];
        end
    case 4
        un = complex(ones(4, 1));
        switch cbIdx
            case 0
                un = [1 -1 -1 -1].';
            case 1
                un = [1 -j 1 j].';
            case 2
                un = [1 1 -1 1].';
            case 3
                un = [1 j 1 -j].';
            case 4
                un = [1 (-1-j)/sqrt(2) -j (1-j)/sqrt(2)].';
            case 5
                un = [1 (1-j)/sqrt(2) j (-1-j)/sqrt(2)].';
            case 6
                un = [1 (1+j)/sqrt(2) -j (-1+j)/sqrt(2)].';
            case 7
                un = [1 (-1+j)/sqrt(2) j (1+j)/sqrt(2)].';
            case 8
                un = [1 -1 1 1].';
            case 9
                un = [1 -j -1 -j].';
            case 10
                un = [1 1 1 -1].';
            case 11
                un = [1 j -1 j].';
            case 12
                un = [1 -1 -1 1].';
            case 13
                un = [1 -1 1 -1].';
            case 14
                un = [1 1 -1 -1].';
            case 15
                un = [1 1 1 1].';
        end
        Wn = eye(4) - 2*(un*un')./(un'*un);
        switch cbIdx    % order columns, for numLayers = 4 only
            case {2, 3, 14}
                Wn = Wn(:, [3 2 1 4]);
            case {6, 7, 10, 11, 13}
                Wn = Wn(:, [1 3 2 4]);
        end
        Wn = Wn./sqrt(numLayers);
end

% [EOF]
