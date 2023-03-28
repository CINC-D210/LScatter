%lteRMCDLTool Downlink RMC waveform generation
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCDLTool(...) returns WAVEFORM, the
%   generated reference measurement channel waveform, GRID representing the
%   populated resource grid for all the physical channels specified in TS
%   36.101 Annex A.3 and RMCCFGOUT containing information about parameters
%   associated with the generated waveform. The RMC waveform can be
%   configured via a graphical user interface (GUI) or by passing the
%   required input parameters in a function call.
%
%   WAVEFORM is a T-by-P matrix where T is the number of time domain
%   samples and P is the number of antennas. GRID is a 3-dimensional array
%   of resource elements for a number of subframes across all configured
%   antenna ports, as described in the <a href="matlab:web([docroot '/lte/gs/data-structures.html'])">Data Structures</a> documentation.
%   RMCCFGOUT is a structure containing information about the OFDM 
%   modulated waveform as described in <a
%   href="matlab:help('lteOFDMInfo')">lteOFDMInfo</a> and HARQ scheduling
%   in addition to the configuration parameters as described in <a href="matlab:help('lteRMCDL')">lteRMCDL</a>.
%   The internal HARQ scheduling sequence (present in the
%   HARQProcessSequence field of RMCCFGOUT) is optimized according to the
%   transport block sizes, number of HARQ processes, duplex mode and
%   Uplink/Downlink configuration (if duplex mode is TDD).
%
%   lteRMCDLTool launches a GUI for the parameterization, generation and
%   visualization of the RMC waveforms.
%
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCDLTool(RC,TRDATA,DUPLEXMODE,...
%   TOTSUBFRAMES) 
%   returns the WAVEFORM, GRID and RMCCFGOUT for the default reference
%   measurement channel defined by RC (see <a
%   href="matlab:doc('lteRMCDLTool')">lteRMCDLTool</a> for details) using
%   the information bits TRDATA. DUPLEXMODE and TOTSUBFRAMES are optional
%   input parameters and define the duplex mode of the generated waveform
%   and total number of subframes that make up the GRID.
%
%   The input RC identifies the reference measurement channel number as per
%   TS 36.101, must be one of following:
%   ('R.0','R.1','R.2','R.3','R.4','R.5','R.6','R.7','R.8','R.9','R.10',
%   'R.11','R.12','R.13','R.14','R.25','R.26','R.27','R.28','R.31-3A',
%   'R.31-4','R.43','R.44','R.45','R.45-1','R.48','R.50','R.51','R.68-1',
%   'R.105','R.6-27RB','R.12-9RB','R.11-45RB')
%
%   Note that RCs 'R.31-3A', 'R.31-4', 'R.68-1', and 'R.105' are sustained
%   data rate RMCs with user data in subframe 5 and RCs 'R.6-27RB',
%   'R.12-9RB', and 'R.11-45RB' are custom RMCs configured for non-standard
%   bandwidths but with the same code rate as the standardized versions.
%
%   TRDATA is a vector or cell array containing one or two vectors of bit
%   values where each vector is the information bits stream to be coded
%   across the duration of the generation i.e. representing multiple
%   concatenated transport blocks. Internally these vectors are looped if
%   the number of bits required across all subframes of the generation
%   exceeds the length of the vectors provided. This allows for the user to
%   enter a short pattern e.g. [1; 0; 0; 1] that will be repeated as the
%   input to the transport coding. In each subframe of generation, the
%   number of data bits taken from this stream is given by the elements of
%   the TrBlkSizes matrix, a field of the PDSCH substructure of the RMC
%   configuration structure RMCCFGOUT.
%
%   TRDATA can also contain empty vector(s) representing no transport data
%   thus resulting in no PDSCH and its corresponding PDCCH transmission. In
%   other words, the transmission of PDSCH and its corresponding PDCCH can
%   be skipped in the WAVEFORM when the TRDATA contains empty vector(s).
%   The other physical channels and signals are transmitted as normal in
%   generated WAVEFORM.
%
%   DUPLEXMODE is an optional input representing the frame structure type
%   ('FDD'(default),'TDD').
%
%   TOTSUBFRAMES is an optional representing the total number of subframes
%   to be generated (default 10).
%
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCDLTool(RMCCFG,TRDATA) generates the
%   WAVEFORM, GRID and RMCCFGOUT in the same way as above except it takes
%   the user defined reference channel structure RMCCFG (see <a href="matlab:doc('lteRMCDLTool')">lteRMCDLTool</a>
%   for details) as input parameter. The reference configuration structure
%   with default parameters can easily be created with the function 
%   <a href="matlab:help('lteRMCDL')">lteRMCDL</a> which is designed to generate the various RMC configuration 
%   structures as defined in TS 36.101 Annex A.3. This configuration
%   structure can then be modified as required and used in the generation
%   of WAVEFORM. For UE-specific beamforming (i.e. 'Port5', 'Port7-8',
%   'Port8' and 'Port7-14' transmission schemes), the number of transmit
%   antennas is deduced from the precoding matrix W (W is of size
%   NLayers-by-NTxAnts), which is specified in RMCCFG.
%
%   SIB1 messages and the associated PDSCH and PDCCH can be generated if
%   the structure SIB is present in RMCCFG. This structure must include the
%   following fields
%   Data            - SIB1 transport block information bits
%   VRBStart        - Initial virtual RB allocation block
%   VRBLength       - Length of RB allocation
%   Enable          - Optional. Enable SIB generation ('On'(default),'Off')
%   DCIFormat       - Optional. DCI format ('Format1A'(default),'Format1C')
%   AllocationType  - Optional. Localized or distributed allocation 
%                     (0(default),1) respectively
%   Only required for Format1A (DCIFormat = 'Format1A')
%       N1APRB      - Optional. Transport block set selection parameter
%                     (2,3). Indicates the column in Table 7.1.7.2.1-1 of
%                     TS 36.213 for transport block size selection. The
%                     default is the smallest of (2,3) that provides a
%                     valid transport block size bigger than or equal to
%                     the length of the Data field
%   Only required for distributed allocation (AllocationType = 1)
%       Gap         - Distributed allocation gap (0,1) for 1st gap or
%                     2nd gap respectively
%
%   Note that the set of valid transport block sizes is specified in Table
%   7.1.7.2.1-1 of TS 36.213. Only columns 2 and 3 apply to the SIB DL-SCH.
%   The Data field is padded with zeros to the closest valid size from this
%   table.
%
%   The lowest order information bit of the SIB Data field should be mapped
%   to the most significant bit of the SIB1 transport block as defined in
%   section 6.1.1 of TS 36.321.
%
%   As per TS 36.101 Annex A.3, a reference PDSCH transmission will never
%   be scheduled in subframe 5 except for the SIB1 associated PDSCH
%   mentioned above.
%
%   The values for CFI and PRBSet can vary per subframe. If these
%   parameters are arrays then the function cyclically steps through the
%   elements of the array starting with index given by
%   mod(NSubframe,length(parameter)). In the case of PRBSet the parameter
%   will have to be a cell array of column vectors or slot-wise matrices.
%
%   Note that the PHICH symbols carry a single ACK on the first PHICH 
%   instance in each PHICH group.
%
%   Example 1:
%   % Generate a time domain signal txWaveform and a 3-dimensional
%   % array of the resource elements txGrid for R.12 as specified 
%   % in TS 36.101. This transmission will use 16QAM modulation 
%   % scheme instead of QPSK.
%   
%   rmc.RC = 'R.12';
%   rmc.PDSCH.Modulation = '16QAM';
%   [txWaveform,txGrid,rmcCfgOut] = lteRMCDLTool(rmc,[1;0;0;1]);
%
%   Example 2:
%   % Generate a time domain signal txWaveform and a 3-dimensional
%   % array of the resource elements txGrid for R.31-4 FDD as specified 
%   % in TS 36.101 Annex A.3.9.1-1. R.31-4 FDD is 20MHz, 64QAM, variable 
%   % code rate and has user data scheduled in subframe 5:
%   
%   [txWaveform,txGrid,rmcCfgOut] = lteRMCDLTool('R.31-4',{[1;0] [1;0]});
%
%   Example 3:
%   % Enable SIB transmission using DCI format 1A and localized
%   % allocation. 
%   
%   rmc = lteRMCDL('R.3');
%   rmc.SIB.Enable = 'On'; 
%   rmc.SIB.DCIFormat = 'Format1A';
%   rmc.SIB.AllocationType = 0; % 0 localized
%   rmc.SIB.VRBStart = 8;
%   rmc.SIB.VRBLength = 8;
%   rmc.SIB.Data = randi([0 1],144,1);
%   [txWaveform,txGrid,~] = lteRMCDLTool(rmc,[1;0;0;1]);
%
%   See also lteRMCDL, lteDLConformanceTestTool, lteRMCULTool,
%   lteTestModelTool.

%   Copyright 2009-2020 The MathWorks, Inc.

function [timeDomainSig,frame,rmc,freqGrid] = lteRMCDLTool(varargin)

% Preserve global random number generator state when function exits 
s = rng;
c = onCleanup(@()rng(s));
% localSeed initializes the random number generator seeds used for OCNG generation
localSeed = 0; 

if(isempty(varargin))
    wirelessWaveformGenerator('Downlink RMC');
else
    if((ischar(varargin{1}) || isstring(varargin{1})) && (nargin >= 2))
        % RMC number
        rmcNo = varargin{1};
       
        % Input transport data
        trData = varargin{2};
        
        % Acquiring "TotSubframes" if supplied
        if(nargin == 4)
            totSubframes = varargin{4};
        else
            totSubframes = 10;
        end
        
        % Deducing number of codewords from input transport block size
        if(iscell(trData))
            [rows,cols] = size(trData);
            Ncodewords = max(rows,cols);            % Number of codewords
        elseif(isvector(trData) || isempty(trData))
            Ncodewords = 1;
        end
        if(Ncodewords == 0 || isempty(trData))
            Ncodewords = 1;
        end
        
        % Acquiring "DuplexMode" if supplied
        if(nargin > 2)
            enbConfig.DuplexMode = varargin{3};
        end
        
        % Generating config structure for given RC
        enbConfig.RC = rmcNo;
        enbConfig.TotSubframes = totSubframes;
        enbConfig = lteRMCDL(enbConfig,Ncodewords);
        PDSCH = enbConfig.PDSCH;
    elseif(isstruct(varargin{1}) && (nargin == 2))
        % Input configuration structure
        enbConfig = varargin{1};
        % Input transport data
        trData = varargin{2};                   
        
        % Deducing number of codewords from input transport block size
        if(iscell(trData))
            [rows,cols] = size(trData);
            Ncodewords = max(rows,cols);            % Number of codewords
        elseif(isvector(trData)|| isempty(trData))
            Ncodewords = 1;
        end
        if(Ncodewords == 0 || isempty(trData))
            Ncodewords = 1;
        end
        
        % Checking if input config structure is missing any parameter
        completeStruct = checkFields(enbConfig);
        
        % Defaulting parameter values if the input config structure is not
        % complete
        if(~completeStruct)
            enbConfig = lteRMCDL(enbConfig,Ncodewords);
        end

        PDSCH = enbConfig.PDSCH;
        % For a single codeword transmission, ensure that there is only
        % one modulation scheme specified, else ltePDSCHIndices would
        % treat it as a 2 codeword transmission.
        if(numel(PDSCH.Modulation)==2 && Ncodewords==1)
            PDSCH.Modulation = PDSCH.Modulation(1);
        end
    else
        if nargin==1
            if ischar(varargin{1}) || isstring(varargin{1})
                lteRMCDL(varargin{1}); % Check if RMC is valid
            end
            error('lte:error','The function call resulted in an error: Information bit stream not specified.');
        else
            error('lte:error','The function call resulted in an error: Unrecognized function signature.');
        end
    end

    % Setting up data generator source(s).
    if(iscell(trData))
        idx = cellfun(@isempty,trData);
        if(isempty(trData))
            dataSource1=saVectorDataSource([]);
            dataSource2=saVectorDataSource([]);
            PDSCH.TrBlkSizes(:) = 0;
            PDSCH.CodedTrBlkSizes(:) = 0;
        else
            dataSource1=saVectorDataSource(trData{1});
            if(numel(trData)>1)
                dataSource2=saVectorDataSource(trData{2});
            elseif(idx)
                PDSCH.TrBlkSizes(:) = 0;
                PDSCH.CodedTrBlkSizes(:) = 0;
            end
        end
    else
        if(isempty(trData))
            PDSCH.TrBlkSizes(:) = 0;
            PDSCH.CodedTrBlkSizes(:) = 0;
        end
        dataSource1=saVectorDataSource(trData);
    end
    
    frame = [];
    reGrid = [];
        
    % Setting default value to SerialCat if its not defined.
    if(~isfield(enbConfig,'SerialCat'))
        enbConfig.SerialCat = true;
    end
    serialCat = enbConfig.SerialCat;
    
    if ~any(strcmpi(enbConfig.OCNGPDSCHEnable,{'On','Off'}))
        error('lte:error','The function call resulted in an error: The valid values for OCNGPDSCHEnable are ''On'' and ''Off''');
    end
    if ~(isnumeric(enbConfig.OCNGPDSCHPower) && isscalar(enbConfig.OCNGPDSCHPower))
        error('lte:error','The function call resulted in an error: OCNGPDSCHPower must be a numeric scalar value');
    end
    if ~any(strcmpi(enbConfig.OCNGPDCCHEnable,{'On','Off'}))
        error('lte:error','The function call resulted in an error: The valid values for OCNGPDCCHEnable are ''On'' and ''Off''');
    end
    if ~(isnumeric(enbConfig.OCNGPDCCHPower) && isscalar(enbConfig.OCNGPDCCHPower))
        error('lte:error','The function call resulted in an error: OCNGPDCCHPower must be a numeric scalar value');
    end
    
    % Empty subframe construction for given configuration filled with all
    % zeros.
    griddims=lteDLResourceGridSize(enbConfig);
    % For the UE specific beamforming, size the output grid for maximum
    % dimensions. Depending on the specific configuration, the grid's 3rd
    % dimension can be CellRefP, NTxAnts, NLayers or CSIRefP.
    if (isUeRS(PDSCH.TxScheme))
        if(isfield(PDSCH,'W'))
            PDSCH.NTxAnts = size(PDSCH.W,2);
            if (PDSCH.NTxAnts==0)
                griddims(3)= max([PDSCH.NLayers,griddims(3)]) ;
            else
                griddims(3)= max([PDSCH.NTxAnts,griddims(3)]) ;
            end
            
            if (isfield(enbConfig,'CSIRefP') && (enbConfig.CSIRefP>griddims(3)))
                if ~((isfield(enbConfig,'CSIRSPeriod') && strcmpi(enbConfig.CSIRSPeriod,'Off')))
                    griddims(3) =  enbConfig.CSIRefP;
                end
            end
        end
    end
    subframe = zeros(griddims);
    
    % Getting total number of subframes to be generated
    totalSubframes = enbConfig.TotSubframes;
    
    % Validate total number of subframes
    if(~isnumeric(totalSubframes) || isempty(totalSubframes) || totalSubframes<0)
        error('lte:error','The function call resulted in an error: The total number of subframes must be a positive integer');
    end
    
    % Setting modulation format
    if(~iscell(PDSCH.Modulation) && ischar(PDSCH.Modulation))
        PDSCH.Modulation = {PDSCH.Modulation};
    end
    
    % Setting transport channel parameters to be used in channel coding.   
    if(~isfield(PDSCH,'RVSeq') || (isfield(PDSCH,'RVSeq') && isempty(PDSCH.RVSeq)))
        error('lte:error','The function call resulted in an error: Could not find a structure field called RVSeq or it is empty.');
    end
	
    % HARQ setup initialization
    if(~isfield(PDSCH,'NHARQProcesses'))
        error('lte:error','The function call resulted in an error: Could not find a structure field called NHARQProcesses.');
    end
    harqTable = getHARQTable(enbConfig);
    harqprocess.data = struct('blk1',[],'blk2',[]);
    harqprocess.RVIdx = ones(Ncodewords,1);
    harqprocess.NDI = [false false];
    harqProcesses(1:max(harqTable)) = harqprocess;
    newData = [1 1];
    
    % Set up for correct coded transport block sizes as per available
    % resources
    cBlkSizeIdx=0;
    PDSCH.CodedTrBlkSizes=0;
    
    [nRows,nColms] = size(PDSCH.TrBlkSizes);
    tempBlkSize = PDSCH.TrBlkSizes;
    
    if(Ncodewords==2)
        if(nColms==2)
            tempBlkSize = tempBlkSize.';
        end
    else
        if(nRows>2)
            tempBlkSize = tempBlkSize.';
        end
    end

    if isempty(enbConfig.OCNGPDSCH.RNTI)
        % If PDSCH OCNG RNTI is specified as empty, use 0
        enbConfig.OCNGPDSCH.RNTI = 0;
    end
    
    % Input configuration structure is passed at output with NTxAnts
    % deduced from W so that the waveform can be regenerated with the
    % function output configuration and the low-level functions requiring
    % NTxAnts have the value defined correctly.
	rmc = enbConfig;
    if (isfield(PDSCH,'NTxAnts') && isfield(rmc.PDSCH,'NTxAnts'))
        rmc.PDSCH.NTxAnts = PDSCH.NTxAnts;
    end
    
    % Store CFI to extract values per subframe if it is a vector
    cfi = enbConfig.CFI;
    % Store PRBSet to extract values per subframe if it is a a cell array
    prbset = PDSCH.PRBSet;
        
    % Get the absolute subframe number from NSubframe and NFrame
    NSubframe = enbConfig.NFrame*10+enbConfig.NSubframe;
    
    for subframeIdx=NSubframe:(NSubframe+totalSubframes)-1
        
        % Update subframe number and clear subframe
        subframe(:)=0;
        enbConfig.NSubframe = mod(subframeIdx,10);
        enbConfig.NFrame = floor(subframeIdx/10);

        harqIdx = harqTable(mod(subframeIdx,length(harqTable)) + 1);
        
        info=lteDuplexingInfo(enbConfig);
        cBlkSizeIdx = cBlkSizeIdx + 1;
        if (info.NSymbolsDL)
            % This function provides transport block sizes and coded transport
            % block sizes for given RMC and subframe number, as specified in TS
            % 36.101 A.3
            if(Ncodewords==2)
                transportBlkSize = tempBlkSize(:,mod(enbConfig.NSubframe,size(tempBlkSize,2)) + 1);
                if(size(transportBlkSize,1)==2)
                    transportBlkSize = transportBlkSize.';
                elseif(size(transportBlkSize,1)==1 && size(transportBlkSize,2)==1)
                    transportBlkSize = [transportBlkSize transportBlkSize]; %#ok<AGROW>
                end
            else
                transportBlkSize = tempBlkSize(1,mod(enbConfig.NSubframe,size(tempBlkSize,2)) + 1);
            end
            
            % Indices needs to be generated to map the modulated symbols on resource
            % grid. PDSCH indices generator is being used to generate indices to
            % map PDSCH symbols on resource grid.
            if iscell(prbset)
                % If varying PRBSet, choose the value corresponding to
                % the active subframe
                PDSCH.PRBSet = prbset{mod(subframeIdx,numel(prbset))+1};
                if ~isempty(PDSCH.PRBSet) && (size(PDSCH.PRBSet,2)~=1) && (size(PDSCH.PRBSet,2)~=2)
                    error('lte:error',"The contents of cell number " + ...
                        string(mod(subframeIdx,numel(prbset))+1) + ...
                        " in PDSCH.PRBSet must be a column vector or a two-column matrix.");
                end
            else
                if ~isempty(prbset) && (size(prbset,2)~=1) && (size(prbset,2)~=2)
                    error('lte:error','The field PDSCH.PRBSet must be a column vector, a two-column matrix or a cell array.');
                end
            end
            prbs = PDSCH.PRBSet;
            
            % Get the correct CFI value according to subframe
            enbConfig.CFI = cfi(mod(subframeIdx,numel(cfi))+1);
            % These indices are made '1based' for direct mapping on resource grid
            [pdschIndices,pdschinfo] = ltePDSCHIndices(enbConfig,PDSCH,prbs,{'1based'});
            
            % If Ncodewords is two and modulation is a single element cell
            % array then only pick the first element of vector G to treat
            % this as a one codeword transmission.
            if(numel(pdschinfo.G)==Ncodewords)
                % Get the DLSCH capacities per active codeword 
                codedTrBlkSize = pdschinfo.G;
            else
                codedTrBlkSize = pdschinfo.G(1);
            end
            
            if (Ncodewords==2) && (sum(idx)==1) 
                codedTrBlkSize = getcodedTrBlkSize(PDSCH,pdschinfo.Gd,find(idx==0));
            end
            codedTrBlkSize(transportBlkSize==0) = 0;
            PDSCH.CodedTrBlkSizes(1:Ncodewords,cBlkSizeIdx) = codedTrBlkSize;

            % Updating RV sequence
            if(harqIdx>0)
                tempRv(1,1) = PDSCH.RVSeq(1,harqProcesses(harqIdx).RVIdx(1));
                if(Ncodewords==2)
                    tempRv(2,1) = PDSCH.RVSeq(end,harqProcesses(harqIdx).RVIdx(2));
                end
                PDSCH.RV = tempRv(:,1);
                newData(:) = (harqProcesses(harqIdx).RVIdx == 1);
                % Toggle the NDI bits when the RV index goes back to 1
                harqProcesses(harqIdx).NDI(harqProcesses(harqIdx).RVIdx == 1) = ~harqProcesses(harqIdx).NDI(harqProcesses(harqIdx).RVIdx == 1);
                harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) = mod(harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) + 1,size(PDSCH.RVSeq,2)+1);
                harqProcesses(harqIdx).RVIdx(harqProcesses(harqIdx).RVIdx == 0) = 1;
                
                % DL-SCH transport block size of configured RMC as in TS 36.101.
                if(Ncodewords==1 && (newData(1) || isempty(harqProcesses(harqIdx).data.blk1)))
                    harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize);
                elseif(Ncodewords == 2)
                    transportBlkSize(idx) = 0;
                    codedTrBlkSize(idx) = 0;
                    if(isempty(PDSCH.RV))
                        PDSCH.RV = 0;
                    end
                    if(newData(1) || isempty(harqProcesses(harqIdx).data.blk1))
                        harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize(1));
                    end
                    if(newData(2) || isempty(harqProcesses(harqIdx).data.blk2))
                        harqProcesses(harqIdx).data.blk2 = dataSource2.getData(transportBlkSize(2));
                    end
                end
                if(transportBlkSize(1)~=0)
                    dlschTransportBlk1 = harqProcesses(harqIdx).data.blk1;
                else
                    dlschTransportBlk1 = [];
                end
                if(Ncodewords == 2 && transportBlkSize(2)~=0)
                    dlschTransportBlk2 = harqProcesses(harqIdx).data.blk2;
                else
                    dlschTransportBlk2 = [];
                end
            else
                dlschTransportBlk1 = [];
                dlschTransportBlk2 = [];
            end
            
            % Generating DL-SCH coded bits by performing complete channel coding including
            % CRC calculation, code block segmentation and CRC attachment, turbo
            % coding, rate matching and code block concatenation.
            trData = {dlschTransportBlk1};
            if Ncodewords==2
                 trData{1,2} = dlschTransportBlk2;
            end
            codedTrBlock = lteDLSCH(enbConfig,PDSCH,codedTrBlkSize,trData);

            if(~all(codedTrBlkSize==0) && (harqIdx>0))
                
                if (isUeRS(PDSCH.TxScheme))
                    if(~isfield(PDSCH,'W'))
                        error('lte:error','The function call resulted in an error: For UE specific beamforming the precoding matrix W must be present.');
                    end
                end
                
                if(Ncodewords==2 && isempty(codedTrBlock{1}) && ~isempty(codedTrBlock{2}))
                    % Swapping codeword and corresponding modulation scheme
                    codedTrBlock{1} = codedTrBlock{2};
                    codedTrBlock{2} = [];
                    tempMod = PDSCH.Modulation;
                    PDSCH.Modulation{1} = tempMod{2};
                    PDSCH.Modulation{2} = tempMod{1};

                    % Complex-valued modulated symbol generation for PDSCH. This involves
                    % scrambling, modulation, layering and precoding processes.
                    pdschSymbols = ltePDSCH(enbConfig,PDSCH,codedTrBlock);

                    % Setting the modulation scheme to original order
                    PDSCH.Modulation = tempMod;
                else
                    % Complex-valued modulated symbol generation for PDSCH. This involves
                    % scrambling, modulation, layering and precoding processes.
                    pdschSymbols = ltePDSCH(enbConfig,PDSCH,codedTrBlock);
                end
                
                % PDSCH symbols mapping on resource grid is performed using PDSCH indices.
                powerAdjPerRE = 10^((PDSCH.Rho)/20);
                subframe(pdschIndices) = pdschSymbols*powerAdjPerRE;    
                
                % Transmit UE-specific reference signal (DM-RS) if applicable
                if (isUeRS(PDSCH.TxScheme))
                    ueRSIndices = lteDMRSIndices(enbConfig,PDSCH);
                    ueRSSymbols = lteDMRS(enbConfig,PDSCH);
                    subframe(ueRSIndices) = ueRSSymbols*powerAdjPerRE;    
                end
                
                enbConfig.PDSCH = PDSCH;
                
                % Generating proper DCI message for configured PDSCH
                % transmission.
                [~,dciMsgBits] = getDCI(enbConfig,harqIdx,harqProcesses(harqIdx).NDI,Ncodewords,prbs);
                
                % Calculating PDCCH region resource information, e.g. total number
                % of bits associated with PDCCH (Mtot), total number of CCEs, total 
                % number of REs and REGs.
                pdcchInfo = ltePDCCHInfo(enbConfig);
                
                % The following parameters are required by the DCI channel coding stages: 
                % number of downlink resource blocks, UE-specific mask (16-bit RNTI value) 
                % and PDCCH format. 
                pdcchConfig.NDLRB = enbConfig.NDLRB;  % No of DL-RB in total BW
                pdcchConfig.RNTI = PDSCH.RNTI;        % 16-bit value number
                pdcchConfig.PDCCHFormat = PDSCH.PDCCHFormat; % Aggregation level
                
                % Performing DCI message bits coding to form coded DCI bits.
                codedDciBits = lteDCIEncode(pdcchConfig,dciMsgBits);
                
                % Not all the available bits in the PDCCH region are necessarily used. 
                % Therefore the convention we have adopted is to set unused bits to -1, 
                % while bit locations with values 0 or 1 are used. Initially, all elements 
                % are initialized with -1 to indicate that all the bits are unused.
                pdcchBits = -1*ones(1,pdcchInfo.MTot);

                % Performing search space for UE-specific control channel candidates.
                candidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});

                % Mapping PDCCH payload on available UE-specific candidate. In this example
                % we use the first available candidate to map the coded DCI bits.
                if(~isempty(candidates))
                    pdcchBits ( candidates(1,1) : candidates(1,2) ) = codedDciBits;
                end
                % PDCCH complex-valued modulated symbol generation.
                pdcchSymbols = ltePDCCH(enbConfig, pdcchBits);
            else
                % If no data to map to PDSCH no symbols and no PRBs
                pdschSymbols = [];
                prbs = zeros(0,1);
                % No PDCCH symbols
                pdcchSymbols = [];
            end
            
            % PDCCH indices generation for symbol mapping on resource grid.
            pdcchIndices = ltePDCCHIndices(enbConfig,{'1based'});

            % PDCCH symbols mapping on resource grid.
            if ~isempty(pdcchSymbols)
                subframe(pdcchIndices) = pdcchSymbols * 10^((enbConfig.PDSCH.PDCCHPower)/20);
            end
            
            % Adding SIB1            
            SFN = enbConfig.NFrame;
            % If SIB substructure is present and for the right subframes
            if (isfield(enbConfig,'SIB')) && (mod(subframeIdx,10) == 5) && (mod(SFN,2)==0)
                
                % If Enable field not specified default to On
                if ~isfield(enbConfig.SIB,'Enable')
                    enbConfig.SIB.Enable = 'On';
                    lte.internal.defaultValueWarning('SIB.Enable','On');
                end
                
                if any(strcmpi(enbConfig.SIB.Enable,{'On','Enable'}))                    
                    [subframe, SIBPRBSet] = addSIB1(enbConfig,subframe);
                    % Avoid OCNG from overwriting SIB PDSCH by adding it to
                    % prbs (these RBs are avoided when adding OCNG)
                    prbs = vertcat(prbs, SIBPRBSet); %#ok<AGROW>
                    prbs = sort(prbs); % sort in ascending order
                end
            end
            
            % Fill unused PRB with OCNG
            if strcmpi(enbConfig.OCNGPDSCHEnable,'On')
                refPRBS = (0:enbConfig.NDLRB-1).';
                if(size(prbs,2) == 1)
                    ocngPRBS = setdiff(refPRBS,prbs);
                else
                    ocngPRBS1 = setdiff(refPRBS,prbs(:,1));
                    ocngPRBS2 = setdiff(refPRBS,prbs(:,2));
                    ocngPRBS = [ocngPRBS1 ocngPRBS2];
                end
                [ocngSym,ocngInd] = getPDSCHOCNG(enbConfig,PDSCH,ocngPRBS,localSeed);
                subframe(ocngInd) = ocngSym;
            end
            
            % Fill unused PDCCH REGs (if present) with OCNG
            if strcmpi(enbConfig.OCNGPDCCHEnable,'On') 
                
                % Get the used PDCCH Indices
                usedPdcchIndices = find(pdcchSymbols~=0);

                % If there are unused PDCCH REGs, fill with OCNG
                if ~isequal(usedPdcchIndices, pdcchIndices)
                    % For the current subframe, set the PDCCH OCNG bits according 
                    % to the absolute subframe (frame & subframe) number. 
                    % To do this we reinitialize the random number generator using
                    % the absolute subframe number, localSeed and PDCCH channel ID
                    pdcchInfo = ltePDCCHInfo(enbConfig);
                    initializeRandomGenerator(enbConfig,localSeed,1);
                    pdcchOcngBits = randi([0 1],1,pdcchInfo.MTot);

                    % PDCCH complex-valued modulated symbol generation.
                    pdcchIncOcngSymbols = ltePDCCH(enbConfig, pdcchOcngBits) * 10^((enbConfig.OCNGPDCCHPower)/20);
                    
                    % Get the used PDCCH symbols
                    usedPdcchSymbols = pdcchSymbols(usedPdcchIndices);

                    % Now write the PDCCH symbols 
                    pdcchIncOcngSymbols(usedPdcchIndices) = usedPdcchSymbols;

                    % PDCCH symbols mapping on resource grid.
                    subframe(pdcchIndices) = pdcchIncOcngSymbols;
                end
            end
            
            % Generating coded CFI bits using CFI value. Given CFI value is used to 
            % indicate number of OFDM symbol used in PDCCH transmission.
            cfiBits = lteCFI(enbConfig);
            
            % PCFICH complex-valued modulated symbol generation.
            pcfichSymbols = ltePCFICH(enbConfig,cfiBits)* 10^((enbConfig.PCFICHPower)/20);
            
            % PCFICH indices generation for mapping on resource grid.
            pcfichIndices = ltePCFICHIndices(enbConfig);
            
            % Mapping PCFICH symbols on physical downlink resource elements.
            subframe(pcfichIndices) = pcfichSymbols;

            % PHICH symbol generation and mapping on REs of a resource grid.
            phichSymbols = ltePHICH(enbConfig,enbConfig.HISet) * 10^((enbConfig.PHICHPower)/20);
            phichIndices = ltePHICHIndices(enbConfig,{'1based'});
            subframe(phichIndices) = phichSymbols;
            
            % PSCH symbol generation and mapping on REs of a resource grid using
            % linear indices. These indices are made '1based' for direct mapping on
            % resource grid.
            pssSymbols = ltePSS(enbConfig);
            pssIndices = ltePSSIndices(enbConfig,0,{'1based'});
            subframe(pssIndices) = pssSymbols;
            
            % SSCH symbol generation and mapping on REs of a resource grid.
            sssSymbols = lteSSS(enbConfig);
            sssIndices = lteSSSIndices(enbConfig,0,{'1based'});
            subframe(sssIndices) = sssSymbols;
            
            % PBCH symbol generation and mapping on REs of a resource grid.
            % Note: PBCH symbols are only generated for first Subframe and
            % will be empty set for rest of the subframes in a 10ms frame.
            if(enbConfig.NSubframe==0)
                 % Generate symbols if its the first simulated frame or
                 % when mod(NFrame,4) is 0;
                 if (subframeIdx<(NSubframe+10)) || (mod(enbConfig.NFrame,4) == 0)
                    bchBits = lteMIB(enbConfig);
                    pbchBits = lteBCH(enbConfig,bchBits);
                    pbchSymbols = ltePBCH(enbConfig,pbchBits);
                end
                pbchIndices = ltePBCHIndices(enbConfig,{'1based'});
                noOfPbchSymbols = size(pbchIndices,1);
                startIdx = mod(enbConfig.NFrame,4)*noOfPbchSymbols + 1;
                endIdx = (mod(enbConfig.NFrame,4) + 1)*noOfPbchSymbols;
                subframe(pbchIndices) = pbchSymbols(startIdx:endIdx,:);
            end
            
            % Generate and map CSI-RS into grid
            if (~isempty(pdschSymbols) && any(strcmpi(PDSCH.TxScheme,'Port7-14')))
                csiRsIndices=lteCSIRSIndices(enbConfig);
                csiRsSymbols=lteCSIRS(enbConfig);
                % Remove the zero power CSI-RS symbols and indices as the
                % Zero power CSIRS symbols are always 4 ports, we only need
                % to transmit CSIRefP ports, but ensure that all RE
                % locations corresponding to zero power CSIRS symbols in
                % all ports (4) are zeroed.
                csiRsIndices(csiRsSymbols==0)=[];
                csiRsSymbols(csiRsSymbols==0)=[];
                % Extract the RE indices corresponding to CSIRefP for all
                % used and unused CSI-RS locations
                csiRsIndicesall=lteCSIRSIndices(enbConfig,'rs+unused');
                [~,reind] = lteExtractResources(csiRsIndicesall,subframe);
                % Add zeros for all CSIRS RE positions in the active
                % CSIRefP ports
                subframe(reind) = 0;
                % Map the CSI-RS symbols onto the resource elements
                subframe(csiRsIndices)=csiRsSymbols;
            end
            
            % Ensure that all subcarrier+symbol positions in all planes for
            % CellRS symbols are set to zero so that if any conflict
            % occurs, only CellRS gets transmitted
            cellRSindices = lteCellRSIndices(enbConfig);
            [~,reind] = lteExtractResources(cellRSindices,subframe);
            subframe(reind) = 0;
            
            % Cell-specific Reference Signal (CellRS) symbol generation.
            % These symbols are then mapped onto Resource Elements (RE's)
            % of a resource grid with the help of linear indices.
            cellRSsymbols = lteCellRS(enbConfig);
            subframe(cellRSindices) = cellRSsymbols;
        else
            PDSCH.CodedTrBlkSizes(1:Ncodewords,cBlkSizeIdx) = deal(0);
        end
        
        % Concatenating subframes to form a complete frame
        frame = cat(2,frame,subframe);
        if(~serialCat)
            reGrid(:,:,:,mod(subframeIdx-NSubframe,NSubframe+totalSubframes)+1) = subframe; %#ok<AGROW>
        end
    end
    
    % Update the CodedTrBlkSizes of the output config structure with the
    % one used in waveform generation
    if(length(PDSCH.CodedTrBlkSizes)>10)
        rmc.PDSCH.CodedTrBlkSizes = PDSCH.CodedTrBlkSizes(:,1:10);
    else
        rmc.PDSCH.CodedTrBlkSizes = PDSCH.CodedTrBlkSizes;
    end
    
    if(~isfield(enbConfig,'Windowing'))
        enbConfig.Windowing = 0;
    end
    
    % Time Domain mapping by performing OFDM modulation for downlink
    % symbols.
    %由时频资源frame生成时域信号timeDomainSig（FFT关系）
    [timeDomainSig,ofdmInfo,freqGrid] = lteOFDMModulate_fullfreq(enbConfig,frame); 
    if(~serialCat)
        frame = reGrid;
    end
    rmc.SamplingRate = ofdmInfo.SamplingRate;
    rmc.Nfft = ofdmInfo.Nfft;
    rmc.Windowing = ofdmInfo.Windowing;
    
    % Output the HARQ pattern used in waveform generation
    rmc.PDSCH.HARQProcessSequence = harqTable;
end
end

% Function to calculate the coded transport block sizes for the cases of 2
% CW transmission with one of the CW empty
function codedTrBlkSize = getcodedTrBlkSize(PDSCH,Gd,activeCW)
    bitspersym = find(strcmpi(PDSCH.Modulation{activeCW},{'QPSK','16QAM','64QAM','256QAM','1024QAM'}))*2;
    codedTrBlkSize(activeCW) = Gd*bitspersym*PDSCH.NLayers;
end

function [dciMsg,dciMsgBits] = getDCI(enbConfig,harqNo,NDI,ncw,prbs)

    % Create baseline DCI message structure
    istr = struct('DCIFormat',enbConfig.PDSCH.DCIFormat);
    if any(strcmpi(istr.DCIFormat,{'Format1','Format2','Format2A','Format2B','Format2C','Format2D'})) 
        % If only one PRB, set the Allocation Type to 1 as one bit represents
        % one RB in the bitmap, unless NDLRB <=10 where Allocation type must be
        % 0 as per TS 36.212 Sections 5.3.3.1.2, 5.3.3.1.5, 5.3.3.1.5A,
        % 5.3.3.1.5B, 5.3.3.1.5C and 5.3.3.1.5D
        if (size(prbs,1) == 1)
            if enbConfig.NDLRB <= 10
                istr.AllocationType = 0;
            else
                istr.AllocationType = 1;
            end
        end
        istr = lteDCI(enbConfig,istr);
        % Set the Bitmap according to the PRBSet. The Bitmap is initialized
        % with all zeros, so set all bits corresponding to active RBG to
        % '1' as per the specified PRBSet
        NRBG = numel(istr.Allocation.Bitmap); % Total RBGs for the NDLRB
        istr.Allocation.Bitmap(unique(ceil(double(prbs+1)/ceil(enbConfig.NDLRB/NRBG)))) = '1';
    end
    % Creating DCI message for given configuration
    dci = lteDCI(enbConfig,istr);

    mcsIdx = getMCSIndex(enbConfig,ncw);
    harqNo = harqNo - 1; % making HARQNo 0-based
    preCodingFieldVal = 0;

    if(strcmpi(dci.DCIFormat,'Format1'))
        % DCI format 1 is used for the scheduling of one PDSCH codeword in one cell
        dci.ModCoding = mcsIdx(1);
        dci.HARQNo = harqNo;
        dci.NewData = NDI(1);
        dci.RV = enbConfig.PDSCH.RV(1);
        dci.TPCPUCCH = 0;
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
        end
  
    elseif(strcmpi(dci.DCIFormat,'Format1A')) 
        % Compact scheduling of one PDSCH codeword in one cell, valid for
        % Port0, Port7 and TxDiversity
        dci = getFormat1A1B1DCommonSettings(enbConfig,dci,mcsIdx,harqNo,NDI);
        
    elseif(strcmpi(dci.DCIFormat,'Format1B'))
        % DCI format 1B is used for the compact scheduling of one PDSCH
        % codeword in one cell with precoding information.
        dci = getFormat1A1B1DCommonSettings(enbConfig,dci,mcsIdx,harqNo,NDI);
        dci.TPMI = enbConfig.PDSCH.PMISet;
        dci.PMI = 0; % Precoding according to TPMI

    elseif(strcmpi(dci.DCIFormat,'Format1D'))
        % DCI format 1D is used for the compact scheduling of one PDSCH
        % codeword in one cell with precoding and power offset information.
        dci = getFormat1A1B1DCommonSettings(enbConfig,dci,mcsIdx,harqNo,NDI);
        dci.TPMI = enbConfig.PDSCH.PMISet;
        dci.DlPowerOffset = 0; 

    elseif(strcmpi(dci.DCIFormat,'Format2'))
        if((enbConfig.CellRefP == 2))
            if(ncw == 2)
                if(enbConfig.PDSCH.NLayers == 2)
                    preCodingFieldVal = 2;
                end
            elseif(ncw == 1)
                if(enbConfig.PDSCH.NLayers == 1)
                    preCodingFieldVal = 5;
                end
            end
        elseif(enbConfig.CellRefP == 4)
            if(ncw == 2)
                switch(enbConfig.PDSCH.NLayers)
                    case 2
                        preCodingFieldVal = 16;
                    case 3
                        preCodingFieldVal = 33;
                    case 4
                        preCodingFieldVal = 50;
                end
            elseif(ncw == 1)
                switch(enbConfig.PDSCH.NLayers)
                    case 1
                        preCodingFieldVal = 17;
                    case 2
                        preCodingFieldVal = 34;
                end
            end        
        end

        dci.SwapFlag = 0;
        dci.HARQNo = harqNo;
        dci.ModCoding1 = mcsIdx(1);
        dci.NewData1 = NDI(1);
        dci.RV1 = enbConfig.PDSCH.RV(1);
        if(ncw==2)
            dci.ModCoding2 = mcsIdx(2);
            dci.NewData2 = NDI(2);
            dci.RV2 = enbConfig.PDSCH.RV(2);
        else
            % According to TS 36.213 Section 7.1.7.2, for DCI formats 2, 2A,
            % 2B, 2C and 2D a transport block is disabled if IMCS = 0 and if
            % rvidx = 1. Otherwise the transport block is enabled.
            dci.ModCoding2 = 0;
            dci.RV2 = 1;
        end
        dci.PrecodingInfo = preCodingFieldVal;
        dci.TPCPUCCH = 0;
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
        end
 
    elseif(strcmpi(dci.DCIFormat,'Format2A'))
        if(enbConfig.CellRefP == 4 && ncw == 2)
            switch(enbConfig.PDSCH.NLayers)
                case 3
                    preCodingFieldVal = 1;
                case 4
                    preCodingFieldVal = 2;
            end
        end

        dci.SwapFlag = 0;
        dci.HARQNo = harqNo;
        dci.ModCoding1 = mcsIdx(1);
        dci.NewData1 = NDI(1);
        dci.RV1 = enbConfig.PDSCH.RV(1);
        if(ncw==2)
            dci.ModCoding2 = mcsIdx(2);
            dci.NewData2 = NDI(2);
            dci.RV2 = enbConfig.PDSCH.RV(2);
        else
            % According to TS 36.213 Section 7.1.7.2, for DCI formats 2, 2A,
            % 2B, 2C and 2D a transport block is disabled if IMCS = 0 and if
            % rvidx = 1. Otherwise the transport block is enabled.
            dci.ModCoding2 = 0;
            dci.RV2 = 1;
        end
        dci.TPCPUCCH = 0;
        if(enbConfig.CellRefP == 4)
            dci.PrecodingInfo = preCodingFieldVal;
        end
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
        end

    elseif(strcmpi(dci.DCIFormat,'Format2B'))
        dci.TPCPUCCH = 0;
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
            dci.SRSRequest = 0;
        end
        dci.HARQNo = harqNo;
        dci.ScramblingId = enbConfig.PDSCH.NSCID;
        dci.ModCoding1 = mcsIdx(1);
        dci.NewData1 = NDI(1);
        dci.RV1 = enbConfig.PDSCH.RV(1);
        if(ncw==2)
            dci.ModCoding2 = mcsIdx(2);
            dci.NewData2 = NDI(2);
            dci.RV2 = enbConfig.PDSCH.RV(2);
        else
            % According to TS 36.213 Section 7.1.7.2, for DCI formats 2, 2A,
            % 2B, 2C and 2D a transport block is disabled if IMCS = 0 and if
            % rvidx = 1. Otherwise the transport block is enabled.
            dci.ModCoding2 = 0;
            dci.RV2 = 1;
        end

    elseif(strcmpi(dci.DCIFormat,'Format2C'))
        dci.AllocationType = 0;
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
        end
        dci.TPCPUCCH = 0;
        dci.HARQNo = harqNo;
        dci.TxIndication = getTxIndication(enbConfig,ncw);
        dci.ModCoding1 = mcsIdx(1);
        dci.NewData1 = NDI(1);
        dci.RV1 = enbConfig.PDSCH.RV(1);
        if(ncw==2)
            dci.ModCoding2 = mcsIdx(2);
            dci.NewData2 = NDI(2);
            dci.RV2 = enbConfig.PDSCH.RV(2);
        else
            % According to TS 36.213 Section 7.1.7.2, for DCI formats 2, 2A,
            % 2B, 2C and 2D a transport block is disabled if IMCS = 0 and if
            % rvidx = 1. Otherwise the transport block is enabled.
            dci.ModCoding2 = 0;
            dci.RV2 = 1;
        end

    elseif(strcmpi(dci.DCIFormat,'Format2D'))
        % DCI format used for the scheduling of up to 8 layer transmission
        % (antenna port 7 to 14) using TM10. Only applicable to the Port7, Port8
        % and Port7-14 schemes
        dci.TPCPUCCH = 0;
        if(strcmpi(enbConfig.DuplexMode,'TDD'))
            dci.TDDIndex = 1;
            dci.SRSRequest = 0;
        end
        dci.HARQNo = harqNo;
        dci.TxIndication = getTxIndication(enbConfig,ncw);
        dci.ModCoding1 = mcsIdx(1);
        dci.NewData1 = NDI(1);
        dci.RV1 = enbConfig.PDSCH.RV(1);
        if(ncw==2)
            dci.ModCoding2 = mcsIdx(2);
            dci.NewData2 = NDI(2);
            dci.RV2 = enbConfig.PDSCH.RV(2);
        else
            % According to TS 36.213 Section 7.1.7.2, for DCI formats 2, 2A,
            % 2B, 2C and 2D a transport block is disabled if IMCS = 0 and if
            % rvidx = 1. Otherwise the transport block is enabled.
            dci.ModCoding2 = 0;
            dci.RV2 = 1;
        end
        dci.REMappingAndQCL = 0;
    end
    [dciMsg,dciMsgBits] = lteDCI(enbConfig,dci);
end

function dci = getFormat1A1B1DCommonSettings(enbConfig,dci,mcsIdx,harqNo,NDI)
    dci.AllocationType = 0; % Localized allocation for Type 2
    % Get the first RB
    prbSorted = sort(enbConfig.PDSCH.PRBSet,1);
    dci.VRBStart = prbSorted(1,1); 
    dci.VRBLength = size(enbConfig.PDSCH.PRBSet,1); % 
    % Calculate the RIV
    dci.Allocation.RIV = calculateRIV(enbConfig.NDLRB, dci,false);
    dci.ModCoding = mcsIdx(1);
    dci.HARQNo = harqNo;
    dci.NewData = NDI(1);
    dci.RV = enbConfig.PDSCH.RV(1);
    dci.TPCPUCCH = 0;
    if(strcmpi(enbConfig.DuplexMode,'TDD'))
        dci.TDDIndex = 1;
    end
end

function txIndication = getTxIndication(enbConfig,ncw)   
    % Mapping TxIndication as per Table 5.3.3.1.5C-1 in TS 36.212
    switch double(enbConfig.PDSCH.NLayers)
        case 1
            txIndication = enbConfig.PDSCH.NSCID;
            if (strcmpi(enbConfig.PDSCH.TxScheme,'Port8'))
                % TxIndication = 2 if NCSID = 0 and TxIndication = 3 if
                % NCSID = 1
                txIndication = double(enbConfig.PDSCH.NSCID) + 2;
            end
        case 2
            if(ncw==1)
                txIndication = 4;
            elseif(ncw==2)
                txIndication = enbConfig.PDSCH.NSCID;
            end
        case 3
            if(ncw==1)
                txIndication = 5;
            elseif(ncw==2)
                txIndication = 2;
            end
        case 4
            if(ncw==1)
                txIndication = 6;
            elseif(ncw==2)
                txIndication = 3;
            end
        case {5, 6, 7, 8}
            txIndication = double(enbConfig.PDSCH.NLayers) - 1;
    end
end

% Select the 'best fit' MCS indices from PDSCH MCS table 1 or MCS table 2.
% For multiple codewords the indices must be selected from the same table,
% with a preference for table 1 if matches are found in both.
% Select the retransmission indices for the given modulation order if no 
% transport block length match is found.
function mcsIdx = getMCSIndex(enbConfig,ncw)
    
    % Calculate the NPRB to be used in the subsequent calculations
    ndashprb = size(enbConfig.PDSCH.PRBSet,1);
    nprb = ndashprb;
    dupinfo = lteDuplexingInfo(enbConfig);
    if (ndashprb>0) && strcmpi(dupinfo.SubframeType,'Special') && strcmpi(enbConfig.DuplexMode,'TDD')
        % Apply rules in 36.213 Section 7.1.7 for selecting the TB Size
        if ((enbConfig.SSC == 9) && strcmpi(enbConfig.CyclicPrefix,'Normal')) || ((enbConfig.SSC == 7) && strcmpi(enbConfig.CyclicPrefix,'Extended'))
            nprb = max(floor(ndashprb*0.375),1);
        else
            nprb = max(floor(ndashprb*0.75),1);
        end
    end

    % Get transport block size of current subframe being processed
    tbSize = zeros(1,ncw);  
    for n = 1:ncw
        tbSize(n) = enbConfig.PDSCH.TrBlkSizes(n,mod(enbConfig.NSubframe,size(enbConfig.PDSCH.TrBlkSizes,2))+1); 
    end
    
    % Calculate the number of spatial layers per codeword for use in the TBS lookup
    slayers = ones(1,ncw);
    % If using a transmission scheme capable of supporting multiple spatial layers 
    % then divide the total number of spatial layers across the codewords
    if any(strcmpi(enbConfig.PDSCH.TxScheme,{'CDD','SpatialMux','Port7-8','Port7-14'}))
        slayers = enbConfig.PDSCH.NLayers*slayers;
        % Divide the total number of layers between the codewords, rounding 
        % up on the second codeword     
        if ncw==2
            slayers(1) = fix(slayers(1)/2);
            slayers(2) = slayers(2)-slayers(1);
        end
    end
    
    % Identify the most appropriate PDSCH MCS table for the I_MCS
    mcsIdx = zeros(1,ncw);
    candidates = zeros(1,ncw);
    tables = {'PDSCH','PDSCHTable2','PDSCHTable3'};   % I_MCS table types
    maxtms = 0;
    % Walk through the pair of tables
    for t = 1:length(tables)
        [itbs,modulation] = lteMCS(tables{t});
        itbs(isnan(itbs)) = -1;         % Reassign the reserved I_TBS for use with lteTBS
        % Find possible I_MCS matches for each codeword
        tms = 0;             % Weighting measure of the I_MCS match
        for n = 1:ncw
            tbSizes = lteTBS(nprb,itbs,slayers(n));     % Get possible TBS sizes (for I_TBS = 0...33) that can be signaled for this allocation bandwidth  and layers        
            % Identify matches by modulation scheme and either target TBS or TBS equal to 0 (reserved) 
            matches = find(and(strcmpi(enbConfig.PDSCH.Modulation{n}, modulation),or(tbSize(n) == tbSizes, tbSizes==0)));      
            if isempty(matches)
                tms = 0;     % Penalize the case that there is no match at all
                break;       % then escape from the inner loop across codewords early
            else
                tms = tms + numel(matches);     % Increase the 'best fit' measure by the number of matches (1 or 2)
                candidates(n) = matches(1)-1;   % Cache the first match (there can only be 1 or 2, where one of then is the retransmission index)
            end
        end
        % Test whether these I_MCS candidates are better than the alternative
        if tms > maxtms
            maxtms = tms;
            mcsIdx = candidates;
        end
    end
end    

% Predicate test of whether the transmission scheme supports UE-specific DM-RS
function uers = isUeRS(txScheme)
    uers=any(strcmpi(txScheme,{'Port5' 'Port7-8' 'Port8' 'Port7-14'}));
end

% OCNG symbols as defined in TS 36.101 Annex A.5
function [sym,indices] = getPDSCHOCNG(enb,ocngPdsch,prbset,localSeed)
    
    sym = [];    
    indices = [];
    if(~isempty(prbset))

        if strcmpi(enb.OCNGPDSCH.TxScheme,'Port7-14')
           % Use the default values for CSIRSPeriod and
           % ZeroPowerCSIRSPeriod if not present at input
           if ~isfield(enb,'CSIRSPeriod')
               enb.CSIRSPeriod = 'Off';
           end
           if ~isfield(enb,'ZeroPowerCSIRSPeriod')
               enb.ZeroPowerCSIRSPeriod = 'Off';
           end
        end
        
        ocngPdsch.TxScheme = enb.OCNGPDSCH.TxScheme;
        ocngPdsch.Modulation = enb.OCNGPDSCH.Modulation;
        ocngPdsch.RNTI = enb.OCNGPDSCH.RNTI;
        
        [indices,info] = ltePDSCHIndices(enb,ocngPdsch,prbset);

        capacity = info.G;
        
        % For the current subframe, set the PDSCH OCNG bits according 
        % to the absolute subframe (frame & subframe) number. 
        % To do this we reinitialize the random number generator using
        % the absolute subframe number, localSeed and PDSCH channel ID
        initializeRandomGenerator(enb,localSeed,0);
        ocngData = randi([0 1],1,capacity);

        powerAdjPerRE = 10^((enb.OCNGPDSCHPower)/20); 

        sym = ltePDSCH(enb,ocngPdsch,ocngData)*powerAdjPerRE;
    end    
end

function completeStruct = checkFields(enbConfig)
    completeStruct = true;
    fields1 = isfield(enbConfig,{'NDLRB','CellRefP','NCellID','CyclicPrefix',...
        'CFI','PCFICHPower','Ng','PHICHDuration','HISet','PHICHPower','NFrame','NSubframe','TotSubframes',...
        'DuplexMode','PDSCH','OCNGPDCCHEnable','OCNGPDCCHPower','OCNGPDSCHEnable','OCNGPDSCHPower','OCNGPDSCH',});
    if (isfield(enbConfig,'PDSCH'))
        fields2 = isfield(enbConfig.PDSCH,{'TxScheme','Modulation','NLayers','Rho','RNTI',...
            'RVSeq','RV','NHARQProcesses','NTurboDecIts','PRBSet','TrBlkSizes','DCIFormat','PDCCHFormat','PDCCHPower'});
    end
    if (isfield(enbConfig,'OCNGPDSCH'))
        fields3 = isfield(enbConfig.OCNGPDSCH,{'TxScheme','Modulation','RNTI'});
    end
    if isfield(enbConfig,'OCNG')
        % If the old interface is provided, it should overwrite all fields
        % of the new interface with the values required for backward
        % compatibility 
        fields3 = false;
    end
    if(~all(fields1) || ~all(fields2) || ~all(fields3))
        completeStruct = false;
    end
end

% Calculate the HARQ table according to the procedure detailed in
% R5-095777 ("Scheduling of retransmissions and number of active HARQ
% processes for DL performance RMC-s")
function harqProcessSequence  = getHARQTable(rmc)

    noHarqProcesses = rmc.PDSCH.NHARQProcesses;
        
    % Identify all active subframes (from downlink and special subframes and non-zero TrBlkSizes)
    info = arrayfun(@(x)lteDuplexingInfo(setfield(rmc,'NSubframe',x)),0:9); 
    dlandspecialsfs = arrayfun(@(x)any(strcmpi(x.SubframeType,{'Downlink','Special'})),info);
    
    % Expand the TrBlkSizes to its period if length is not full frames
    trBlkSizes = trBlkExpand(rmc.PDSCH, dlandspecialsfs);
    trBlkSizes = trBlkSizes(1,:); % Need only 1 CW for HARQ calculation
 
    % Initialize harqProcessSequence of all zeros for a frame
    harqProcessSequence = zeros(1,10);
    
    % If there are no non-zero subframes, return
    if ~any(trBlkSizes)
        return
    end
    
    % Special case of 1 HARQ process, only the first TrBlk is transmitted
    if noHarqProcesses == 1
        harqProcessSequence(1) = 1;
        return
    else
        % If the number of unique non-zero TrBlkSizes in downlink & special
        % subframes exceeds number of HARQ processes this can't be
        % transmitted
        uniqueSizes = numel(unique(trBlkSizes)) - (~all(trBlkSizes));
        if uniqueSizes > noHarqProcesses 
             error('lte:error', 'The number of unique values in TrBlkSizes (%d) cannot exceed NHARQProcesses (%d)',uniqueSizes,noHarqProcesses);
        end

        % Initialize frame counter and loop variable
        nframeSet = 0;
        tableComplete = false;
        % Initialize HARQ process sequence and index
        harqProcessQueue = 1:noHarqProcesses;
        harqTable = zeros(numel(trBlkSizes),1);
        while (~tableComplete)
            % Assign HARQ processes according to the TrBlkSize for a
            % frameSet (frameSet contains 1 or more frames)
            for sf = 1:numel(trBlkSizes)
                [harqTable(sf,nframeSet+1), harqProcessQueue] = getHARQId(harqTable,sf,harqProcessQueue,trBlkSizes);
            end
            nframeSet = nframeSet+1; % Increment the frameSet counter
            
            % When the frame starts repeating, create a properly formatted
            % sequence and exit the loop
            if nframeSet> 1  && isequal(harqTable(:,1),harqTable(:,nframeSet))
                harqTable(:,nframeSet) = []; % Remove the repeating frameSet
                harqProcessSequence = harqTable(:)';
                tableComplete = true; % Set the flag to exit the while loop
            end
        end
    end
end

% Function to calculate the next HARQ ID in the sequence. 
% If the trBlkSize for the current location(subframe) is same as that for a
% previous assignment of the HARQId or there was no previous assignment of
% the HARQId, assign the HARQId, otherwise move on to the next HARQId till
% the condition is met
function [harqId,harqProcessQueue] = getHARQId(harqTable,sfidx,harqProcessQueue,trBlkSizes)
    
    if trBlkSizes(sfidx) == 0
        harqId = 0; % No HARQ process for this subframe
        return
    else
        harqId = [];
        tempharqIdx = 1;
    end

    while isempty(harqId)
        % Get the next HARQId from the queue to check if that is the HARQId
        % to be used for the subframe
        tempharqId = harqProcessQueue(tempharqIdx);
        % Check if there was a previous assignment of this HARQId
        prevTrBlkIdx = find(harqTable==tempharqId);
        if isempty(prevTrBlkIdx)
            % No previous assignment, so last TrBlkSize assigned is current
            prevTrBlkWithThisHARQId = trBlkSizes(sfidx);
        else
            prevTrBlkWithThisHARQId = trBlkSizes(mod(prevTrBlkIdx(1)-1, numel(trBlkSizes)) +1);
        end
        % Now check if the position corresponding to the Idx has
        % the same transport block size (sfidx is modulo valid TrBlkSize)
        if (trBlkSizes(sfidx) == prevTrBlkWithThisHARQId)
            harqId = tempharqId;
            % Now move this ID to the back of the queue so that other IDs
            % with the same TrBlkSize will be scheduled next
            harqProcessQueue(harqProcessQueue==harqId) = [];
            harqProcessQueue(end+1) = harqId; %#ok<AGROW>
        else
            % No need for modulo as there will be one valid HARQId
            tempharqIdx = tempharqIdx+1;
        end
        
    end
end


% Function to calculate TrBlkSizes including any cyclic repetition.
% If the TrBlkSizes vector is not an integer multiple of frame size, we
% need to expand TrBlkSizes up to its maximum possible period (lcm of the
% number of values in TrBlkSizes and subframes in a frame). The actual
% period might be less if the contents of TrBlkSizes have same values.
function trBlkSize = trBlkExpand(PDSCH, dlandspecialsfs)
    % Initialize the number of frames and block sizes
    nframes = lcm(numel(PDSCH.TrBlkSizes(1,:)),10)/10;
    trBlkSize = zeros(size(PDSCH.TrBlkSizes,1),nframes*10);
    for sf = 0:nframes*10-1
        if dlandspecialsfs(mod(sf,10)+1)
            trBlkSize(:,sf+1) =  PDSCH.TrBlkSizes(:,mod(sf,size(PDSCH.TrBlkSizes,2)) + 1);
        end
        % If current frame is same as first frame, adjust nframes and exit early 
        if ~mod(sf+1,10) && (sf>10)
           if isequal(trBlkSize(:,1:10),trBlkSize(:,sf-9+1:sf+1)) 
               nframes = floor(sf/10); % correct the nframes
               trBlkSize = trBlkSize(:,1:nframes*10); % truncate the TrBlkSize
               break;
           end
        end    
    end
end

    % Add SIB1 to a subframe
function [ subframe, PRBSet ] = addSIB1( enbConfig, subframe )
% This function adds the SIB1 to the provided subframe.
% Parameters:
%   - subframe: subframe to add the SIB1 to
%   - enbConfig: enb configuration structure. This structure must contain a
%   SIB substructure with the SIB1 parameters. These are:
%       - Data: vector of bits to map to the SIB1
%       - DCIFormat: DCI format used, can be 'Format1A' and 'Format1C'
%       - AllocationType: indicates whether the resources allocated are
%       localized (0) or distributed (1)
%       - VRBStart: virtual RB allocation starting resource block
%       - VRBLength: length of virtually contiguously allocated resource block
%       Only required for Format1A
%           - N1APRB: Number of physical resource blocks for Format1A
%       Only required for distributed allocation AllocationType = 1
%           - Gap: Distributed allocation gap, can take the values 0 (gap1
%           and 1 (gap2)

if ~isfield(enbConfig.SIB,'Data')
    error('lte:error','The function call resulted in an error: Could not find a SIB structure field called Data.');
end

if ~isfield(enbConfig.SIB,'DCIFormat')
    enbConfig.SIB.DCIFormat = 'Format1A';
    lte.internal.defaultValueWarning('SIB.DCIFormat','Format1A');    
end

if ~any(strcmpi(enbConfig.SIB.DCIFormat,{'Format1A' 'Format1C'}) )
    error('lte:error','The function call resulted in an error: DCIFormat in SIB structure should be one of {''Format1A'', ''Format1C''}.');
end

if ~isfield(enbConfig.SIB,'AllocationType')
    enbConfig.SIB.AllocationType = 0;
    lte.internal.defaultValueWarning('SIB.AllocationType','0');    
end

if ~any(enbConfig.SIB.AllocationType == [0 1])
    error('lte:error','The function call resulted in an error: The SIB field AllocationType should take one of the following values {0, 1}');
end

if enbConfig.SIB.AllocationType == 1 % distributed
    if ~isfield(enbConfig.SIB,'Gap')
        error('lte:error','The function call resulted in an error: Could not find a SIB structure field called Gap.');
    end
end

if ~isfield(enbConfig.SIB,'VRBStart')
    error('lte:error','The function call resulted in an error: Could not find a SIB structure field called VRBStart.');
end

if ~isfield(enbConfig.SIB,'VRBLength')
    error('lte:error','The function call resulted in an error: Could not find a SIB structure field called VRBLength.');
end
  
% Create and populate PDSCH structure for SIB1
PDSCHsib1.RNTI=65535; % SI-RNTI, 36.213
if (enbConfig.CellRefP==1)
    PDSCHsib1.TxScheme='Port0';
    PDSCHsib1.NLayers=1;
else
    PDSCHsib1.TxScheme='TxDiversity';
    PDSCHsib1.NLayers=enbConfig.CellRefP;
end
PDSCHsib1.Modulation='QPSK'; %as per 36.213 Section 7.1.7.1

% RV setup for SIB1 per TS 36.321 Section 5.3.1
SFN = enbConfig.NFrame;
k = mod(floor(SFN/2),4);
RVK = mod(ceil(3/2*k),4);
PDSCHsib1.RV = RVK;

% Get TBS
sib1Length = length(enbConfig.SIB.Data);

if (strcmpi(enbConfig.SIB.DCIFormat,'Format1A'))
    % Check for N1APRB only for Format1A
    if ~isfield(enbConfig.SIB,'N1APRB')
        enbConfig.SIB.N1APRB = 2;
        tbSizes = lteTBS(enbConfig.SIB.N1APRB,0:26);
        if sib1Length>max(tbSizes) 
            enbConfig.SIB.N1APRB = 3;
        end
    end
    % N1APRB: 36.212 Section 5.3.3.1.3 
    N1APRB = enbConfig.SIB.N1APRB;    
    if ~any(N1APRB == [2 3])
        error('lte:error','The function call resulted in an error: SIB1 N1APRB field should take one of the following values {2, 3}');
    end
    % Tables from 36.213 Section 7.1.7.2.1    
    tbSizes = lteTBS(N1APRB,0:26);    
else % format 1C
    % From 36.213 Sect 7.1.7.2.3
    tbSizes = [40 56 72 120 136 144 176 208 224 256 280 296 328 336 392 488 552 600 632 696 776 840 904 1000 1064 1128 1224 1288 1384 1480 1608 1736];
end

% Find nearest allowed TBS length
if sib1Length>max(tbSizes) 
    error('lte:error',['The function call resulted in an error: SIB1 data too long, maximum allowed length is ' num2str(max(tbSizes)) ' bits.']);
else
    % Find closer allowed TBS length
    iTBS = find(tbSizes>=sib1Length,1)-1; %I_TBS is zero based
    tbSize = tbSizes(iTBS+1);
    % Pad with zeros
    paddedSib1Bits = enbConfig.SIB.Data(:);
    paddedSib1Bits=[paddedSib1Bits; zeros(tbSize-sib1Length,1)];
end

% Populate dcistr structure.
dcistr.AllocationType = enbConfig.SIB.AllocationType;  
dcistr.NDLRB = enbConfig.NDLRB;
dcistr.ModCoding = iTBS; % for SIB: I_TBS = I_MCS as per 36.213 Section 7.1.7.2 and Section 7.1.7    
if (strcmpi(enbConfig.SIB.DCIFormat,'Format1A'))
    dcistr.DCIFormat='Format1A';    
    dcistr.RV = PDSCHsib1.RV; % RV setup for SIB1 per TS 36.321 Section 5.3.1
    if N1APRB == 2  % 36.212, Section 5.3.3.1.3
        dcistr.TPCPUCCH = 0;
    else
        dcistr.TPCPUCCH = 1;
    end
else    
    dcistr.DCIFormat='Format1C';
    if dcistr.AllocationType ~= 1
        error('lte:error', 'The function call resulted in an error: The SIB field AllocationType can only take the value 1 (distributed) for DCIFormat Format1C');
    end    
end

% gap field
if enbConfig.SIB.AllocationType % distributed
    if ~any(enbConfig.SIB.Gap == [0 1])
        error('lte:error', 'The function call resulted in an error: The SIB field Gap must be 0 (gap1) or 1 (gap2)');
    end
    % As per TS 36.212, Section 5.3.3.1.3 (DCI Format1A), the Gap value is
    % signaled via "New data indicator" field if the DCI message Format1A
    % is scrambled by SI-RNTI.
    if(strcmpi(dcistr.DCIFormat,'Format1A'))
        dcistr.NewData = enbConfig.SIB.Gap;
    else
        dcistr.Allocation.Gap = enbConfig.SIB.Gap;
    end
end

% RIV calculation
dcistr.Allocation.RIV = calculateRIV(enbConfig.NDLRB,enbConfig.SIB,true);

% Calculate SIB PRBSet from DCI
% As mentioned, the Gap value is signaled via "New data indicator" field
% if the DCI message Format1A is scrambled by SI-RNTI (TS 36.212 Section
% 5.3.3.1.3 (DCI Format1A)). However,the lteDCIResourceAllocation function
% always expects the Gap signaled in the "Allocation.Gap" field. Create a
% copy of dcistr to pass to lteDCIResourceAllocation(). Set the
% "AllocationType.Gap" field to carry the Gap value. Note that the actual
% mapped DCI message is correctly carrying the Gap value within the "New
% data indicator" field (structure dcistr).
tmpdcistr = dcistr;
if enbConfig.SIB.AllocationType % distributed
    if(strcmpi(dcistr.DCIFormat,'Format1A'))
        tmpdcistr.Allocation.Gap = enbConfig.SIB.Gap;
    end
end
PRBSet = lteDCIResourceAllocation(enbConfig,tmpdcistr);

% SIB1 DL-SCH/PDSCH transmission
PDSCHsib1.PRBSet = PRBSet;

% SIB1 PDSCH indices
[sib1Indices, sibPDSCHInfo] = ltePDSCHIndices(enbConfig,PDSCHsib1,PDSCHsib1.PRBSet);

% Channel code SIB1 bits
codedSib1Bits = lteDLSCH(enbConfig,PDSCHsib1,sibPDSCHInfo.G,paddedSib1Bits);

% SIB1 PDSCH symbols, map to grid
pdschSymbols = ltePDSCH(enbConfig,PDSCHsib1,codedSib1Bits); 
subframe(sib1Indices) = pdschSymbols;

% Add PDCCH
pdcchDims = ltePDCCHInfo(enbConfig);                        
pdcchBits = -1*ones(1,pdcchDims.MTot);
pdcchConfig.NDLRB = enbConfig.NDLRB;
pdcchConfig.PDCCHFormat = 2;
candidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});

pdcchIndices = ltePDCCHIndices(enbConfig);
pdcchConfig.RNTI = PDSCHsib1.RNTI;

[~,dciBits] = lteDCI(enbConfig,dcistr);
codedDciBits = lteDCIEncode(pdcchConfig,dciBits);
pdcchBits ( candidates(1,1) : candidates(1,2) ) = codedDciBits;
pdcchSymbols = ltePDCCH(enbConfig, pdcchBits);
% Assign PDCCH symbols to the subframe grid, but only where they are
% non-zero, to avoid overwriting any reference PDSCH DCI that might be
% present
subframe(pdcchIndices(abs(pdcchSymbols)~=0)) = pdcchSymbols(abs(pdcchSymbols)~=0);
end

% Calculates the RIV based on TS 36.213 Section 7.1.6.3
function RIV = calculateRIV(nRB, dciConfig, isSIB) 
%   Calculates the RIV based on TS 36.213 Section 7.1.6.3
%   Parameters:
%       - nRB: is the bandwidth in resource blocks
%       - dciConfig: SIB/DCI structure, the following fields are required:
%           - AllocationType (0: localized; 1 distributed)
%           - DCIFormat (Format1A, format1C)
%           - VRBStart
%           - VRBLength
%           - Gap (0: gap1; 1: gap2)
%       - isSIB: Flag to indicate SIB transmission where a SIB specific
%                error is thrown for incorrect configuration

    VRBStart = dciConfig.VRBStart;
    VRBLength = dciConfig.VRBLength;
    dciFormat = dciConfig.DCIFormat;
    allocationType = dciConfig.AllocationType;
    
    RIV = [];
    
    % Determine Nvrb
    if(allocationType == 1) % distributed
        % According to 36.211 Section 6.2.3.2
        gap = dciConfig.Gap;
        Nvrb = distributedNvrb(nRB,gap);          
    elseif(allocationType == 0) % localized
        % According to 36.211 Section 6.2.3.1
        Nvrb = nRB;
    end
    
    if(any(strcmpi(dciFormat, {'Format1A','Format1B','Format1D'})))        
        % From 36.213 Section 7.1.6.3
        if (VRBStart<0) || (VRBStart>(Nvrb-1))
            error('lte:error',['The function call resulted in an error: The SIB VRBStart parameter should be in the range [0, ' num2str(Nvrb-1) '] for the provided set of parameters'])
        end
        if (VRBLength<=(Nvrb-VRBStart) && VRBLength>=1)
            if((VRBLength-1)<= floor(nRB/2))
                RIV = nRB*(VRBLength-1) + VRBStart;
            else
                RIV = nRB*(nRB-VRBLength+1) + (nRB-1-VRBStart);
            end
        else
            if isSIB
                % If SIB creation, return error
                error('lte:error',['The function call resulted in an error: The SIB VRBLength parameter should be in the range [1, ' num2str(Nvrb-VRBStart) '] for the provided set of parameters'])
            else
                % If DCI creation, return empty RIV
                return;
            end    
        end
    elseif(strcmpi(dciFormat, 'Format1C'))        
        % Set Nstep: 36.213 Table 7.1.6.3-1
        if (nRB>5 && nRB<50)
            Nstep = 2;
        elseif(nRB>49 && nRB<111)
            Nstep = 4;
        end

        % From 36.213 Section 7.1.6.3
        % VRBStart valid values for Format1C
        VRBStartUpperLimit = (floor(Nvrb/Nstep)-1)*Nstep;
        if (VRBStart<0) || ( VRBStart>VRBStartUpperLimit ) || mod(VRBStart,Nstep)~=0
            VRBStartValidSet = 0:Nstep:VRBStartUpperLimit; % valid set of values for VRBStart for Format1C
            error('lte:error',['The function call resulted in an error: The SIB VRBStart parameter should be one of {' num2str(VRBStartValidSet) '} for the provided set of parameters'])
        end

        % Lcrb valid values for Format1C
        % Lcrb upper limits:
        %    - According to 36.213 Section 7.1.6.3 Lcrb <= floor(Nvrb/Nstep)*Nstep
        %    - However VRBLength_prime <= Nvrb_prime - VRBStart_prime
        % Considering the definitions od VRBLength_prime, Nvrb_prime and
        % VRBStart_prime it follows that:
        %    VRBLength <= Nstep*(Nvrb_prime - VRBStart_prime) = floor(Nvrb/Nstep)*Nstep - VRBStart
        LcrbUpperLimit = floor(Nvrb/Nstep)*Nstep - VRBStart;
        if (VRBLength<Nstep) || (VRBLength>LcrbUpperLimit) || (mod(VRBLength,Nstep)~=0)
            VRBLengthValidSet = Nstep:Nstep:LcrbUpperLimit; % valid set of values for VRBLength for Format1C           
            error('lte:error',['The function call resulted in an error: The SIB VRBLength parameter should be one of {' num2str(VRBLengthValidSet) '} for the provided set of parameters'])
        end
        
        VRBLength_prime = VRBLength/Nstep;
        VRBStart_prime = VRBStart/Nstep;
        Nvrb_prime = floor(Nvrb/Nstep);
        
        if(VRBLength_prime<=(Nvrb_prime-VRBStart_prime))
            if((VRBLength_prime-1) <= floor(Nvrb_prime/2))
                RIV = Nvrb_prime*(VRBLength_prime - 1) + VRBStart_prime;
            else
                RIV = Nvrb_prime*(Nvrb_prime - VRBLength_prime + 1) + (Nvrb_prime - 1 - VRBStart_prime);
            end
%         else % Condition not needed as already checked with the range of VRBLength
%             error('lte:error',['The SIB Lcrb parameter should be <=' num2str(Nvrb_prime*Nstep - VRBStart) ' for the provided bandwidth'])
        end
    end
end

% distributedNvrb(nRB,gap) calculates Nvrb (DL) as specified in 36.211,
% Section 6.2.3.2
function [Nvrb] = distributedNvrb(nRB,gap)
% Calculate the DL Nvrb parameter. 
%   - nRB is the bandwidth in resource blocks
%   - gap indicates whether the 1st or second gap is used:
%       - gap = 0: gap1
%       - gap = 1: gap2
    
    if (nRB<=49) && (gap == 1) % gap2. Table 6.2.3.2-1 in 36.211
         error('lte:error', 'The function call resulted in an error: SIB field Gap = 1 (gap2) is only defined for bandwidths > 49 RBs');
    end

    if(nRB>5 && nRB<11)
        idx = 1;
    elseif(nRB==11)
        idx = 2;
    elseif(nRB>11 && nRB<20)
        idx = 3;
    elseif(nRB>19 && nRB<27)
        idx = 4;
    elseif(nRB>26 && nRB<45)
        idx = 5;
    elseif(nRB>44 && nRB<50)
        idx = 6;
    elseif(nRB>49 && nRB<64)
        idx = 7;
    elseif(nRB>63 && nRB<80)
        idx = 8;
    elseif(nRB>79 && nRB<111)
        idx = 9;
    end

    % Table 6.2.3.2-1 from 36.211
    g1 = [ceil(nRB/2) 4 8 12 18 27 27 32 48];
    g2 = [0 0 0 0 0 0 9 16 16];

    if(gap) %   gap = 1 (gap2)
        Ngap = g2(idx);
        Nvrb = floor(nRB/(2*Ngap))*(2*Ngap);
    else    %   gap = 0 (gap1)
        Ngap = g1(idx);
        Nvrb = 2*min(Ngap,(nRB-Ngap));
    end
end

% This function initializes the rng seed depending on localSeed, absolute
% subframe number and channel ID(0 for PDSCH, 1 for PDCCH)
function initializeRandomGenerator(enb,localSeed,channelID)
    rng(mod(enb.NFrame*10+enb.NSubframe+localSeed+channelID*100,2^32)); 
end

