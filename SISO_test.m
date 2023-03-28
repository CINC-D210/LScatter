clear all; clc;close all;
clear functions
disp('Simulating the LTE Mode 1: Single Tx and Rx antrenna');
%% LTE的参数设置
% PDSCH   生成资源块的信息
numTx          = 1;    % Number of transmit antennas
numRx          = 1;    % Number of receive antennas
chanBW       = 6;    % Index to chanel bandwidth used [1,....6]
contReg       = 1;    % No. of OFDM symbols dedictaed to control information [1,...,3]
modType      =  2;   % Modulation type [1, 2, 3] for ['QPSK,'16QAM','64QAM']
% DLSCH   Turbo解码的信息
cRate            = 1/3; % Rate matching target coding rate 
maxIter         = 6;     % Maximum number of turbo decoding terations  
fullDecode    = 0;    % Whether "full" or "early stopping" turbo decoding is performed
% Channel model
chanMdl        =  'frequency-selective-low-mobility'; 
% chanMdl        =  'flat-high-mobility'; 
corrLvl           = 'Low'; 
% Simulation parametrs
Eqmode        = 2;      % Type of equalizer used [1,2] for ['ZF', 'MMSE']
chEstOn        = 1;     % Whether channel estimation is done or ideal channel model used
maxNumErrs = 1e7; % Maximum number of errors found before simulation stops
maxNumBits = 1e7;  % Maximum number of bits processed before simulation stops
visualsOn     = 1;      % Whether to visualize channel response and constellations
snrdB            = 18;   % Value of SNR used in this experiment
%% Set simulation parametrs & initialize parameter structures
[prmLTEPDSCH, prmLTEDLSCH, prmMdl] = commlteSISO_initialize( chanBW, contReg,  modType, Eqmode,...
    cRate,maxIter, fullDecode, chanMdl, corrLvl, chEstOn, maxNumErrs, maxNumBits);
clear chanBW contReg numTx numRx modType Eqmode cRate maxIter fullDecode chanMdl corrLvl chEstOn maxNumErrs maxNumBits;
%至此 程序的变量仅有3个结构体变量： prmLTEPDSCH资源块生成相关  prmLTEDLSCH解码相关   prmMdl信道相关
% 2个常数： snrdB 好像没用到？？   visualsOn：可视化
%%
zReport_data_rate(prmLTEDLSCH, prmLTEPDSCH);  %打印参数信息
hPBer = comm.ErrorRate;  %hPBer统计误码率
% snrdB=prmMdl.snrdBs(end);
maxNumErrs=prmMdl.maxNumErrs;
maxNumBits=prmMdl.maxNumBits;
%% 
nS = 0; % Slot number, one of [0:2:18]
Measures = zeros(3,1); %initialize BER output

%%  
% 1 TX
% 1.1 产生用户bit Generate payload  按照ns（时隙号）产生随机比特
dataIn = genPayload(nS,  prmLTEDLSCH.TBLenVec);
% 1.2 加CRC  Transport block CRC generation  %加24位CRC校验
tbCrcOut1 =CRCgenerator(dataIn);
% 1.3 信道编码 Channel coding includes - CB segmentation, turbo coding, rate matching,
% bit selection, CB concatenation - per codeword
[data, Kplus1, C1] = lteTbChannelCoding(tbCrcOut1, nS, prmLTEDLSCH, prmLTEPDSCH);
%1.4 Scramble codeword    加绕 位数不变
scramOut = lteScramble(data, nS, 0, prmLTEPDSCH.maxG);
%1.5 Modulate  调制
modOut = Modulator(scramOut, prmLTEPDSCH.modType);
%1.6 Generate Cell-Specific Reference (CSR) signals  产生参考信号
csr = CSRgenerator(nS, prmLTEPDSCH.numTx);%CSR  调制星座 200*2*2？
%1.7 将用户数据、CSR及其他信令开销映射到时频资源格上 Resource grid filling  
E=8*prmLTEPDSCH.Nrb;  %？？？*8？？
csr_ref=reshape(csr(1:E),2*prmLTEPDSCH.Nrb,4);%重新排列CSR
txGrid = REmapper_1Tx(modOut, csr_ref, nS, prmLTEPDSCH);  %************把用户数据和CSR映射到资源块上*********
%1.8 OFDM transmitter
txSig = OFDMTx(txGrid, prmLTEPDSCH); %发送的串行OFDM
% 2 Channel
% SISO Fading channel
[rxFade, chPathG] = MIMOFadingChan(txSig, prmLTEPDSCH, prmMdl);
[idealhD,idealhD_grid]= lteIdChEst(prmLTEPDSCH,  prmMdl, chPathG, nS);  %%*********估计了理想的频域信道 方法？？*********

% Add AWG noise
nVar = 10.^(0.1.*(-snrdB));
rxSig =  AWGNChannel(rxFade, nVar);
% 3 RX
% 3.1 把接收到的信号重新排列，IFFT到频域栅格  OFDM Rx
rxGrid = OFDMRx(rxSig, prmLTEPDSCH);
% 3.2 按层提取用户数据和CSR updated for numLayers -> numTx
[dataRx, csrRx, idx_data] = REdemapper_1Tx(rxGrid, nS, prmLTEPDSCH);
% 3.3 用CSR进行信道估计  MIMO channel estimation
if prmMdl.chEstOn
    %**********用接收到的CSR，和其原本值估计信道，并插值出完整的视频时频栅格上对应的频域信道**************
    chEst = ChanEstimate_1Tx(prmLTEPDSCH, csrRx,  csr_ref, 'interpolate'); %整个栅格上信道的估计值
%     chEst = ChanEstimate_1Tx(prmLTEPDSCH, csrRx,  csr_ref, 'hybrid'); %整个栅格上信道的估计值
%     chEst = ChanEstimate_1Tx(prmLTEPDSCH, csrRx,  csr_ref, 'average'); %整个栅格上信道的估计值
    hD=chEst(idx_data); %取出用户子载波上的频域信道值
else
    hD = idealhD;  %各个子载波上、basic time内的频域信道
end
% 3.3 信道均衡（频域） Frequency-domain equalizer  用得到的信道对数据做均衡
yRec = Equalizer(dataRx, hD, nVar, prmLTEPDSCH.Eqmode);
% 3.4 解调用户数据  Demodulate
demodOut = DemodulatorSoft(yRec, prmLTEPDSCH.modType, nVar);
% 3.5 解绕 Descramble both received codewords
rxCW=lteDescramble(demodOut,nS,0,prmLTEPDSCH.maxG);
% 3.6 信道编码解码  Channel decoding includes - CB segmentation, turbo decoding, rate dematching
[decTbData1, ~,~] = lteTbChannelDecoding(nS, rxCW, Kplus1, C1,  prmLTEDLSCH, prmLTEPDSCH);
% 3.7 CRC校验 Transport block CRC detection
[dataOut, ~] = CRCdetector(decTbData1);

% 信道估计误差
chl_error=idealhD_grid-chEst;
figure;
surf(abs((chl_error)./(idealhD_grid)).^2)

