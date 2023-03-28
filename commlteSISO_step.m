function [dataIn, dataOut, txSig, rxSig, dataRx, yRec, csr_ref]...
    = commlteSISO_step(nS, snrdB, prmLTEDLSCH, prmLTEPDSCH, prmMdl)
%dataIn原始比特, dataOut, 
%txSig发送的串行数据, rxSig接收的串行数据, 
%dataRx接收到资源格上的用户数据的星座点（未均衡）, yRec（对dataRx进行信道均衡后的值）, 
%csr_ref 产生的参考信号
%
clc
%% 1 TX
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
%% 2 Channel
% SISO Fading channel
[rxFade, chPathG] = MIMOFadingChan(txSig, prmLTEPDSCH, prmMdl);
[idealhD,~]= lteIdChEst(prmLTEPDSCH,  prmMdl, chPathG, nS);  %%*********估计了理想的频域信道 方法？？*********

% Add AWG noise
nVar = 10.^(0.1.*(-snrdB));
rxSig =  AWGNChannel(rxFade, nVar);
%% 3 RX
% 3.1 把接收到的信号重新排列，IFFT到频域栅格  OFDM Rx
rxGrid = OFDMRx(rxSig, prmLTEPDSCH);
% 3.2 按层提取用户数据和CSR updated for numLayers -> numTx
[dataRx, csrRx, idx_data] = REdemapper_1Tx(rxGrid, nS, prmLTEPDSCH);
% 3.3 用CSR进行信道估计  MIMO channel estimation
if prmMdl.chEstOn
    %**********用接收到的CSR，和其原本值估计信道，并插值出完整的视频时频栅格上对应的频域信道**************
    chEst = ChanEstimate_1Tx(prmLTEPDSCH, csrRx,  csr_ref, 'interpolate'); %整个栅格上信道的估计值
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
% chl_error=idealhD_grid-chEst;

%% 
end

% %% 
% figure
% plot(abs(real(rxSig).^2))

