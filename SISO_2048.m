%2022.12.1  1.
% 1. 现在TAG调制后的东西 仅仅过了和直射路径一样的信道h，要改 
% 2. 现在假设接收机完美拿到发送的 时x频X信息    
%    实际在解调后可以逆向一步一步恢复（把接收解调的bit恢复 走一遍发送的流程）
clear all; clc;close all;
clear functions
disp('Simulating the LTE Mode 1: Single Tx and Rx antrenna');
%% LTE的参数设置
% PDSCH   生成资源块的信息
numTx          = 1;    % Number of transmit antennas
numRx          = 1;    % Number of receive antennas
chanBW       = 6;    % Index to chanel bandwidth used [1,....6]
contReg       = 1;    % No. of OFDM symbols dedictaed to control information [1,...,3]
modType      =  1;   % Modulation type [1, 2, 3] for ['QPSK,'16QAM','64QAM']
% DLSCH   Turbo解码的信息
cRate            = 1/3; % Rate matching target coding rate 
maxIter         = 6;     % Maximum number of turbo decoding terations  
fullDecode    = 0;    % Whether "full" or "early stopping" turbo decoding is performed
% Channel model
chanMdl        =  'frequency-selective-low-mobility'; 
% chanMdl        =  'flat-low-mobility'; 
corrLvl           = 'Low'; 
% Simulation parametrs
Eqmode        = 2;      % Type of equalizer used [1,2] for ['ZF', 'MMSE']
chEstOn        = 1;     % Whether channel estimation is done or ideal channel model used
maxNumErrs = 1e7; % Maximum number of errors found before simulation stops
maxNumBits = 1e7;  % Maximum number of bits processed before simulation stops
visualsOn     = 1;      % Whether to visualize channel response and constellations
snrdB            = 50;   % Value of SNR used in this experiment

%tag信道
prmMdl_Tag_h.chanMdl = chanMdl;
prmMdl_Tag_h.corrLevel = 'frequency-selective-low-mobility';
prmMdl_Tag_h.chEstOn = 1;  %是否估计信道

prmMdl_Tag_g.chanMdl = chanMdl;
prmMdl_Tag_g.corrLevel = 'frequency-selective-low-mobility';
prmMdl_Tag_g.chEstOn = 1;  %是否估计信道
%% PDSCH   生成资源块的信息
qpsk_box=[0.707+0.707*1j, 0.707-0.707*1j,-0.707-0.707*1j,-0.707+0.707*1j];
for temp=1:1000
    pdcch(temp)=qpsk_box(randi(4));
end
pss=[ones(1,5).*[0.8+0.6j],PSS_generator(),ones(1,5).*[0.8+0.6j]];
sss=j*pss;
%  pss=ones(1,72).*0.99
% sss=ones(1,72).*0.88;

for temp=1:276
    bch(temp)=qpsk_box(randi(4));%276
end
%% Set simulation parametrs & initialize parameter structures
[prmLTEPDSCH, prmLTEDLSCH, prmMdl] = commlteSISO_initialize( chanBW, contReg,  modType, Eqmode,...
    cRate,maxIter, fullDecode, chanMdl, corrLvl, chEstOn, maxNumErrs, maxNumBits);
% clear chanBW contReg numTx numRx modType Eqmode cRate maxIter fullDecode chanMdl corrLvl chEstOn maxNumErrs maxNumBits;
%至此 程序的变量仅有3个结构体变量： prmLTEPDSCH资源块生成相关  prmLTEDLSCH解码相关   prmMdl信道相关
%%
zReport_data_rate(prmLTEDLSCH, prmLTEPDSCH);  %打印参数信息
hPBer = comm.ErrorRate;  %hPBer统计误码率
% snrdB=prmMdl.snrdBs(end);
maxNumErrs=prmMdl.maxNumErrs;
maxNumBits=prmMdl.maxNumBits;
%% 
nS = 0; % Slot number, one of [0:2:18]
% Measures = zeros(3,1); %initialize BER output

%  1子帧=1*14=14个OFDM符号=(2048*15)*1=30720采样点
%%  1 TX
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
txGrid = REmapper_1Tx(modOut, csr_ref, nS, prmLTEPDSCH,pdcch,pss,sss,bch);  %************把用户数据和CSR映射到资源块上*********
%1.8 OFDM transmitter
txSig = OFDMTx(txGrid, prmLTEPDSCH); %发送的串行OFDM

%1.9 把产生的OFDM符号经过TAG调制（14个OFDM+CP）
[txSig_tag,p,offset_phi,preamble_p,BIT,Bit01]=Modbytag(txSig);  %p为 offset_int 
%%  2.1  Channel-h
% SISO Fading channel
[rxFade, chPathG] = MIMOFadingChan(txSig, prmLTEPDSCH, prmMdl);
% [idealhD,idealhD_grid]= lteIdChEst(prmLTEPDSCH,  prmMdl, chPathG, nS);  %%*********估计了理想的频域信道 方法？？*********

% Add AWG noise
nVar = 10.^(0.1.*(-snrdB));
rxSig =  AWGNChannel(rxFade, nVar);

%%  2.1  Channel-f-tag-g 
% 被TAG调制的信号过一次信道（假设fg中有一段信道是常数（LOS））
[rxFade_tag, chPathG_tag] = MIMOFadingChan(txSig_tag, prmLTEPDSCH, prmMdl);

% [idealhD_tag,idealhD_grid_tag]= lteIdChEst(prmLTEPDSCH,  prmMdl, chPathG, nS);   % 估信道 应该是错的
rxSig_tag =  AWGNChannel(rxFade_tag, nVar);
%% 3 RX-LTE 
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


%% 4 RX-TAG

TxGrid = OFDMRx(txSig, prmLTEPDSCH);  %假设完全恢复了原信号，这里得到

% 4.1 把接收到的经过TAG二次调制的信号重新排列，IFFT到频域栅格  OFDM Rx
[rxGrid_1200,rxGrid_2048] = OFDMRx_2048(rxSig_tag, prmLTEPDSCH);
% 4.2 利用OFDM1 估计信道 扩展为14个OFDM符号的信道
H1_est=rxGrid_1200(:,1)./txGrid(:,1);%被TAG调制后  等效信道的频域值 （认为其14个OFDM符号内不变）
H_est=repmat(H1_est,1,14);
% 4.3 信道均衡(1200) 消除了信道的影响 
yRec_tag = Equalizer(rxGrid_1200, H_est, nVar, prmLTEPDSCH.Eqmode);
% 4.4 估计整数倍偏移p，及2048子载波信道
 [p_est,idx,H2048_est]=findp(yRec_tag(:,2),preamble_p,txGrid(:,2),rxGrid_2048(:,2)); 
% 4.5 信道均衡（2048）  （前两列可以不计算）
yRec_tag2048 = Equalizer(rxGrid_2048, H2048_est, nVar, prmLTEPDSCH.Eqmode);

%需要TAG均衡后的频域值，发送端栅格放置的原本频域值， 前导码（未来要固定）的时域值

%4.6 解调TAG的1200bit
%现在每次解调1200个
bit_seq_matrix=[];

for cnt=3:2+12
bit_seq=tagdecode2048(p_est, yRec_tag2048(:,cnt), txGrid(:,cnt));
bit_seq_matrix=[bit_seq_matrix,bit_seq];
end
% txGrid为LTE解调的频域数据，用于恢复在没有TAG调制情况下的时域数据

% Bit10=BIT;
% Bit10(find(BIT==1))=0;
% Bit10(find(BIT==-1))=1;
[number,ratio,loc] = symerr(bit_seq_matrix,Bit01(:,1:12));
% [number,ratio] = biterr(bit_seq_matrix,Bit10)
% %% 验证估计信道的准确性：
% 'phi:'
% offset_phi
% 'e_jfai'
% exp(j*offset_phi)
% % t=H1_est./(exp(j*offset_phi));   %估计的信道/ 相位的修正
% t=H1_est./idealhD_grid(:,1);    %级联信道/f
% '信道修正均值'
% mean(t)