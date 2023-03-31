clc;clear all ;close all
% 测试PSS信号的捕捉
%% 各种控制信令的生成
num_subframe=20  %仿真的子帧数

% PDSCH   生成资源块的信息
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
%% LTE的参数设置
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

signal_2frame=[];
signal_Grid=[];

for n_subframe=1:num_subframe
%每次循环生成 1个子帧（）
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
txGrid = REmapper_1Tx(modOut, csr_ref, nS, prmLTEPDSCH,pdcch,pss,sss,bch);  %************把用户数据和CSR映射到资源块上*********
%1.8 OFDM transmitter
txSig = OFDMTx(txGrid, prmLTEPDSCH); %发送的串行OFDM（14个，1RB,1ms）
signal_2frame=[signal_2frame;txSig];
signal_Grid=[signal_Grid,txGrid];
nS = nS + 2; if nS > 19, nS = mod(nS, 20); 
end
end

figure;
plot(abs(txSig))
figure;
plot(abs(signal_2frame))

% %查看几个PSS
% pss_where=find(signal_Grid==0.88);
% num_sss=length(pss_where)/72