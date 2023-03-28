% Script for SISO LTE (mode 1)
%
% Single codeword transmission only,
%
clear all; clc;close all;
clear functions
disp('Simulating the LTE Mode 1: Single Tx and Rx antrenna');
%% Set simulation parametrs & initialize parameter structures
commlteSISO_params;  %运行函数  读入仿真参数
[prmLTEPDSCH, prmLTEDLSCH, prmMdl] = commlteSISO_initialize( chanBW, contReg,  modType, Eqmode,...
    cRate,maxIter, fullDecode, chanMdl, corrLvl, chEstOn, maxNumErrs, maxNumBits);
clear chanBW contReg numTx numRx modType Eqmode cRate maxIter fullDecode chanMdl corrLvl chEstOn maxNumErrs maxNumBits;
%至此 程序的变量仅有3个结构体变量： prmLTEPDSCH资源块生成相关  prmLTEDLSCH解码相关   prmMdl信道相关
% 2个常数： snrdB 好像没用到？？   visualsOn：可视化
%%
zReport_data_rate(prmLTEDLSCH, prmLTEPDSCH);  %打印参数信息
hPBer = comm.ErrorRate;   %hPBer统计误码率
% snrdB=prmMdl.snrdBs(end);
maxNumErrs=prmMdl.maxNumErrs;
maxNumBits=prmMdl.maxNumBits;
%% Simulation loop
nS = 0; % Slot number, one of [0:2:18]
Measures = zeros(3,1); %initialize BER output
while (( Measures(2)< maxNumErrs) && (Measures(3) < maxNumBits))
    [dataIn, dataOut, txSig, rxSig, dataRx, yRec, csr] = ...
        commlteSISO_step(nS, snrdB, prmLTEDLSCH, prmLTEPDSCH, prmMdl);
%dataIn原始比特, dataOut, 
%txSig发送的串行数据, rxSig接收的串行数据, 
%dataRx接收到资源格上的的星座点（未均衡）, yRec（对dataRx进行信道均衡后的值）, 
%csr_ref 产生的参考信号

    % Calculate  bit errors
    Measures = step(hPBer, dataIn, dataOut);
    % Visualize constellations and spectrum
    if visualsOn
        zVisualize( prmLTEPDSCH, txSig, rxSig, yRec, dataRx, csr, nS);
    end;
    % Update subframe number 0-18循环 每次+2
    nS = nS + 2; if nS > 19, nS = mod(nS, 20); end;
end
disp(Measures);