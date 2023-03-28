%128 子载波 1.92M采样率 （0.015Mhz*128）

%% Powered  by Boboo 2023.3.9
%修改了MATLAB内置的LTE库函数底层的DSP代码，修改后的代码后缀加了"_128" 
%修改后可以实现128解调  注意：要仔细看TAG调制的函数Modbytag_128
% N=128, len=72 若之后修改子载波数，修改子函数内的N和len应该即可
clc
clear
%% Input an image file and convert to binary stream
load('Picture_all.mat');
index = 1;
fData = Picture_all(index).data;   % Read image data from file
scale = 0.02;                      % Image scaling factor
origSize = size(fData);            % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Resize image
imsize = size(fData);              % Store new image size
binData = dec2bin(fData(:),8);     % Convert to 8 bit unsigned binary
trData = reshape((binData-'0').',1,[]).'; % Create binary stream  %二进制比特

%% 产生OFDM
Global_Parameters;
[eNodeBOutput,txGrid,rmcc,freqGrid] = lteRMCDLTool_fullfreq(rmc3,trData);  %freqGrid是IFFT前的频域所有值
%freqGrid 为什么不只有128个非0？？？？

%% 信道
rxWaveform=eNodeBOutput;

%% 解调LTE链路的OFDM
enb = rmc3; 

frame=txGrid;
grid=txGrid;
waveform=eNodeBOutput;
%SFO
frequencyOffset = lteFrequencyOffset(enb,rxWaveform);%%%
rxWaveform = lteFrequencyCorrect(enb,rxWaveform,frequencyOffset);
fprintf('\nCorrected a frequency offset of %i Hz.\n',frequencyOffset)

[frameOffset,corr] = lteDLFrameOffset(enb,rxWaveform);  %用PSS/SSS进行帧同步

%解调原LTE链路
[rxGrid,freqGrid_r] = lteOFDMDemodulate_fullfreq(enb,rxWaveform);


%% TAG调制
Sig_f=rxWaveform(1:1920);   %过信道 %用了pss的？？？
[Sig_f_tag,p,offset_phi,preamble_p,BIT,Bit01]=Modbytag_128(Sig_f);  %p为 offset_int 
Noise = 10^(-30/20)*rand(1920,1);
Sig_f_tag_and_Noise=Sig_f_tag;
%% 接收机解调TAG

% 4.1 把接收到的经过TAG二次调制的信号重新排列，IFFT到频域栅格  OFDM Rx
[rxGrid_72,rxGrid_128] = lteOFDMDemodulate_fullfreq(enb,Sig_f_tag_and_Noise);

% 4.2 利用OFDM1 估计信道 扩展为14个OFDM符号的信道
H1_est=rxGrid_72(:,1)./txGrid(:,1);%被TAG调制后  等效信道的频域值 （认为其14个OFDM符号内不变）
for i=1:size(H1_est,1)
    if(isfinite(H1_est(i))==0)  
        H1_est(i)=(H1_est(i-1)+H1_est(i+1))/2;
    end
end

H_est=repmat(H1_est,1,14);%被TAG调制后  等效信道的频域值 （认为其14个OFDM符号内不变）


% 4.3 信道均衡(72) 消除了信道的影响 
nVar=1e-5;
yRec_tag = Equalizer(rxGrid_72, H_est, nVar, 1);  %zero-forcing

% 4.4 估计整数倍偏移p，及128子载波信道
 [p_est,H128_est]=findp_128(yRec_tag(:,2),preamble_p,txGrid(:,2),rxGrid_128(:,2)); 
% 4.5 信道均衡（128）  （前两列可以不计算）
yRec_tag128 = Equalizer(rxGrid_128, H128_est, nVar, 1);  %zero-forcing

%4.6 解调TAG的72bit
%现在每次解调72个
%在频域解调  （矩阵求逆）
bit_seq_matrix=[];
for cnt=3:2+12
bit_seq=tagdecode128(p_est, yRec_tag128(:,cnt), txGrid(:,cnt));
bit_seq_matrix=[bit_seq_matrix,bit_seq];

end
[number,ratio,loc] = symerr(bit_seq_matrix,Bit01(:,1:12));
sum(loc);
ratio

%IFFT到时域解调   （矩阵乘法）
bit_seq_matrix2=[];
for cnt=3:2+12
bit_seq2=tagdecode128_ifft(p_est, yRec_tag128(:,cnt), txGrid(:,cnt));
bit_seq_matrix2=[bit_seq_matrix2,bit_seq2];
end
[number2,ratio2,loc2] = symerr(bit_seq_matrix2,Bit01(:,1:12));
sum(loc2);
ratio2