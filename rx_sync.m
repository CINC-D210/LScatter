clc
clear
fc=852.5e6;
NumFramesInBurst = 100;

loc = 0;  %%0:本地  1：实时接收

if loc == 1
    radio = comm.SDRuReceiver('Platform','X310','IPAddress', '192.168.10.2');
    radio.MasterClockRate = 200e6;
    R_bps=2.5e6;%  基带速率
    radio.DecimationFactor =radio.MasterClockRate/R_bps ; %分频偶数最好  
    radio.ChannelMapping = 1;     % Receive signals from both channels
    radio.CenterFrequency =fc ; %455.5e6
    radio.Gain = 0;
    radio.SamplesPerFrame = 19200; % Sampling rate is 1.92 MHz. LTE frames are 10 ms long
    radio.OutputDataType = 'double';
    radio.EnableBurstMode = true;
    radio.NumFramesInBurst = NumFramesInBurst;
    radio.OverrunOutputPort = true;
    radio
    
    % % Capture Signal
    burstCaptures = zeros(19200,NumFramesInBurst);
    len = 0;
    for frame = 1:NumFramesInBurst
        while len == 0
            [data,len,lostSamples] = step(radio);
            burstCaptures(:,frame) = data;
        end
        len = 0;
    end
    release(radio);
    save('rx.mat',"burstCaptures");
    %%usrp接收end
else
    %%%%%%%%本地数据start
    burst = load('rx.mat');
    burstCaptures = burst.burstCaptures;
    %%%%%%%%%本地数据end
end
%% tx lte signal
load('Picture_all.mat');
index = 1;
fData = Picture_all(index).data;      % Read image data from file
scale = 0.02;                      % Image scaling factor
origSize = size(fData);            % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Resize image
imsize = size(fData);              % Store new image size
binData = dec2bin(fData(:),8);     % Convert to 8 bit unsigned binary
trData = reshape((binData-'0').',1,[]).'; % Create binary stream  %二进制比特
Global_Parameters;
[eNodeBOutput,txGrid0,rmcc,freqGrid] = lteRMCDLTool_fullfreq(rmc3,trData);
tagWave = reshape(burstCaptures,[],1);
Global_Parameters;
enb = rmc3;
h128s = zeros(128,4,80);
err=[];
err_72=[];
err_72_over=[];
%% OFFSET
for offset = 11:90  %%80个5ms  11-90
    rxWaveform = tagWave(19200*offset:19200*offset+19200*2);
    rxWaveformCorrected = waveformCorrect(rxWaveform,enb);
    %% 接收机解调TAG
%     tag_signal_1pss = rxWaveformCorrected(961:9600-960);%%%%只取前5m
    tag_signal_1pss = rxWaveformCorrected(9600+961:19200-960);%%%%只取后5ms
    tag_e_rates = zeros(10,4);
    tag_e_rates2 = zeros(10,4);
    tag_e_rates3 = zeros(10,4);
    tag_e_rates4 = zeros(10,4);
    preamble_p = load("tag_code.mat").tag_preamble;
    for repeat=1:72
        if preamble_p(repeat)==0
            break
        end
    end
    code_len = 72/(repeat-1);
    Bit01 = load("tag_code.mat").tag_data;
%     preamble_p = vertcat(zeros(12,1),load("tag_code.mat").tag_preamble,zeros(12,1));
%     Bit01 = vertcat(zeros(12,10),load("tag_code.mat").tag_data,zeros(12,10));
    
    for index=0:3 %%%一个5ms里的，pss后的4次14个ofdm信息，（2+12）
        % 4.1 把接收到的经过TAG二次调制的信号重新排列，IFFT到频域栅格  OFDM Rx
        Sig_f_tag_and_Noise = tag_signal_1pss(1920*index+1:1920*index+1920);
%         txGrid = txGrid0(:,index*14+8:index*14+21);
        txGrid = txGrid0(:,index*14+8+70:index*14+21+70);
        [rxGrid_72,rxGrid_128] = lteOFDMDemodulate_fullfreq(enb,Sig_f_tag_and_Noise);
        
        % 4.2 利用OFDM1 估计信道 扩展为14个OFDM符号的信道
        H1_est=rxGrid_72(:,3)./txGrid(:,3);%被TAG调制后  等效信道的频域值 （认为其14个OFDM符号内不变）
%         plot(abs(H1_est(:,1)))
%         for i=1:size(H1_est,1)
%             if(isfinite(H1_est(i))==0)  
%                 H1_est(i)=(H1_est(i-1)+H1_est(i+1))/2;
%             end
%         end
        H_est=repmat(H1_est,1,14);%被TAG调制后  等效信道的频域值 （认为其14个OFDM符号内不变
        % 4.3 信道均衡(72) 消除了信道的影响 
        nVar=1e-5;
        yRec_tag = Equalizer(rxGrid_72, H_est, nVar, 1);  %zero-forcing %%m???????
%          preamble_p = zeros(72,1);
        % 4.4 估计整数倍偏移p，及128子载波信道
        [p_est,H128_est,idx]=findp_128(yRec_tag(:,4),preamble_p,txGrid(:,4),rxGrid_128(:,4)); 

          ps(index+1,offset)=p_est;
        % 4.5 信道均衡（128）  （前两列可以不计算）
        yRec_tag128 = Equalizer(rxGrid_128, H128_est, nVar, 1);  %zero-forcing
  
        %4.6 解调TAG的72bit
        %现在每次解调72个
        %在频域解调  （矩阵求逆）
%         bit_seq_matrix=[];
%         for cnt=5:14
%               bit_seq=tagdecode128(p_est, yRec_tag128(:,cnt), txGrid(:,cnt));            
% %               bit_seq=mydecodeFreq72(p_est, yRec_tag(:,3), txGrid(:,3),code_len);
%             [number,ratio,loc] = symerr(bit_seq,Bit01(:,cnt-4));
%             tag_e_rates(cnt-4,index+1)=ratio;
%         end
        %FFT 100xsampple 72subcarrier
        bit_seq_matrix2=[];
        for cnt=5:14
            bit_seq2=mydecodefft100(p_est, yRec_tag(:,cnt), txGrid(:,cnt),code_len);
            [number2,ratio2,loc2] = symerr(bit_seq2,Bit01(:,cnt-4));
            tag_e_rates2(cnt-4,index+1)=ratio2;
        end
        %IFFT 72subcarrier
        bit_seq_matrix3=[];
        for cnt=5:14
            bit_seq3=mydecode72(p_est, yRec_tag(:,cnt), txGrid(:,cnt),code_len);
            [number3,ratio3,loc3] = symerr(bit_seq3,Bit01(:,cnt-4));
            tag_e_rates3(cnt-4,index+1)=ratio3;
        end
        %FFT 1xsampple 72subcarrier
        bit_seq_matrix4=[];
        for cnt=5:14
            bit_seq4=mydecodefft1(p_est, yRec_tag(:,cnt), txGrid(:,cnt),code_len);
            [number4,ratio4,loc4] = symerr(bit_seq4,Bit01(:,cnt-4));
            tag_e_rates4(cnt-4,index+1)=ratio4;
        end

        nRepeat=72/code_len;
       for cnt=5:14
            bit_seq5=decode72(p_est, yRec_tag(:,cnt), txGrid(:,cnt),idx,nRepeat);
            Bit_down=downsample(Bit01(:,cnt-4),nRepeat);
            [number5,ratio5,loc5] = symerr(bit_seq5,Bit_down);
            tag_e_rates5(cnt-4,index+1)=ratio5;
       end

        nOvers=10;
       for cnt=5:14
            bit_seq6=decode72_overSample(p_est, yRec_tag(:,cnt), txGrid(:,cnt),idx,nRepeat,nOvers);
            Bit_down=downsample(Bit01(:,cnt-4),nRepeat);
            [number6,ratio6,loc6] = symerr(bit_seq6,Bit_down);
            tag_e_rates6(cnt-4,index+1)=ratio6;
        end
    end
      err=[err mean(tag_e_rates2(:))];
      err_72=[err_72 mean(tag_e_rates5(:))];
      err_72_over=[err_72_over mean(tag_e_rates6(:))];
end
mean(err)
mean(err_72)
mean(err_72_over)