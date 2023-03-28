per=0;
nRepeat = 100; 
vari=0;
ofdm_time = eNodeBOutput(11:10+128);
m1=[];
m2=[];
m3=[];
repeat = 6;
snrs = -45:2:-35;
for snr = snrs
    ber=[];
    ber2=[];
    ber3=[];
    for i=1:50
        % % %%%test1
    %     data_num = 20;
    %     tag_time=randi([0 1],data_num,1);
    %     tag_time(find(tag_time==0)) = -1;
    %     tag_time = vertcat(ones((128-data_num)/2,1),tag_time,ones((128-data_num)/2,1));
        % %%%test
        data_num = 12;
        tag_time=randi([0 1],data_num,1);
        tag_time(find(tag_time==0)) = -1;
        tag_time=repmat(tag_time,1,repeat).';%方波调制信号
        tag_time=tag_time(:);
        tag_time = vertcat(ones((128-data_num*repeat)/2,1),tag_time,ones((128-data_num*repeat)/2,1));
        % tag_time=ones(128,1);
        % tag_time(64)=1;
        % %%%test
        
        ofdm_freq = fft(ofdm_time);
        %理想低通滤波器
        ofdm_freq_n_freq = vertcat(ofdm_freq(1:64),zeros((nRepeat-1)*128,1),ofdm_freq(65:128));
        ofdm_freq_n_time = ifft(ofdm_freq_n_freq);%空中的usrp信号（没过信道）
        
        tag_repeat=repmat(tag_time,1,nRepeat).';%方波调制信号
        tag_repeat2=tag_repeat(:);
        
    %     snr = 35;
        modu_n_time = ofdm_freq_n_time.*tag_repeat2+0.707*0.3e-3*rand(12800,1)*10^(-snr/20)+0.0707i*0.3e-3*rand(12800,1)*10^(-snr/20);
        modu_n_freq = fft(modu_n_time);
        
        freq2_5 = [[0];modu_n_freq(2:37);zeros(55,1);modu_n_freq(end-35:end)];%理论上TAG修正后频域的值（没过信道
        %%% method1
        decod_bit = (ifft(freq2_5))./ofdm_time;
        decod_bit = decod_bit(29:100);
        for index=0:data_num-1
            if(sum(real(decod_bit(index*repeat+1:index*repeat+repeat))))<0
                decod_bit(index*repeat+1:index*repeat+repeat)=-1;
            else
                decod_bit(index*repeat+1:index*repeat+repeat)=1;
            end
        end
        [~,ratio,~] = symerr(decod_bit,tag_time(29:100));
        ber=[ber ratio];

        %%% method2
        freq_72s = zeros(72,2^data_num);
        ZEROS = zeros(28,1);
        for code_index=0:2^data_num-1
            code_bin = dec2bin(code_index);
            code_24 = ones(data_num,1)*-1;
            code_24(end-length(code_bin)+1:end) = 2*(code_bin-'0')-1;
            code_temp=repmat(code_24,1,repeat).';%方波调制信号
            code_72=code_temp(:);
            [~,freq_72s(:,code_index+1)]=multiAndFilter(ofdm_time,[ZEROS; code_72; ZEROS]);
        end
        freq72=[modu_n_freq(end-35:end);modu_n_freq(2:37)];
        freq72_rp=repmat(freq72,1,2^data_num);
        error_vec=abs(freq_72s-freq72_rp);  %这里误差=|实际值-真实值|  
        error_sum=sum(error_vec,1);%将误差累加
        [~,code_est]=min(error_sum); 
        code_bin = dec2bin(code_est-1);
        code_24 = ones(data_num,1)*-1;
        code_24(end-length(code_bin)+1:end) = 2*(code_bin-'0')-1;
        code_temp=repmat(code_24,1,repeat).';%方波调制信号
        code_72=code_temp(:);
        [~,ratio,~] = symerr(code_72,tag_time(29:100));
        ber2=[ber2 ratio];
    %     power_per = sum(abs(freq2_5).*abs(freq2_5))/sum(abs(modu_n_freq).*abs(modu_n_freq));
    %     per=per+power_per;
        
        %%% method3
        freq_72s2 = zeros(72,2^data_num);
        ZEROS = zeros(28,1);
        for code_index2=0:2^data_num-1
            code_bin2 = dec2bin(code_index2);
            code_24 = ones(data_num,1)*-1;
            code_24(end-length(code_bin2)+1:end) = 2*(code_bin2-'0')-1;
            code_temp=repmat(code_24,1,repeat).';%方波调制信号
            code_722=code_temp(:);
            fft_temp=fft(ofdm_time.*[ZEROS; code_722; ZEROS]);
            freq_72s2(:,code_index2+1) = [fft_temp(end-35:end);fft_temp(2:37)];
        end
        freq72=[modu_n_freq(end-35:end);modu_n_freq(2:37)];
        freq72_rp=repmat(freq72,1,2^data_num);
        error_vec=abs(freq_72s2-freq72_rp);  %这里误差=|实际值-真实值|  
        error_sum=sum(error_vec,1);%将误差累加
        [~,code_est2]=min(error_sum); 
        code_bin2 = dec2bin(code_est2-1);
        code_24 = ones(data_num,1)*-1;
        code_24(end-length(code_bin2)+1:end) = 2*(code_bin2-'0')-1;
        code_temp=repmat(code_24,1,repeat).';%方波调制信号
        code_722=code_temp(:);
        [~,ratio,~] = symerr(code_722,tag_time(29:100));
        ber3=[ber3 ratio];
    end
    % per
    m1=[m1 mean(ber)];
    m2=[m2 mean(ber2)];
    m3=[m3 mean(ber3)];
end
figure(1)
hold on
plot(snrs,m1);
plot(snrs,m2);
plot(snrs,m3);
hold off
xlabel('SIR')
ylabel('BER')
legend('ifft','fft100','fft1');
title('12bit*6/128xn')
